!==================================================================================!
!                                      ld                                          !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This module performs the lattice dynamics                                        !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module ld

  use constants
  use rng
  use cell
  use ddp
  use spglib_f08

  implicit none

  private

  public :: ld_calc

  real(kind=gp),     public :: ld_time=0.0_gp
  real(kind=gp),     public :: dm_time=0.0_gp
  real(kind=gp),     public :: cell_time=0.0_gp
  real(kind=gp),     public :: therm_time=0.0_gp
  real(kind=gp),     public :: dsyev_time=0.0_gp
  real(kind=gp),     public :: map_time=0.0_gp
  real(kind=gp),     public :: kpath_time=0.0_gp
  real(kind=gp),     public :: dos_time=0.0_gp
  
  real(kind=gp),     public :: tmax
  integer,           public :: ncells
  integer,           public :: natoms
  integer,           public :: ntemp
  integer,           public :: nsamp
  integer,           public :: nkpts
  character(len=10), public :: unit
  logical,           public :: thermo
  logical,           public :: dos
  logical,           public :: dispersion
  logical,           public :: dryrun=.false.

  !---------------------------------------------------!

  integer, parameter :: unit_disp=33
  integer, parameter :: unit_dos=34

  real(kind=gp) :: kb=8.617333262e-5_gp
  real(kind=gp) :: toev=0.06465415133595565911_gp
  real(kind=gp) :: tounit

contains

  subroutine ld_calc(str)

    type(structure), intent(in) :: str

    !-------------------------------------------------------------!

    type(structure) :: super

    real(kind=gp),    allocatable, dimension(:,:) :: dmat

    real(kind=gp) :: cst,start_time,end_time

    integer :: i,j,k,mat(3,3)

    select case(unit)
    case('meV')
       tounit=1000.0_gp
    case('cm-1')
       tounit=8065.541154_gp
    case('THz')
        tounit=241.79884076620228_gp
    case default
       stop 'ld_calc: unit not recognised - allowed: meV, cm-1, THz'
    end select
    
    call cpu_time(start_time)

    if(all(str%supercell_matrix.eq.0)) then

       if(ncells.lt.0) ncells=(natoms-1)/str%num_ions+1

       call shake_nearly_cubic_supercell(str,ncells,mat,cst)

       write (stderr,'(a,f8.4)') ' deviation from cubic: ',cst

    else

       mat=str%supercell_matrix
       
    end if

    ncells=0
    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       ncells=ncells+mat(1,i)*(mat(2,j)*mat(3,k)-mat(2,k)*mat(3,j))
    end do

    super=supercell(mat,str)
    call write_cell(super)

    write (stderr,'(a,i0,a,i0,a)') ' supercell matrix generating ',ncells,' cells containing ',super%num_ions,' atoms:'
    write (stderr,*)
    write (stderr,'(3i5)') mat
    write (stderr,*)
    write (stderr,'(1x,3f8.3)') super%lattice_abc(1:3)
    write (stderr,'(1x,3f8.2)') super%lattice_abc(4:6)
    write (stderr,*)

    call cpu_time(end_time)
    
    cell_time=cell_time+end_time-start_time

    if(dryrun) return

    allocate(dmat(3*str%num_ions,3*super%num_ions))

    call get_dynamical_matrix(super,dmat)
    
    if(thermo) call calc_thermodynamics(str,super,dmat)
    
    if(dispersion) call calc_dispersion(str,super,dmat)

    if(dos) call calc_dos(str,super,dmat)

    call cpu_time(end_time)

    ld_time=ld_time+end_time-start_time

    return

  end subroutine ld_calc

  subroutine get_dynamical_matrix(str,dmat)

    type(structure),               intent(inout) :: str
    real(kind=gp), dimension(:,:), intent(out)   :: dmat

    real(kind=gp) :: e,g(3,str%num_ions),gv0(3*str%num_ions),gvp(3*str%num_ions)
    real(kind=gp) :: gvm(3*str%num_ions),pos0(3,str%num_ions),massvec(3*str%num_ions)
    real(kind=gp) :: step=1e-4_gp,start_time,end_time

    integer :: ni,i,n

    call cpu_time(start_time)
    
    if(str%num_symm.gt.1) stop 'get_dynamical_matrix: num_symm>1'

    n=0
    do ni=1,str%num_ions
       do i=1,3
          n=n+1
          massvec(n)=str%ion_mass(ni)
       end do
    end do
    
    pos0=str%ion_positions

    call eval_ddp(str,e,g)

    gv0=reshape(g,(/3*str%num_ions/))

    dmat=0.0_gp
    
    n=0
    do ni=1,str%num_ions/str%num_cells
       
       do i=1,3
          n=n+1

          ! Plus

          str%ion_positions=pos0
          str%ion_positions(i,ni)=str%ion_positions(i,ni)+step

          call eval_ddp(str,e,g)

          gvp=reshape(g,(/3*str%num_ions/))

          ! Minus

          str%ion_positions=pos0
          str%ion_positions(i,ni)=str%ion_positions(i,ni)-step

          call eval_ddp(str,e,g)

          gvm=reshape(g,(/3*str%num_ions/))

          ! Hessian

          dmat(n,:)=-(gvp-gvm)/2/step/sqrt(str%ion_mass(ni))/sqrt(massvec)

       end do
    end do
        
    str%ion_positions=pos0

    call cpu_time(end_time)

    dm_time=dm_time+end_time-start_time
    
  end subroutine get_dynamical_matrix

  subroutine get_dynamical_matrix_k(str,super,kvec,dmat,dmk)

    type(structure),                   intent(in) :: str
    type(structure),                   intent(in) :: super
    real(kind=gp),    dimension(:),    intent(in) :: kvec
    real(kind=gp),    dimension(:,:),  intent(in) :: dmat
    complex(kind=gp), dimension(:,:), intent(out) :: dmk

    integer :: i,j,k,ii,jj,kk,num_phase,nc,nna,nnb,nnc

    real(kind=gp) :: qvec(3),xvec(3),h(3),rsphere

    complex(kind=gp) :: phase

    rsphere=max(rcut,minval(super%lattice_abc(1:3))/2)+0.01_gp
    
    ! Determine the supercell required

    h=car2habc(super%lattice_car)

    nna = int(2*rsphere/h(1))+1 !! Do we need the +1
    nnb = int(2*rsphere/h(2))+1
    nnc = int(2*rsphere/h(3))+1
    
    qvec=matmul(transpose(str%lattice_rec),kvec)*tpi
    
    dmk=(0.0_gp,0.0_gp)
    
    do i=1,3*str%num_ions

       k=0
       do nc=1,super%num_cells

          phase=0.0_gp
          num_phase=0

          do ii=-nna,nna
             do jj=-nnb,nnb
                do kk=-nnc,nnc
                
                   xvec=matmul(super%lattice_car,(/ii,jj,kk/))+super%subcell_vec(:,nc)-super%subcell_vec(:,1)

                   if(norm2(xvec).lt.rsphere) then
                      
                      num_phase=num_phase+1
                      phase=phase+exp(-cmplx(0.0_gp,1.0_gp,gp)*dot_product(qvec,xvec))
                      
                   end if

                end do
             end do
          end do
          
          do j=1,3*str%num_ions

             k=k+1

             if(num_phase.gt.0) dmk(i,j)=dmk(i,j)+dmat(i,k)*phase/real(num_phase,gp)

          end do

       end do
       
    end do
    
  end subroutine get_dynamical_matrix_k

  subroutine calc_thermodynamics(str,super,dmat)

    type(structure),                   intent(in) :: str
    type(structure),                   intent(in) :: super
    real(kind=gp),    dimension(:,:),  intent(in) :: dmat


    real(kind=gp),    allocatable, dimension(:,:)   :: amat
    real(kind=gp),    allocatable, dimension(:)     :: eval,work,hw
    real(kind=gp)                                   :: temp
    real(kind=gp)                                   :: e,fmu,g(3,str%num_ions)

    integer, allocatable, dimension(:,:) :: isc 

    integer :: lwork,info,n,nt,nc,ncc,ni,nii,nccc

    integer :: npts=2001

    integer :: ii,jj,kk,lvec(3,27),ivec(3),incncc(3),incnccnccc(3)

    real(kind=gp) :: start_time,end_time,extra_time,emax,sigma

    real(kind=gp), allocatable, dimension(:) :: egrid,phdos

    call cpu_time(start_time)

    lwork=10*3*super%num_ions

    allocate(amat(3*super%num_ions,3*super%num_ions),eval(3*super%num_ions),&
         hw(3*super%num_ions),work(lwork),isc(3,super%num_cells),egrid(npts),phdos(npts))

    n=0
    do ii=-1,1
       do jj=-1,1
          do kk=-1,1
             n=n+1

             lvec(:,n) = nint(matmul(str%lattice_rec,matmul(super%lattice_car,(/ii,jj,kk/))))

          end do
       end do
    end do

    do nc=1,super%num_cells
       isc(:,nc)=nint(matmul(str%lattice_rec,super%subcell_vec(:,nc)))
    end do

    !$omp parallel do private(ivec,ncc,n,nccc,ni,nii,incncc,incnccnccc) schedule(dynamic)
    do nc=1,super%num_cells
       do ncc=nc,super%num_cells
          incncc(:)=isc(:,nc)-isc(:,ncc)

          loop: do nccc=1,super%num_cells

             incnccnccc(:)=isc(:,nccc)+incncc(:)

             do n=1,27
                ivec(:)=lvec(:,n)+incnccnccc(:)
                if(all(ivec.eq.0)) exit loop
             end do

          end do loop

          do ni=1,3*str%num_ions
             do nii=1,3*str%num_ions

                amat((ncc-1)*3*str%num_ions+nii,(nc-1)*3*str%num_ions+ni)=dmat(ni,(nccc-1)*3*str%num_ions+nii)

             end do
          end do

       end do
    end do
    !$omp end parallel do

    call cpu_time(extra_time)

    map_time=map_time+extra_time-start_time
    
    if(gp.eq.sp) then
       call ssyev('N','L',3*super%num_ions,amat,3*super%num_ions,eval,work,lwork,info)
    else
       call dsyev('N','L',3*super%num_ions,amat,3*super%num_ions,eval,work,lwork,info)
    end if
    
    call cpu_time(end_time)

    dsyev_time=dsyev_time+end_time-extra_time

    write (stderr,'(a,i0,a,f0.5,a)') ' negative eigenvalues: ',count(eval<0.0_gp),' minimum: ',minval(eval)*tounit,' '//trim(unit)
    write (stderr,*)

    hw=sqrt(max(0.0_gp,eval))*toev

    call eval_ddp(str,e,g)

    write(stderr,'(a,f0.5,a)') ' zero point energy: ',sum(hw)/super%num_ions*str%num_ions/2*tounit,' '//trim(unit)
    write(stderr,*)
    write (stderr,*) '                                 U/eV'
    
    write (stderr,'(15x,f25.7)') e
    write (stderr,*)
    write (stderr,*) '          T/K                  F-U/eV'
    do nt=1,ntemp
       if(ntemp.gt.1) then
          temp=real(nt-1,gp)/real(ntemp-1,gp)*tmax
       else
          temp=0.0_gp
       end if

       fmu=0.0_gp
       do n=4,size(eval)
          fmu=fmu+hw(n)/2
          if(temp.gt.0.0_gp) fmu=fmu+kb*temp*log(1-exp(-hw(n)/(kb*temp)))
       end do
       write (stdout,'(f15.3,f25.7)') temp,fmu*str%num_ions/super%num_ions
    end do

    call cpu_time(end_time)

    therm_time=therm_time+end_time-start_time

    if(.not.dos) then

       phdos=0.0_gp

       emax=maxval(hw)*tounit*1.1_gp

       emax=(int(emax/5)+1)*5

       sigma=emax/real(npts-1,gp)*4.0_gp

       do n=1,npts
          egrid(n)=real(n-1,gp)/real(npts-1,gp)*emax
       end do

       do n=1,size(hw)
          call add_gaussian(hw(n)*tounit,egrid,phdos,sigma,1.0_gp/real(ncells,gp))
       end do

       call write_dos(egrid,phdos)

    end if

  end subroutine calc_thermodynamics
  
  subroutine calc_dispersion(str,super,dmat)

    type(structure),                   intent(in) :: str
    type(structure),                   intent(in) :: super
    real(kind=gp),    dimension(:,:),  intent(in) :: dmat
    
    real(kind=gp) :: kvec(3),start_time,end_time

    complex(kind=gp), allocatable, dimension(:)     :: cwork
    complex(kind=gp), allocatable, dimension(:,:)   :: dmk
    real(kind=gp),    allocatable, dimension(:)     :: rwork
    real(kind=gp),    allocatable, dimension(:,:,:) :: ek
    
    integer :: np,nk,lwork,info

    if(str%num_path.eq.0) stop 'calc_dispersion: set PHONON_FINE_KPOINT_PATH in <seed>.cell'
    
    call cpu_time(start_time)

    lwork=10*3*str%num_ions
    allocate(dmk(3*str%num_ions,3*str%num_ions),ek(3*str%num_ions,nkpts,str%num_path-1),cwork(lwork),rwork(2*lwork))
    
    do np=1,str%num_path-1

       !$omp parallel do private(kvec,cwork,rwork,dmk) schedule(dynamic)
       do nk=1,nkpts

          kvec(:)=str%kpoint_path(:,np)+real(nk-1,gp)/real(nkpts-1,gp)*(str%kpoint_path(:,np+1)-str%kpoint_path(:,np))
          
          call get_dynamical_matrix_k(str,super,kvec,dmat,dmk)
          
          if(gp.eq.sp) then
             call cheev('N','U',3*str%num_ions,dmk,3*str%num_ions,ek(:,nk,np),cwork,lwork,rwork,info)
          else
             call zheev('N','U',3*str%num_ions,dmk,3*str%num_ions,ek(:,nk,np),cwork,lwork,rwork,info)
          end if
          
          where (ek(:,nk,np).ge.0.0_gp)
             ek(:,nk,np)=sqrt(ek(:,nk,np))*toev
          elsewhere
             ek(:,nk,np)=-sqrt(abs(ek(:,nk,np)))*toev
          end where
          
       end do
       !$omp end parallel do
       
    end do

    call cpu_time(end_time)

    kpath_time=kpath_time+end_time-start_time

    call write_dispersion(str,ek)
    
  end subroutine calc_dispersion

  subroutine calc_dos(str,super,dmat)

    type(structure),                   intent(in) :: str
    type(structure),                   intent(in) :: super
    real(kind=gp),    dimension(:,:),  intent(in) :: dmat

    real(kind=gp) :: kvec(3),start_time,end_time

    complex(kind=gp), allocatable, dimension(:)   :: cwork
    complex(kind=gp), allocatable, dimension(:,:) :: dmk
    real(kind=gp),    allocatable, dimension(:)   :: rwork,e,egrid,phdos

    real(kind=gp) :: emax,sigma,zpe

    integer :: npts=2001,nsmall=100

    integer :: lwork,info,i,n

    call cpu_time(start_time)

    lwork=10*3*str%num_ions
    allocate(e(3*str%num_ions),dmk(3*str%num_ions,3*str%num_ions),cwork(lwork),rwork(2*lwork),egrid(npts),phdos(npts))

    emax=0.0_gp
    !$omp parallel do private(kvec,cwork,rwork,dmk,e) schedule(dynamic)
    do i=1,nsmall

       kvec(:)=random_triple()

       call get_dynamical_matrix_k(str,super,kvec,dmat,dmk)

       if(gp.eq.sp) then
          call cheev('N','U',3*str%num_ions,dmk,3*str%num_ions,e,cwork,lwork,rwork,info)          
       else
          call zheev('N','U',3*str%num_ions,dmk,3*str%num_ions,e,cwork,lwork,rwork,info)
       end if
       
       if(maxval(e).gt.emax) then
          emax=maxval(e)
       end if

    end do
    !$omp end parallel do

    emax=sqrt(abs(emax))*toev*tounit*1.1_gp

    emax=(int(emax/5)+1)*5

    sigma=emax/real(npts-1,gp)*4.0_gp
    
    do n=1,npts
       egrid(n)=real(n-1,gp)/real(npts-1,gp)*emax
    end do

    phdos=0.0_gp
    zpe=0.0_gp
    !$omp parallel do private(kvec,cwork,rwork,dmk,e) reduction(+:phdos,zpe) schedule(dynamic)
    do i=1,nsamp

       kvec(:)=random_triple()

       call get_dynamical_matrix_k(str,super,kvec,dmat,dmk)

       if(gp.eq.sp) then
          call cheev('N','U',3*str%num_ions,dmk,3*str%num_ions,e,cwork,lwork,rwork,info)          
       else
          call zheev('N','U',3*str%num_ions,dmk,3*str%num_ions,e,cwork,lwork,rwork,info)
       end if
       
       where (e.ge.0.0_gp)
          e=sqrt(e)*toev
       elsewhere
          e=-sqrt(abs(e))*toev
       end where

       do n=1,size(e)
          call add_gaussian(e(n)*tounit,egrid,phdos,sigma,1.0_gp/real(nsamp,gp))
       end do

       zpe=zpe+sum(e)/2/real(nsamp,gp)
       
    end do
    !$omp end parallel do

    write (stderr,*)
    write (stderr,'(a,f0.3,a)') ' zero point energy (dos): ',zpe*tounit,' '//trim(unit)
    
    call write_dos(egrid,phdos)

    call cpu_time(end_time)

    dos_time=dos_time+end_time-start_time

  end subroutine calc_dos
  
  subroutine write_dispersion(str,ek)
    
    type(structure),                 intent(in) :: str
    real(kind=gp), dimension(:,:,:), intent(in) :: ek

    integer :: i,np,nk,nb,stat

    real(kind=gp) :: q,dq,qmax,tick,emax
    
    character(len=240) :: ctemp2
    
    open(unit=unit_disp,file=trim(seedname)//'-disp.agr',form='formatted',status='unknown',err=99)

    emax=maxval(ek)*tounit
    
    emax=int(emax)+1
    
    tick=int(emax/10.0_gp)+1

    qmax=0
    do np=1,str%num_path-1
       
       dq=norm2(real(1,gp)/real(nkpts-1,gp)*(matmul(transpose(str%lattice_rec),str%kpoint_path(:,np+1)-str%kpoint_path(:,np))))
       qmax=qmax+dq*(nkpts-1)
       
    end do
    
    call agr_head_disp(0.0_gp,qmax,0.0_gp,emax,tick,tick)

    do i=1,str%num_path-1
       write (unit_disp,'(a,i0,a)') '@    s',i-1,' line color 3'
       write (unit_disp,'(a,i0,a)') '@    s',i-1,' line linewidth 1.0'
    end do

    do i=str%num_path,str%num_path+size(ek,1)
       write (unit_disp,'(a,i0,a)') '@    s',i-1,' line color 2'
       write (unit_disp,'(a,i0,a)') '@    s',i-1,' line linewidth 2.0'
    end do

    write (unit_disp,'(a)') '@with string'
    write (unit_disp,'(a)') '@    string on'
    write (unit_disp,'(a)') '@    string loctype view'
    write (unit_disp,'(a,f10.5,a,f10.5)') '@    string ',0.11_gp,',',0.0125_gp
    write (unit_disp,'(a)') '@    string color 1'
    write (unit_disp,'(a)') '@    string rot 0'
    write (unit_disp,'(a)') '@    string font 4'
    write (unit_disp,'(a)') '@    string just 0'
    write (unit_disp,'(a)') '@    string char size 1.20000'
    if(str%path_label(1).eq.'G') then
       write (unit_disp,'(a)') '@    string def "\x G \f{}"'
    else
       write (unit_disp,'(a)') '@    string def "'//str%path_label(1)//'"'
    end if
 
    q=0
    do np=1,str%num_path-1
       
       dq=norm2(real(1,gp)/real(nkpts-1,gp)*(matmul(transpose(str%lattice_rec),str%kpoint_path(:,np+1)-str%kpoint_path(:,np))))
       q=q+dq*(nkpts-1)
      
       write (unit_disp,'(a)') '@with string'
       write (unit_disp,'(a)') '@    string on'
       write (unit_disp,'(a)') '@    string loctype view'
       write (unit_disp,'(a,f10.5,a,f10.5)') '@    string ',0.11+q/qmax*1.1275_gp,',',0.0125_gp
       write (unit_disp,'(a)') '@    string color 1'
       write (unit_disp,'(a)') '@    string rot 0'
       write (unit_disp,'(a)') '@    string font 4'
       write (unit_disp,'(a)') '@    string just 0'
       write (unit_disp,'(a)') '@    string char size 1.20000'
       if(str%path_label(np+1).eq.'G') then
          write (unit_disp,'(a)') '@    string def "\x G \f{}"'
       else
          write (unit_disp,'(a)') '@    string def "'//str%path_label(np+1)//'"'
       end if
       
    end do

    q=0
    do np=1,str%num_path-1
       
       dq=norm2(real(1,gp)/real(nkpts-1,gp)*(matmul(transpose(str%lattice_rec),str%kpoint_path(:,np+1)-str%kpoint_path(:,np))))
       q=q+dq*(nkpts-1)
       write (unit_disp,*) q,0.0_gp
       write (unit_disp,*) q,emax
       write (unit_disp,*) '&'
    end do
    
    do nb=1,size(ek,1)

       q=0.0_gp
       do np=1,str%num_path-1

          dq=norm2(real(1,gp)/real(nkpts-1,gp)*(matmul(transpose(str%lattice_rec),str%kpoint_path(:,np+1)-str%kpoint_path(:,np))))

          
          do nk=1,nkpts

             write (unit_disp,*) q,ek(nb,nk,np)*tounit
             
             q=q+dq

          end do
          
          q=q-dq

       end do
       
       write (unit_disp,*) '&'

    end do

    close(unit_disp)
    
    ctemp2="gracebat "//trim(seedname)//"-disp.agr -hdevice PDF -hardcopy -printfile "//trim(seedname)//"-disp.pdf"
    stat=0
    call system(ctemp2,stat)
    ctemp2='Problem executing external command :: '//trim(ctemp2)
    if (stat.ne.0) then
       write (stderr,'(a)') trim(ctemp2)
    end if
    
    return

99  stop 'write_dispersion: problem opening file'

  end subroutine write_dispersion

  subroutine write_dos(egrid,phdos)

    real(kind=gp), dimension(:), intent(in) :: egrid
    real(kind=gp), dimension(:), intent(in) :: phdos

    integer :: n,stat

    real(kind=gp) :: tick

    character(len=240) :: ctemp2

    open(unit=unit_dos,file=trim(seedname)//'-dos.agr',form='formatted',status='unknown',err=99)

    tick=maxval(egrid)/5.0_gp

    call agr_head_dos(0.0_gp,maxval(egrid),0.0_gp,maxval(phdos)*1.1_gp,tick,tick)

    do n=1,size(egrid)
       write(unit_dos,*) egrid(n),phdos(n)
    end do
    write (unit_dos,*) '&'

    close(unit_disp)

    ctemp2="gracebat "//trim(seedname)//"-dos.agr -hdevice PDF -hardcopy -printfile "//trim(seedname)//"-dos.pdf"
    stat=0
    call system(ctemp2,stat)
    ctemp2='Problem executing external command :: '//trim(ctemp2)
    if (stat.ne.0) then
       write (stderr,'(a)') trim(ctemp2)
    end if
    
    return

99  stop 'write_dos: problem opening file'

  end subroutine write_dos
  
  subroutine agr_head_disp(xmin,xmax,ymin,ymax,xtick,ytick)

    real(kind=gp), intent(in) :: xmin,xmax,ymin,ymax,xtick,ytick

    write(unit_disp,'(a)') '@version 50109'
    write(unit_disp,'(a)') '@default linewidth 2.0'
    write(unit_disp,'(a)') '@default linestyle 1'
    write(unit_disp,'(a)') '@g0 on'
    write(unit_disp,'(a)') '@with g0'
    write(unit_disp,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_disp,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_disp,'(a)') '@map color 0 to (255, 255, 255), "white"'
    write(unit_disp,'(a)') '@map color 1 to (0, 0, 0), "black"'
    write(unit_disp,'(a)') '@map color 2 to (228, 26, 28), "red"'
    write(unit_disp,'(a)') '@map color 3 to (55, 126, 184), "blue"'
    write(unit_disp,'(a)') '@map color 4 to (77, 175, 74), "green"'
    write(unit_disp,'(a)') '@map color 5 to (152, 78, 163), "purple"'
    write(unit_disp,'(a)') '@map color 6 to (255, 127, 0), "orange"'
    write(unit_disp,'(a)') '@map color 7 to (255, 255, 51), "yellow"'
    write(unit_disp,'(a)') '@map color 8 to (166, 86, 40), "brown"'
    write(unit_disp,'(a)') '@map color 9 to (247, 129, 191), "pink"'
    write(unit_disp,'(a)') '@map color 10 to (153, 153, 153), "grey"'
    write(unit_disp,'(a)') '@map color 11 to (166, 206, 227), "lightblue"'
    write(unit_disp,'(a)') '@map color 12 to (178, 223, 138), "lightgreen"'
    write(unit_disp,'(a)') '@map color 13 to (251, 154, 153), "lightred"'
    write(unit_disp,'(a)') '@map color 14 to (253, 191, 111), "lightorange"'
    write(unit_disp,'(a)') '@map color 15 to (202, 178, 214), "lightpurple"'
    write(unit_disp,'(a,f12.3)') '@    world xmin',xmin
    write(unit_disp,'(a,f12.3)') '@    world xmax',xmax
    write(unit_disp,'(a,f12.3)') '@    world ymin',ymin
    write(unit_disp,'(a,f12.3)') '@    world ymax',ymax
    write(unit_disp,'(a)') '@    view xmin 0.12500000'
    write(unit_disp,'(a)') '@    view xmax 1.2500000'
    write(unit_disp,'(a)') '@    view ymin 0.0500000'
    write(unit_disp,'(a)') '@    view ymax 0.9500000'
    write(unit_disp,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_disp,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_disp,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_disp,'(a,f10.2)') '@    xaxis  tick major',xtick
!!$    write(unit_disp,'(a)') '@    xaxis  label "Kpoint"'
!!$    write(unit_disp,'(a)') '@    xaxis  label font 4'
    write(unit_disp,'(a)') '@    xaxis  tick off'
    write(unit_disp,'(a)') '@    xaxis  ticklabel off'
    write(unit_disp,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_disp,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_disp,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_disp,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_disp,'(a)') '@    yaxis  ticklabel prec 0'
    write(unit_disp,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_disp,'(a)') '@    yaxis  label "Energy ('//trim(unit)//')"'
    write(unit_disp,'(a)') '@    yaxis  label font 4'
    write(unit_disp,'(a)') '@    yaxis  ticklabel font 4'

  end subroutine agr_head_disp

  subroutine agr_head_dos(xmin,xmax,ymin,ymax,xtick,ytick)

    real(kind=gp), intent(in) :: xmin,xmax,ymin,ymax,xtick,ytick

    write(unit_dos,'(a)') '@version 50109'
    write(unit_dos,'(a)') '@default linewidth 2.0'
    write(unit_dos,'(a)') '@default linestyle 1'
    write(unit_dos,'(a)') '@g0 on'
    write(unit_dos,'(a)') '@with g0'
    write(unit_dos,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_dos,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_dos,'(a)') '@map color 0 to (255, 255, 255), "white"'
    write(unit_dos,'(a)') '@map color 1 to (0, 0, 0), "black"'
    write(unit_dos,'(a)') '@map color 2 to (228, 26, 28), "red"'
    write(unit_dos,'(a)') '@map color 3 to (55, 126, 184), "blue"'
    write(unit_dos,'(a)') '@map color 4 to (77, 175, 74), "green"'
    write(unit_dos,'(a)') '@map color 5 to (152, 78, 163), "purple"'
    write(unit_dos,'(a)') '@map color 6 to (255, 127, 0), "orange"'
    write(unit_dos,'(a)') '@map color 7 to (255, 255, 51), "yellow"'
    write(unit_dos,'(a)') '@map color 8 to (166, 86, 40), "brown"'
    write(unit_dos,'(a)') '@map color 9 to (247, 129, 191), "pink"'
    write(unit_dos,'(a)') '@map color 10 to (153, 153, 153), "grey"'
    write(unit_dos,'(a)') '@map color 11 to (166, 206, 227), "lightblue"'
    write(unit_dos,'(a)') '@map color 12 to (178, 223, 138), "lightgreen"'
    write(unit_dos,'(a)') '@map color 13 to (251, 154, 153), "lightred"'
    write(unit_dos,'(a)') '@map color 14 to (253, 191, 111), "lightorange"'
    write(unit_dos,'(a)') '@map color 15 to (202, 178, 214), "lightpurple"'
    write(unit_dos,'(a,f12.3)') '@    world xmin',xmin
    write(unit_dos,'(a,f12.3)') '@    world xmax',xmax
    write(unit_dos,'(a,f12.3)') '@    world ymin',ymin
    write(unit_dos,'(a,f12.3)') '@    world ymax',ymax
    write(unit_dos,'(a)') '@    view xmin 0.12500000'
    write(unit_dos,'(a)') '@    view xmax 1.2500000'
    write(unit_dos,'(a)') '@    view ymin 0.1000000'
    write(unit_dos,'(a)') '@    view ymax 0.9500000'
    write(unit_dos,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_dos,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_dos,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_dos,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_dos,'(a)') '@    xaxis  label "Energy ('//trim(unit)//')"'
    write(unit_dos,'(a)') '@    xaxis  label font 4'
    write(unit_dos,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_dos,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_dos,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_dos,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_dos,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_dos,'(a)') '@    yaxis  ticklabel prec 0'
    write(unit_dos,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_dos,'(a)') '@    yaxis  label "Density of States"'
    write(unit_dos,'(a)') '@    yaxis  label font 4'
    write(unit_dos,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_dos,'(a)') '@    yaxis  tick off'
    write(unit_dos,'(a)') '@    yaxis  ticklabel off'
    write(unit_dos,'(a)') '@    yaxis  label place spec'
    write(unit_dos,'(a)') '@    yaxis  label place 0.000000, 0.050000'
    write(unit_dos,'(a)') '@    s0 line type 1'
    write(unit_dos,'(a)') '@    s0 line linestyle 1'
    write(unit_dos,'(a)') '@    s0 line linewidth 2.0'
    write(unit_dos,'(a)') '@    s0 line color 4'
    write(unit_dos,'(a)') '@    s0 line pattern 1'
    write(unit_dos,'(a)') '@    s0 baseline type 0'
    write(unit_dos,'(a)') '@    s0 baseline off'
    write(unit_dos,'(a)') '@    s0 dropline off'
    write(unit_dos,'(a)') '@    s0 fill type 2'
    write(unit_dos,'(a)') '@    s0 fill rule 0'
    write(unit_dos,'(a)') '@    s0 fill color 12'
    write(unit_dos,'(a)') '@    s0 fill pattern 1'
    
  end subroutine agr_head_dos

  subroutine add_gaussian(x0,x,f,sig,wgt)

    real(kind=gp),               intent(in)    :: x0
    real(kind=gp), dimension(:), intent(in)    :: x
    real(kind=gp), dimension(:), intent(inout) :: f
    real(kind=gp),               intent(in)    :: sig
    real(kind=gp),               intent(in)    :: wgt

    integer :: i,n

    n=size(x)

    do i=1,n
       f(i)=f(i)+exp(-(x(i)-x0)**2/2.0_gp/sig**2)/sig*wgt/sqrt(tpi)
    end do

  end subroutine add_gaussian
  
end module ld
