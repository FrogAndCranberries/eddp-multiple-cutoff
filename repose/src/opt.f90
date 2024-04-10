!==================================================================================!
!                                      opt                                         !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This module performs the geometry optimisation                                   !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module opt

  use constants
  use rng
  use cell
  use ddp

  implicit none

  private

  public :: opt_tpsd
  public :: opt_newton
  public :: opt_test

  logical,       public :: test=.false.
  logical,       public :: optimise=.true.
  logical,       public :: converged
  real(kind=gp), public :: gamma=1.0_gp
  integer,       public :: unit_track=21

  real(kind=gp), public :: prec_time=0.0_gp
  
  !---------------------------------------------------!

  character(len=80)     :: trackfile,xyzefile

  real(kind=gp)         :: h1
  
contains

  subroutine opt_tpsd(str,steps,tol,cell_fix)

    !-------------------------------------------------------------!
    ! Two-point Step Size Gradient Methods - Barzilai and Borwein !
    ! IMA Journal of Numerical Analysis (1988) 8, 141-148         !
    !-------------------------------------------------------------!

    type(structure), intent(inout) :: str
    integer,         intent(inout) :: steps
    real(kind=gp),      intent(in) :: tol
    logical,            intent(in) :: cell_fix

    real(kind=gp) :: step,e,p,sig(3,3),s(3,3),lat0(3,3),strain(3,3),derun,de,strun,var

    real(kind=gp), save, allocatable :: g(:,:),x(:),xp(:),x0(:),xp0(:),pxp(:),pxp0(:)

    real(kind=gp) :: start_time,end_time

    integer :: maxsteps

    maxsteps=steps

    if(.not.allocated(g)) allocate(g(3,str%num_ions),x(3*str%num_ions+3*3),xp(3*str%num_ions+3*3),x0(3*str%num_ions+3*3),&
         xp0(3*str%num_ions+3*3),pxp(3*str%num_ions+3*3),pxp0(3*str%num_ions+3*3))

    ! ** Open a file for the convergence data

    if(track.and.optimise) then

       write(trackfile,'(a)') trim(seedname)//".track"       
       open(unit=unit_track,file=trackfile,status="unknown",err=998)

       if(.not.quiet) then
          write(xyzefile,'(a)') trim(seedname)//".xyze"       
          open(unit=unit_xyze,file=xyzefile,status="unknown",err=998)
          call write_xyze(str)
       end if
       
    end if 
    
    p=(str%external_pressure(1,1)+str%external_pressure(2,2)+str%external_pressure(3,3))/3.0_gp

    derun     = 0.0_gp
    strun     = 0.0_gp
    step      = 1e-4_gp
    converged = .true.

    lat0=str%lattice_car

    strain=ident
    
    x=[reshape(str%ion_positions,(/3*str%num_ions/)),reshape(strain,(/3*3/))]
           
    steps = 0
    do while(((derun.gt.log10(tol)).or.(steps.lt.2)))
       steps=steps+1
       total_steps=total_steps+1

       if(steps.gt.maxsteps) then 
          converged=.false.
          exit
       end if
              
       ! * Compute gradient
       
       if(.not.cell_fix) then
          call eval_ddp(str,e,g,sig)
          s = sig-str%external_pressure*str%volume
       else
          call eval_ddp(str,e,g)
          s = 0.0_gp
       end if
       
       ! * Convert gradient to vector. Eqn 17 of J. Chem. Phys. 144, 164109 (2016)

       xp=[reshape(matmul(strain,g),(/3*str%num_ions/)),reshape(matmul(s,transpose(inv(strain))),(/3*3/))]

       ! * Store initial value of enthalpy, initialise derun
       
       if(total_steps.eq.1) then
          h1=e+p*str%volume
          derun=0.0_gp
       end if
       
       if(.not.optimise) exit

       ! * Precondition gradient
       
       call cpu_time(start_time)
      
       call precondition(str,xp(1:3*str%num_ions),pxp(1:3*str%num_ions),-2.0_gp) ! ** Get L from ddp

       pxp(3*str%num_ions+1:3*str%num_ions+3*3)=xp(3*str%num_ions+1:3*str%num_ions+3*3)/gamma/str%num_ions
       
       call cpu_time(end_time)

       prec_time=prec_time+end_time-start_time
       
       ! * Compute tpsd step
       
       if(steps.gt.1) step = abs(dot_product(x-x0,xp-xp0)/(dot_product(xp-xp0,pxp-pxp0)+tiny(1.0_gp)))
       
       x0=x ; xp0=xp ; pxp0=pxp

       ! * Trust region

       strun=strun+(maxval(abs(step*pxp))-strun)/10
       
       if(steps.gt.1) step=min(step,min(0.1_gp,2*strun)/maxval(abs(pxp)))
       
       ! * Take a step
       
       x=x+step*pxp
       
       ! * Update strain, unit cell and positions

       strain=reshape(x(3*str%num_ions+1:3*str%num_ions+3*3),(/3,3/))

       strain=(strain+transpose(strain))/2
       
       str%lattice_car=matmul(strain,lat0)

       str%ion_positions=matmul(strain,reshape(x(1:3*str%num_ions),(/3,str%num_ions/)))

       call update_cell(str)

       ! * Running average

       de=maxval(abs(step*pxp*xp))
       
       derun=derun+(log10(de)-derun)/10

       ! * Track convergence
       
       if((mod(total_steps,track_every).eq.0).and.track) then
          write (unit_track,'(i6,*(1x,f20.10))') &
               total_steps,((e+p*str%volume)-h1)/str%num_ions,derun
          call flush(unit_track)
          call write_xyze(str)
       end if

    end do
     
    call eval_ddp(str,e,g,sig,var)
    
    if(track) call write_xyze(str)
    
    str%stress=sig/str%volume
    str%energy=e
    str%enthalpy=str%energy+p*str%volume
    str%variance=var

    return

998 stop 'There is a problem opening the track/xyze file. Stopping.'

  end subroutine opt_tpsd

  subroutine opt_newton(str,steps,tol,cell_fix)

    !-------------------------------------------------------------!

    type(structure), intent(inout) :: str
    integer,         intent(inout) :: steps
    real(kind=gp),      intent(in) :: tol
    logical,            intent(in) :: cell_fix

    real(kind=gp) :: e,e0,gn,p,sig(3,3),lambda,lambda_min

    real(kind=gp), save, allocatable :: g0(:,:),g(:,:),hess(:,:),gv(:),hessplus(:,:)

    integer :: maxsteps,n

    if(.not.cell_fix) stop 'variable cell not implemented'

    maxsteps=steps

    if(.not.allocated(g)) allocate(g0(3,str%num_ions),g(3,str%num_ions),hess(3*str%num_ions,3*str%num_ions),&
         hessplus(3*str%num_ions,3*str%num_ions),gv(3*str%num_ions))

    ! ** Open a file for the convergence data

    if(track.and.optimise) then

       write(trackfile,'(a)') trim(seedname)//".track"       
       open(unit=unit_track,file=trackfile,status="unknown",err=998)

       if(.not.quiet) then
          write(xyzefile,'(a)') trim(seedname)//".xyze"       
          open(unit=unit_xyze,file=xyzefile,status="unknown",err=998)
       end if

    end if

    p=(str%external_pressure(1,1)+str%external_pressure(2,2)+str%external_pressure(3,3))/3.0_gp

    gn         = 1.0_gp
    lambda     = 1e2_gp
    lambda_min = 1e-4_gp
    converged  = .true.
    e=0.0_gp
    e0=huge(1.0_gp)

    steps = 0
    do while(((abs(e-e0)/str%num_ions.gt.tol).or.(steps.lt.2)))

       if(steps.gt.maxsteps) then 
          converged=.false.
          exit
       end if

       call eval_ddp(str,e0,g0)

       if(total_steps.eq.0) h1=e+p*str%volume
       
       if(steps.lt.100) then
          hess=0.0_gp
       else
          if(mod(steps,10).eq.0) call get_hessian(str,hess)
       end if
    
       hessplus=real(hess,gp)
       do n=1,3*str%num_ions
          hessplus(n,n)=hess(n,n)+lambda
       end do
       
       if(.not.optimise) exit
       
       g0=reshape(matmul(inv(hessplus),reshape(g0,(/3*str%num_ions/))),(/3,str%num_ions/))

       str%ion_positions = str%ion_positions + g0

       call eval_ddp(str,e,g)
       
       if(e.gt.e0) then
          str%ion_positions=str%ion_positions-g0
          lambda=lambda*1.5_gp
       else
          lambda=max(lambda/2.0_gp,lambda_min)
          steps=steps+1
          total_steps=total_steps+1
          gn=norm2(reshape(g,(/3*str%num_ions/)))
          if(track) then
             write (unit_track,'(i0,*(1x,f0.10))') &
                  total_steps,((e+p*str%volume)-h1)/str%num_ions,log10(abs(e-e0)/str%num_ions),log10(gn/str%num_ions),log10(lambda)
             flush(unit_track)
             call write_xyze(str)         
          end if
       end if

    end do

    if(cell_fix) then
       call eval_ddp(str,e,g,sig)
       if(track) call write_xyze(str)
    end if

    str%stress=sig/str%volume
    str%energy=e
    str%enthalpy=str%energy+p*str%volume

    return

998 stop 'There is a problem opening the track/xyze file. Stopping.'

  end subroutine opt_newton
  
  subroutine opt_test(str)

    type(structure), intent(inout) :: str

    real(kind=gp) :: ion_frac(3,str%num_ions),ion_positions0(3,str%num_ions),lattice_car0(3,3)
    real(kind=gp) :: g(3,str%num_ions),gprime(3,str%num_ions),gm(3,str%num_ions),g0(3,str%num_ions),s(3,3),s0(3,3)

    real(kind=gp) :: step=1e-6_gp,e0,ep,em,p,strain(3,3)

    integer :: n,i,j

    ion_positions0 = str%ion_positions
    lattice_car0   = str%lattice_car

    call eval_ddp(str,e0,g0,s0)

    p=(str%external_pressure(1,1)+str%external_pressure(2,2)+str%external_pressure(3,3))/3.0_gp
    
    str%stress=s0/str%volume
    str%energy=e0
    str%enthalpy=e0+p*str%volume
    
    write (stdout,'(a11,f25.10)') 'Volume:    ',  str%volume
    write (stdout,'(a11,f25.10)') 'Pressure:  ',  (str%stress(1,1)+str%stress(2,2)+str%stress(3,3))/3.0_gp*160.21766208_gp
    write (stdout,'(a11,f25.10)') 'Energy:    ',  str%energy
    write (stdout,'(a11,f25.10)') 'Enthalpy:  ',  str%enthalpy
    write (stdout,*)
    flush (stdout)
    
    ! ** Test forces

    do n=1,str%num_ions
       do i=1,3

          str%ion_positions=ion_positions0

          str%ion_positions(i,n)=ion_positions0(i,n)+step

          call eval_ddp(str,ep,gprime)

          str%ion_positions(i,n)=ion_positions0(i,n)-step

          call eval_ddp(str,em,gm)

          g(i,n)=-(ep-em)/2/step

       end do

       write (stdout,'(3f13.5)') g0(:,n),g(:,n)
       write (stdout,*) '--'

    end do

    str%ion_positions=ion_positions0

    ion_frac=matmul(str%lattice_rec,str%ion_positions)

    ! ** Test stress

    do i=1,3
       do j=1,3

          ! ** Plus

          str%lattice_car=lattice_car0

          strain=ident

          strain(i,j)=strain(i,j)+step

          str%lattice_car=matmul(strain,str%lattice_car)

          call update_cell(str)

          str%ion_positions = matmul(str%lattice_car,ion_frac)

          call eval_ddp(str,ep,gprime)

          ! ** Minus

          str%lattice_car=lattice_car0

          strain=ident

          strain(i,j)=strain(i,j)-step

          str%lattice_car=matmul(strain,str%lattice_car)

          call update_cell(str)

          str%ion_positions = matmul(str%lattice_car,ion_frac)

          call eval_ddp(str,em,gm)

          s(i,j)=-(ep-em)/2/step

       end do
    end do

!!$    do i=1,3
!!$       do j=1,3
!!$
!!$          ! ** Plus
!!$
!!$          str%lattice_car=lattice_car0
!!$
!!$          str%lattice_car(i,j)=str%lattice_car(i,j)+step
!!$
!!$          call update_cell()
!!$
!!$          call eval_ddp(str,ep,gprime)
!!$
!!$          ep=ep+str%volume
!!$          
!!$          ! ** Minus
!!$
!!$          str%lattice_car=lattice_car0
!!$
!!$          str%lattice_car(i,j)=str%lattice_car(i,j)-step
!!$
!!$          call update_cell(str)
!!$
!!$          call eval_ddp(str,em,gm)
!!$
!!$          em=em+str%volume
!!$          
!!$          s(i,j)=-(ep-em)/2/step
!!$
!!$       end do
!!$    end do

    str%lattice_car=lattice_car0

    call update_cell(str)

    write (stdout,'(3f13.5)') s0
    !write (stdout,'(3f13.5)') s0-str%volume*transpose(str%lattice_rec)

    write (stdout,*) '--'

    write (stdout,'(3f13.5)') s

  end subroutine opt_test

  subroutine get_hessian(str,hess)

    type(structure), intent(inout) :: str
    real(kind=gp),     intent(out) :: hess(3*str%num_ions,3*str%num_ions)

    real(kind=gp) :: e,g(3,str%num_ions),gv0(3*str%num_ions),gvp(3*str%num_ions),gvm(3*str%num_ions)
    real(kind=gp) :: step=1e-4_gp,pos0(3,str%num_ions)

    integer :: ni,i,n

    pos0=str%ion_positions
    
    call eval_ddp(str,e,g)

    gv0=reshape(g,(/3*str%num_ions/))

    n=0
    do ni=1,str%num_ions
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

          hess(:,n)=(gvp-gvm)/2/step
          
       end do
    end do

    hess=-(hess+transpose(hess))/2
    
  end subroutine get_hessian
  
  subroutine precondition(str,v,vp,scale)

    type(structure), intent(in) :: str
    real(kind=gp),  intent(in)  :: v(3*str%num_ions)
    real(kind=gp),  intent(out) :: vp(3*str%num_ions)

    real(kind=gp), intent(in)  :: scale

    real(kind=gp) :: vp0(3*str%num_ions),g(3*str%num_ions),g0(3*str%num_ions)
    
    integer :: n,ni,nj

    integer, allocatable, dimension(:,:), save :: A,Lmat
    
    real(kind=gp) :: rad2,rninj2,vec(3),step

    if(.not.allocated(A)) then
       allocate(A(3*str%num_ions,3*str%num_ions),Lmat(str%num_ions,str%num_ions))
    end if
    
    if(scale.lt.0.0_gp) then

       rad2=huge(1.0_gp)
       do ni=1,str%num_ions
          do nj=ni+1,str%num_ions
             vec(:)=str%ion_positions(:,nj)-str%ion_positions(:,ni)
             rninj2=(vec(1)**2+vec(2)**2+vec(3)**2)
             if(rninj2.lt.rad2)rad2=rninj2
          end do
       end do

       rad2=rad2*abs(scale)

    end if
        
    Lmat=0

    do ni=1,str%num_ions
       do nj=ni+1,str%num_ions
          if(scale.gt.0.0_gp) rad2=(str%ion_rad(ni)+str%ion_rad(nj))**2*scale
          Lmat(ni,nj)=-aij(str,ni,nj,rad2)
          Lmat(nj,ni)=Lmat(ni,nj)
       end do
    end do

    do ni=1,str%num_ions
       Lmat(ni,ni)=-sum(Lmat(ni,:))
    end do

    A=universal_P(str,10*Lmat,1)

    vp=sparse_mult(A,v,init=.true.)
        
    step=1e-4_gp
    vp=v
    n=0
    g=0.0_gp
    do while((norm2(g).gt.1e-6_gp).or.(n.eq.0))
       n=n+1
       if(n.gt.10000) then
          vp=v
          exit
       end if
       g0=g       
       g=(sparse_mult(A,vp)/10-v)     
       if(n.gt.1) step=abs(dot_product(vp-vp0,g-g0)/(dot_product(g-g0,g-g0)+tiny(1.0_gp)))
       vp0=vp
       vp=vp-step*g
    end do
    
  end subroutine precondition
  
  function universal_P(str,L,stab)
    
    type(structure), intent(in) :: str
    integer, dimension(str%num_ions,str%num_ions), intent(in) :: L
    integer, intent(in) :: stab

    integer, dimension(3*str%num_ions,3*str%num_ions) :: universal_P

    integer :: ni,nj,i,j,n,m

    universal_P=0
    n=0
    do ni=1,str%num_ions
       do i=1,3
          n=n+1
          m=0
          do nj=1,str%num_ions
             do j=1,3
                m=m+1
                if(i.eq.j) universal_P(m,n)=L(nj,ni)
             end do
          end do
       end do
    end do

    ! ** Add stabilisation

    do n=1,3*str%num_ions
       universal_P(n,n)=universal_P(n,n)+stab
    end do

  end function universal_P

  function aij(str,i,j,r2)

    type(structure), intent(in) :: str

    integer, intent(in) :: i,j
    real(kind=gp), intent(in) :: r2

    integer :: aij

    integer :: n1,n2,n3

    real(kind=gp) :: vec(3),rij2,rsq

    vec(:)=str%ion_positions(:,j)-str%ion_positions(:,i)
    
    rij2=huge(1.0_gp)
    do n1=-1,1
       do n2=-1,1
          do n3=-1,1
             
             rsq=norm2(vec+matmul(str%lattice_car,(/n1,n2,n3/)))**2

             if(rsq.lt.rij2) rij2=rsq
             
          end do
       end do
    end do
    
    if(rij2.lt.r2) then
       aij=1
    else
       aij=0
    end if

  end function aij

  function sparse_mult(M,v,init)

    integer, dimension(:,:), intent(in) :: M
    real(kind=gp), dimension(:), intent(in) :: v
    logical, optional, intent(in) :: init

    real(kind=gp), dimension(size(M,2)) :: sparse_mult

    integer, save, allocatable, dimension(:,:) :: C,indx
    integer, save, allocatable, dimension(:) :: L

    integer :: i,j,k
    
    if(present(init)) then
       if(init) then
          if(allocated(C)) deallocate(C,indx,L)
          allocate(C(size(M,1),size(M,2)),L(size(M,1)),indx(size(M,1),size(M,2)))
          do i=1,size(M,1)
             k=0
             do j=1,size(M,2)
                if(M(j,i).ne.0) then
                   k=k+1
                   C(k,i)=M(j,i)
                   indx(k,i)=j
                end if
             end do
             L(i)=k
          end do
       end if
    else
       if(.not.allocated(C)) stop 'C not allocated'

       sparse_mult=0.0_gp
       do i=1,size(M,1)
          do j=1,L(i)
             sparse_mult(i)=sparse_mult(i)+C(j,i)*v(indx(j,i))
          end do
       end do
       
    end if
    
  end function sparse_mult

  function inv(A) result(Ainv)

    real(kind=gp), dimension(:,:), intent(in) :: A
    real(kind=gp), dimension(size(A,1),size(A,2)) :: Ainv

    real(kind=gp), dimension(size(A,1)) :: work 
    integer,       dimension(size(A,1)) :: ipiv
    integer :: n, info

    external dgetrf,sgetrf
    external dgetri,sgetri

    Ainv=A
    n=size(A,1)

    if(gp.eq.sp) then
       call sgetrf(n,n,Ainv,n,ipiv,info)
    else
       call dgetrf(n,n,Ainv,n,ipiv,info)
    end if

    if (info.ne.0)  stop 'inv: matrix is numerically singular!'

    if(gp.eq.sp) then
       call sgetri(n,Ainv,n,ipiv,work,n,info)
    else
       call dgetri(n,Ainv,n,ipiv,work,n,info)
    end if

    if (info.ne.0) stop 'inv: matrix inversion failed!'

  end function inv
  
end module opt
