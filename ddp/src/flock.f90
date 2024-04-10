!==================================================================================!
!                                    flock                                         !
!==================================================================================!
!                                                                                  !
! This file is part of the ddp package.                                            !
!                                                                                  !
! nn is free software; you can redistribute it and/or modify it under the terms    !
! of the GNU General Public License version 2 as published by the Free Software    !
! Foundation.                                                                      !
!                                                                                  !
! This program is distributed in the hope that it will be useful, but WITHOUT ANY  !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  !
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.        !           
!                                                                                  !
! You should have received a copy of the GNU General Public License along with this!
! program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street,!                   
! Fifth Floor, Boston, MA  02110-1301, USA.                                        !
!                                                                                  !
!----------------------------------------------------------------------------------!
! This program combines ddps using non-negative least squares                      !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
program flock

  use constants
  use omp_lib
  use data
  use nn
#ifdef IFORT
  use ifport, only : getpid
#endif
  
  implicit none

  integer, parameter :: unit_delta=46
  integer, parameter :: unit_agr=47

  type(network), allocatable, dimension(:) :: ddp

  type(dataset) :: training,validation,testing

  type(metadata) :: metatrain,metavalid,metatest

  integer, allocatable, dimension(:) :: seed
  
  integer :: n,nddp,nomp

  real(kind=gp), allocatable, dimension(:) :: weight

  real(kind=gp) :: lambda

  character(len=20)  :: scheme
  character(len=120) :: ddplist(10000)
  character(len=240) :: seedname

  logical :: quiet=.false.,plot=.true.,pdf=.false.,valid=.false.

  call init_pseudorandom()
  
  call get_arguments()

  call omp_set_num_threads(nomp) 

  call banner()

  call load_data()

  call load_ddp()

  allocate(weight(nddp))

  write (stderr,*)
  write (stderr,*) 'combining ..'
  
  if(valid) then
     call flock_combine(ddp,weight,validation,metavalid)
  else
     call flock_combine(ddp,weight,training,metatrain)
  end if
  
  write (stderr,*)
  write (stderr,*) 'testing ..'
  
  call flock_test(' training:   ',ddp,weight,training,metatrain)
  call flock_test(' validation: ',ddp,weight,validation,metavalid)
  call flock_test(' testing:    ',ddp,weight,testing,metatest)

  write (stderr,*)
  write (stderr,*) '--------------------------------------------'

  do n=1,nddp
     if(abs(weight(n)).gt.epsilon(1.0_gp)) write (stdout,'(a,f25.10)') trim(ddplist(n)),weight(n)
  end do

contains

  subroutine load_data()

    if(.not.quiet) write (stderr,'(a)') ' loading training data from: '//'training'
    call data_read('training',training,metatrain)
    if(.not.quiet) write (stderr,'(a,i10)') ' structures found:   ',metatrain%nstruct
    if(.not.quiet) write (stderr,'(a,i10)') ' centres found:      ',training%num
    if(.not.quiet) write (stderr,*)

    if(.not.quiet) write (stderr,'(a)') ' loading validation data from: '//'validation'
    call data_read('validation',validation,metavalid)
    if(.not.quiet) write (stderr,'(a,i10)') ' structures found:   ',metavalid%nstruct
    if(.not.quiet) write (stderr,'(a,i10)') ' centres found:      ',validation%num
    if(.not.quiet) write (stderr,*)
    
    if(.not.quiet) write (stderr,'(a)') ' loading testing data from: '//'testing'
    call data_read('testing',testing,metatest)
    if(.not.quiet) write (stderr,'(a,i10)') ' structures found:   ',metatest%nstruct
    if(.not.quiet) write (stderr,'(a,i10)') ' centres found:      ',testing%num
    if(.not.quiet) write (stderr,*)

  end subroutine load_data

  subroutine load_ddp()

    nddp=0
    do
       nddp=nddp+1
       read(stdin,'(a)',err=999,end=99) ddplist(nddp)
    end do
99  continue

    nddp=nddp-1

    allocate(ddp(nddp))

    do n=1,nddp
       call nn_load(ddp(n),ddplist(n))
    enddo

    write (stderr,'(a,i0,a)') ' loaded ',nddp,' data derived potentials'

    return

999 write (stderr,*) 'problem reading ddp list'

  end subroutine load_ddp

  subroutine flock_combine(net,wgt,dat,met)

    type(network), dimension(:), intent(inout) :: net
    real(kind=gp), dimension(:), intent(inout) :: wgt
    type(dataset),               intent(in)    :: dat
    type(metadata),              intent(in)    :: met

    real(kind=gp), allocatable :: Ep(:,:),E0(:),es(:)

    integer :: ns,nv,nd,nc

    allocate(Ep(met%nstruct,size(net)),E0(met%nstruct),es(1:1))
        
    if(scheme.eq.'mean') then
       wgt=1.0_gp/real(nddp,gp)
       return
    end if
    
    ! ** Set target

    E0=met%enthalpy

    ! ** Set matrix of predictions

    Ep=0.0_gp

    nv=0
    do ns=1,met%nstruct

       do nc=1,met%ncentres(ns) 
          nv=nv+1

          do nd=1,size(net)
             
             es(1:1)=nn_forward(net(nd),dat%in(:,nv),scale=.true.)

             Ep(ns,nd)=Ep(ns,nd)+es(1)+external_pressure*met%volume(ns)/met%ncentres(ns)

          end do

       end do

    end do
    
    select case(scheme)
    case('best')
       call best(wgt,Ep,E0)
    case('linear','ls')
       call linear(wgt,Ep,E0)
    case('nnls')
       call nnls(wgt,Ep,E0)
    case default
       stop 'flock_combine: scheme not known'
    end select
    
  end subroutine flock_combine

  subroutine linear(x,A,b)

    real(kind=gp), dimension(:),  intent(out) :: x
    real(kind=gp), dimension(:,:), intent(in) :: A
    real(kind=gp), dimension(:),   intent(in) :: b
    
    real(kind=gp) :: M(size(x),size(x))

    integer :: n
    
    M=matmul(transpose(A),A)
    do n=1,size(x)
       M(n,n)=M(n,n)+lambda
    end do
    x=nn_solve(M,matmul(transpose(A),b))  
    
  end subroutine linear
  
  subroutine best(x,A,b)

    real(kind=gp), dimension(:),  intent(out) :: x
    real(kind=gp), dimension(:,:), intent(in) :: A
    real(kind=gp), dimension(:),   intent(in) :: b

    real(kind=gp) :: c(size(x))

    do n=1,size(x)

       x=0.0_gp
       x(n)=1.0_gp
       
       c(n)=norm2(matmul(A,x)-b)

    end do

    x=0.0_gp
    x(minloc(c))=1.0_gp
    
  end subroutine best
  
  subroutine nnls(x,A,b)

    real(kind=gp), dimension(:),  intent(out) :: x
    real(kind=gp), dimension(:,:), intent(in) :: A
    real(kind=gp), dimension(:),   intent(in) :: b

    real(kind=gp) :: delta

    real(kind=gp), allocatable :: v(:),g(:),r(:),h(:),step(:),At(:,:),M(:,:),AtA(:,:)
    
    integer :: n,i

    allocate(v(size(x)),g(size(x)),r(size(b)),h(size(x)))
    allocate(step(size(x)),At(size(x),size(b)),M(size(x),size(x)),AtA(size(x),size(x)))
    
    
    At=transpose(A)

    AtA=matmul(At,A)

    call random_number(v)

    v=v/real(size(x),gp)
    
    delta=huge(1.0_gp)
    
    n=0
    
    do while(abs(delta).gt.epsilon(1.0_gp))

       n=n+1
       
       r=matmul(A,v**2)-b
       
       h=matmul(At,r)
       
       g=2*v*h

       M=4*outer(v,v)*AtA
   
       do i=1,size(x)
          M(i,i)=M(i,i)+2*h(i)
       end do

       step=-nn_solve(M,g)

       delta=dot_product(step,g)
       
       v=v-step*sign(1.0_gp,delta)

!!$       write (66,*) n,log10(norm2(r)),log10(abs(delta))
!!$       flush(66)
       
    end do
    
    x=v**2
    
    where(x<1e-13_gp) x=0.0_gp
    
  end subroutine nnls
  
  subroutine banner()

    if(quiet) return

    write (stderr,'(a)') "                                              "                           
    write (stderr,'(a)') "         ▄████  █    ████▄ ▄█▄    █  █▀       "
    write (stderr,'(a)') "         █▀   ▀ █    █   █ █▀ ▀▄  █▄█         "
    write (stderr,'(a)') "         █▀▀    █    █   █ █   ▀  █▀▄         "
    write (stderr,'(a)') "         █      ███▄ ▀████ █▄  ▄▀ █  █        "
    write (stderr,'(a)') "          █         ▀      ▀███▀    █         "
    write (stderr,'(a)') "           ▀                       ▀          "
    write (stderr,'(a)') "                                              "                           
    write (stderr,'(a)') " Author: Chris J. Pickard, Cambridge, 2021 (c)"
    write (stderr,*)
    write (stderr,'(a,*((i0,1x)))')    ' num threads    : ', omp_get_max_threads()
    write (stderr,'(a,*((i0,1x)))')    ' max threads    : ', omp_get_num_procs()
    write (stderr,'(a,a)')             ' optimisation   : ', scheme
    if(valid) then
       write (stderr,'(a,a)')             ' set            : ', 'validation'       
    else
       write (stderr,'(a,a)')             ' set            : ', 'training'
    end if
    write (stderr,'(a,*((g0.3,1x)))')  ' lambda         : ', lambda
    write (stderr,*)

  end subroutine banner

  subroutine get_arguments()

    character(len=80), allocatable,dimension(:) :: cargs

    integer :: i,iargc,num_args

    num_args = iargc()

    allocate(cargs(num_args))

    do i=1,num_args
       call getarg(i,cargs(i))
    end do

    if(any(cargs.eq."-h")) goto 100

    nomp=1
    if(any(cargs.eq."-ompnp")) then               ! * Number of omp threads (note - not mpi)
       do i=1,num_args-1
          if(cargs(i).eq."-ompnp") exit
       end do
       read(cargs(i+1),*) nomp
       nomp=min(nomp,omp_get_num_procs())
    end if

    seedname='flock'
    if(any(cargs.eq."-s")) then                ! * Seedname
       do i=1,num_args-1
          if(cargs(i).eq."-s") exit
       end do
       if(i.lt.num_args) then
          if(cargs(i+1)(1:1).ne.'-') read(cargs(i+1),'(a)') seedname
       end if
    end if

    lambda=1.0_gp
    if(any(cargs.eq."-lambda")) then           ! * Regularisation
       do i=1,num_args-1
          if(cargs(i).eq."-lambda") exit
       end do
       if(i.lt.num_args) then
          if(cargs(i+1)(1:1).ne.'-') read(cargs(i+1),*) lambda
       end if
    end if

    external_pressure=0.0_gp
    if(any(cargs.eq."-p")) then                ! * External pressure
       do i=1,num_args-1
          if(cargs(i).eq."-p") exit
       end do
       if(i.lt.num_args) then
          if(cargs(i+1)(1:1).ne.'-') read(cargs(i+1),*) external_pressure
       end if
    end if
    external_pressure=external_pressure/160.21766208_gp ! * Convert to eV/ang^3
    
    scheme='nnls'
    if(any(cargs.eq."-o")) then                ! * Scheme
       do i=1,num_args-1
          if(cargs(i).eq."-o") exit
       end do
       if(i.lt.num_args) then
          if(cargs(i+1)(1:1).ne.'-') read(cargs(i+1),'(a)') scheme
       end if
    end if

    if(any(cargs.eq."-v")) valid=.true.
    if(any(cargs.eq."-np")) plot=.false.
    if(any(cargs.eq."-pdf")) pdf=.true.
    if(any(cargs.eq."-q")) quiet=.true.

    deallocate(cargs)

    return

100 write (stderr,*) 'Usage: flock [-ompnp] [-s] [-lambda] [-o] [-v] [-p] [-np] [-pdf] [-q] [-h]'
    write (stderr,*)
    write (stderr,*) '-ompnp    : Number of omp threads (not mpi)'
    write (stderr,*) '-s        : Seedname'
    write (stderr,*) '-lambda   : Regularisation parameter'
    write (stderr,*) '-o        : Optimisation scheme'
    write (stderr,*) '-v        : Optimise using validation dataset'
    write (stderr,*) '-p        : External pressure'
    write (stderr,*) '-np       : No plotting'
    write (stderr,*) '-pdf      : Output plots to pdf (requires grace)'
    write (stderr,*) '-q        : Quiet - minimal output'
    write (stderr,*) '-h        : Print this help message'
    write (stderr,*)
    write (stderr,*) 'Note: a list of ddp files is read from STDIN'

    stop

  end subroutine get_arguments

  subroutine flock_test(ctask,net,wgt,dat,met)

    character(len=*), intent(in)  :: ctask
    type(network), dimension(:),   intent(inout) :: net
    real(kind=gp), dimension(:), intent(in) :: wgt
    type(dataset),   intent(in)    :: dat
    type(metadata),  intent(in)    :: met

    real(kind=gp) :: e_predict(met%nstruct),es(1),mae,rmse,emax,abserr,tick,dmax,crmse(met%nstruct),cmae(met%nstruct)

    integer :: ns,nc,nv,nd,stat,indx(met%nstruct),i

    character(len=40)  :: ctemp
    character(len=240) :: ctemp2
    character(len=240) :: cmax

    if(plot) then

       open(unit=unit_delta,file=trim(seedname)//'-'//trim(met%setname)//'.delta',form='formatted',status='unknown')
       open(unit=unit_agr,file=trim(seedname)//'-'//trim(met%setname)//'.agr',form='formatted',status='unknown')

       tick=(maxval(met%enthalpy/met%natoms)-minval(met%enthalpy/met%natoms))/5.0_gp

       call agr_head(minval(met%enthalpy/met%natoms),maxval(met%enthalpy/met%natoms),&
            minval(met%enthalpy/met%natoms),maxval(met%enthalpy/met%natoms),tick,tick)
       
       write (unit_agr,*)
       write (unit_agr,*) minval(met%enthalpy/met%natoms),minval(met%enthalpy/met%natoms)
       write (unit_agr,*) maxval(met%enthalpy/met%natoms),maxval(met%enthalpy/met%natoms)
       write (unit_agr,*)

       write (unit_agr,*) '&'

    end if

    nv=0
    do ns=1,met%nstruct

       es(1:1)=0.0_gp

       do nc=1,met%ncentres(ns) 
          nv=nv+1

          do nd=1,size(net)

             if(abs(wgt(nd)).gt.0.0_gp) es(1:1)=es(1:1)+nn_forward(net(nd),dat%in(:,nv),scale=.true.)*wgt(nd)

          end do

       end do

       e_predict(ns)=es(1)+external_pressure*met%volume(ns)

       if(plot) write (unit_agr,*) met%enthalpy(ns)/met%natoms(ns),e_predict(ns)/met%ncentres(ns)

    end do

    if(plot) then
       write (unit_agr,*) '&'

       dmax=1.1_gp*maxval(abs(met%enthalpy/met%natoms-e_predict/met%ncentres))
       call delta_head(minval(met%enthalpy/met%natoms),maxval(met%enthalpy/met%natoms),&
            0.0_gp,dmax,tick,dmax/5.0_gp)
       
    end if

    do i=1,met%nstruct
       indx(i)=i
    end do
    
    call heap_sort_index(met%nstruct,met%enthalpy/met%natoms,indx)
    
    emax=0.0_gp
    mae=0.0_gp
    rmse=0.0_gp
    do ns=1,met%nstruct
       abserr=abs(e_predict(indx(ns))/met%ncentres(indx(ns))-met%enthalpy(indx(ns))/met%natoms(indx(ns)))
       rmse=rmse+abserr**2
       mae=mae+abserr
       if(plot) write (unit_delta,*) &
            met%enthalpy(indx(ns))/met%natoms(indx(ns)),abs(e_predict(indx(ns))/met%ncentres(indx(ns))-&
            met%enthalpy(indx(ns))/met%natoms(indx(ns)))
       cmae(indx(ns))=mae/ns
       crmse(indx(ns))=sqrt(rmse/ns)
       if(abserr.gt.emax) then
          emax=abserr
          cmax=met%label(indx(ns))
       end if
    end do
    
    rmse=sqrt(rmse/real(met%nstruct,gp))
    mae=mae/real(met%nstruct,gp)
    
    if(plot) then
       write (unit_delta,*) '&'
       do ns=1,met%nstruct
          write (unit_delta,*) met%enthalpy(indx(ns))/met%natoms(indx(ns)),crmse(indx(ns))
       end do
       write (unit_delta,*) '&'
       write (unit_delta,*) minval(met%enthalpy/met%natoms),rmse
       write (unit_delta,*) maxval(met%enthalpy/met%natoms),rmse
       write (unit_delta,*) '&'
       do ns=1,met%nstruct
          write (unit_delta,*) met%enthalpy(indx(ns))/met%natoms(indx(ns)),cmae(indx(ns))
       end do
       write (unit_delta,*) '&'
       write (unit_delta,*) minval(met%enthalpy/met%natoms),mae
       write (unit_delta,*) maxval(met%enthalpy/met%natoms),mae
       write (unit_delta,*) '&'

    end if

    write (stderr,*)
    write (ctemp,'(a)') ctask//trim(met%setname)//' RMSE/MAE:  '
    write (stderr,'(a,2f9.2,a)',advance='no') ctemp,rmse*1000.0_gp,mae*1000.0_gp,' meV'
    write (stderr,'(a,f7.5)',advance='no')  ' Spearman : ',spearman(e_predict/met%ncentres,met%enthalpy/met%natoms)
    write (stderr,'(a,f10.2,a,a)') ' Max : ', emax*1000.0_gp,' meV for ',trim(cmax)

    emax=0.0_gp
    mae=0.0_gp
    rmse=0.0_gp
    do ns=1,met%nstruct
       abserr=abs(e_predict(ns)/met%ncentres(ns)*met%natoms(ns)-met%enthalpy(ns)) 
       rmse=rmse+abserr**2
       mae=mae+abserr
       if(abserr.gt.emax) then
          emax=abserr
          cmax=met%label(ns)
       end if
    end do
    rmse=sqrt(rmse/real(met%nstruct,gp))
    mae=mae/real(met%nstruct,gp)
    write (ctemp,'(a)') '          ( '//trim(met%setname)//' RMSE/MAE:  '
    write (stderr,'(a,2f9.2,a)',advance='no') ctemp,rmse,mae,'  eV'
    write (stderr,'(a,f7.5)',advance='no')  ' Spearman : ',spearman(e_predict/met%ncentres*met%natoms,met%enthalpy)
    write (stderr,'(a,f10.2,a,a,a)') ' Max : ', emax,'  eV for ',trim(cmax),' )'

    if(plot) then

       close(unit_delta)
       close(unit_agr)

       if(pdf) then

          ctemp2="gracebat "//trim(seedname)//'-'//trim(met%setname)//".agr -hdevice PDF -hardcopy -printfile "&
               //trim(seedname)//'-'//trim(met%setname)//".pdf"
          stat=0
          call system(ctemp2,stat)
          ctemp2='Problem executing external command :: '//trim(ctemp2)
          if (stat.ne.0) then
             write (stderr,'(a)') trim(ctemp2)
             stop
          end if
          if(.not.quiet) write (stderr,*) " saved file "//trim(seedname)//'-'//trim(met%setname)//".pdf"

       end if

    end if

  end subroutine flock_test

  subroutine agr_head(xmin,xmax,ymin,ymax,xtick,ytick)

    real(kind=gp), intent(in) :: xmin,xmax,ymin,ymax,xtick,ytick

    write(unit_agr,'(a)') '@version 50109'
    write(unit_agr,'(a)') '@default linewidth 2.0'
    write(unit_agr,'(a)') '@default linestyle 1'
    write(unit_agr,'(a)') '@g0 on'
    write(unit_agr,'(a)') '@with g0'
    write(unit_agr,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_agr,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_agr,'(a)') '@map color 0 to (255, 255, 255), "white"'
    write(unit_agr,'(a)') '@map color 1 to (0, 0, 0), "black"'
    write(unit_agr,'(a)') '@map color 2 to (228, 26, 28), "red"'
    write(unit_agr,'(a)') '@map color 3 to (55, 126, 184), "blue"'
    write(unit_agr,'(a)') '@map color 4 to (77, 175, 74), "green"'
    write(unit_agr,'(a)') '@map color 5 to (152, 78, 163), "purple"'
    write(unit_agr,'(a)') '@map color 6 to (255, 127, 0), "orange"'
    write(unit_agr,'(a)') '@map color 7 to (255, 255, 51), "yellow"'
    write(unit_agr,'(a)') '@map color 8 to (166, 86, 40), "brown"'
    write(unit_agr,'(a)') '@map color 9 to (247, 129, 191), "pink"'
    write(unit_agr,'(a)') '@map color 10 to (153, 153, 153), "grey"'
    write(unit_agr,'(a)') '@map color 11 to (166, 206, 227), "lightblue"'
    write(unit_agr,'(a)') '@map color 12 to (178, 223, 138), "lightgreen"'
    write(unit_agr,'(a)') '@map color 13 to (251, 154, 153), "lightred"'
    write(unit_agr,'(a)') '@map color 14 to (253, 191, 111), "lightorange"'
    write(unit_agr,'(a)') '@map color 15 to (202, 178, 214), "lightpurple"'
    write(unit_agr,'(a,f12.3)') '@    world xmin',xmin
    write(unit_agr,'(a,f12.3)') '@    world xmax',xmax
    write(unit_agr,'(a,f12.3)') '@    world ymin',ymin
    write(unit_agr,'(a,f12.3)') '@    world ymax',ymax
    write(unit_agr,'(a)') '@    view xmin 0.200000'
    write(unit_agr,'(a)') '@    view xmax 0.900000'
    write(unit_agr,'(a)') '@    view ymin 0.200000'
    write(unit_agr,'(a)') '@    view ymax 0.900000'
    write(unit_agr,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_agr,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_agr,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_agr,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_agr,'(a)') '@    xaxis  label "Enthalpy (eV)"'
    write(unit_agr,'(a)') '@    xaxis  label font 4'
    write(unit_agr,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_agr,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_agr,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_agr,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_agr,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_agr,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_agr,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_agr,'(a)') '@    yaxis  label "Predicted Enthalpy (eV)"'
    write(unit_agr,'(a)') '@    yaxis  label font 4'
    write(unit_agr,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_agr,'(a)') '@    s1 symbol 1'
    write(unit_agr,'(a)') '@    s1 symbol size 0.16000'
    write(unit_agr,'(a)') '@    s1 symbol color 1'
    write(unit_agr,'(a)') '@    s1 symbol pattern 1'
    write(unit_agr,'(a)') '@    s1 symbol fill color 2'
    write(unit_agr,'(a)') '@    s1 symbol fill pattern 1'
    write(unit_agr,'(a)') '@    s1 symbol linewidth 0.25'
    write(unit_agr,'(a)') '@    s1 symbol linestyle 1'
    write(unit_agr,'(a)') '@    s1 line type 0'
    write(unit_agr,'(a)') '@    s0 line type 1'
    write(unit_agr,'(a)') '@    s0 line linestyle 1'
    write(unit_agr,'(a)') '@    s0 line linewidth 1.0'
    write(unit_agr,'(a)') '@    s0 line color 2'
    write(unit_agr,'(a)') '@    s0 line pattern 1'


  end subroutine agr_head

  subroutine delta_head(xmin,xmax,ymin,ymax,xtick,ytick)

    real(kind=gp), intent(in) :: xmin,xmax,ymin,ymax,xtick,ytick

    write(unit_delta,'(a)') '@version 50109'
    write(unit_delta,'(a)') '@default linewidth 2.0'
    write(unit_delta,'(a)') '@default linestyle 1'
    write(unit_delta,'(a)') '@g0 on'
    write(unit_delta,'(a)') '@with g0'
    write(unit_delta,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_delta,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_delta,'(a)') '@map color 0 to (255, 255, 255), "white"'
    write(unit_delta,'(a)') '@map color 1 to (0, 0, 0), "black"'
    write(unit_delta,'(a)') '@map color 2 to (228, 26, 28), "red"'
    write(unit_delta,'(a)') '@map color 3 to (55, 126, 184), "blue"'
    write(unit_delta,'(a)') '@map color 4 to (77, 175, 74), "green"'
    write(unit_delta,'(a)') '@map color 5 to (152, 78, 163), "purple"'
    write(unit_delta,'(a)') '@map color 6 to (255, 127, 0), "orange"'
    write(unit_delta,'(a)') '@map color 7 to (255, 255, 51), "yellow"'
    write(unit_delta,'(a)') '@map color 8 to (166, 86, 40), "brown"'
    write(unit_delta,'(a)') '@map color 9 to (247, 129, 191), "pink"'
    write(unit_delta,'(a)') '@map color 10 to (153, 153, 153), "grey"'
    write(unit_delta,'(a)') '@map color 11 to (166, 206, 227), "lightblue"'
    write(unit_delta,'(a)') '@map color 12 to (178, 223, 138), "lightgreen"'
    write(unit_delta,'(a)') '@map color 13 to (251, 154, 153), "lightred"'
    write(unit_delta,'(a)') '@map color 14 to (253, 191, 111), "lightorange"'
    write(unit_delta,'(a)') '@map color 15 to (202, 178, 214), "lightpurple"'
    write(unit_delta,'(a,f12.3)') '@    world xmin',xmin
    write(unit_delta,'(a,f12.3)') '@    world xmax',xmax
    write(unit_delta,'(a,f12.3)') '@    world ymin',ymin
    write(unit_delta,'(a,f12.3)') '@    world ymax',ymax
    write(unit_delta,'(a)') '@    view xmin 0.200000'
    write(unit_delta,'(a)') '@    view xmax 0.900000'
    write(unit_delta,'(a)') '@    view ymin 0.200000'
    write(unit_delta,'(a)') '@    view ymax 0.900000'
    write(unit_delta,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_delta,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_delta,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_delta,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_delta,'(a)') '@    xaxis  label "Enthalpy (eV)"'
    write(unit_delta,'(a)') '@    xaxis  label font 4'
    write(unit_delta,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_delta,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_delta,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_delta,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_delta,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_delta,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_delta,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_delta,'(a)') '@    yaxis  label "Error (eV)"'
    write(unit_delta,'(a)') '@    yaxis  label font 4'
    write(unit_delta,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_delta,'(a)') '@    s0 symbol 1'
    write(unit_delta,'(a)') '@    s0 symbol size 0.16000'
    write(unit_delta,'(a)') '@    s0 symbol color 1'
    write(unit_delta,'(a)') '@    s0 symbol pattern 1'
    write(unit_delta,'(a)') '@    s0 symbol fill color 2'
    write(unit_delta,'(a)') '@    s0 symbol fill pattern 1'
    write(unit_delta,'(a)') '@    s0 symbol linewidth 0.25'
    write(unit_delta,'(a)') '@    s0 symbol linestyle 1'
    write(unit_delta,'(a)') '@    s0 line type 0'
    
    write(unit_delta,'(a)') '@    s1 line type 1'
    write(unit_delta,'(a)') '@    s1 line linestyle 2'
    write(unit_delta,'(a)') '@    s1 line linewidth 1.0'
    write(unit_delta,'(a)') '@    s1 line color 2'
    write(unit_delta,'(a)') '@    s1 line pattern 1'

    write(unit_delta,'(a)') '@    s2 line type 1'
    write(unit_delta,'(a)') '@    s2 line linestyle 1'
    write(unit_delta,'(a)') '@    s2 line linewidth 1.0'
    write(unit_delta,'(a)') '@    s2 line color 2'
    write(unit_delta,'(a)') '@    s2 line pattern 1'

    write(unit_delta,'(a)') '@    s3 line type 1'
    write(unit_delta,'(a)') '@    s3 line linestyle 2'
    write(unit_delta,'(a)') '@    s3 line linewidth 1.0'
    write(unit_delta,'(a)') '@    s3 line color 3'
    write(unit_delta,'(a)') '@    s3 line pattern 1'

    write(unit_delta,'(a)') '@    s4 line type 1'
    write(unit_delta,'(a)') '@    s4 line linestyle 1'
    write(unit_delta,'(a)') '@    s4 line linewidth 1.0'
    write(unit_delta,'(a)') '@    s4 line color 3'
    write(unit_delta,'(a)') '@    s4 line pattern 1'
    

  end subroutine delta_head
  
  ! ---------

  function spearman(x,y)

    real(kind=gp), dimension(:), intent(in) :: x,y

    real(kind=gp) :: spearman

    integer :: n,i,indx(size(x)),indy(size(y)),rnkx(size(x)),rnky(size(y))

    if(size(y).ne.size(x)) stop 'spearman called incorrectly'

    n=size(x)

    do i=1,n
       indx(i)=i
       indy(i)=i
    end do

    call heap_sort_index(n,x,indx)
    call heap_sort_index(n,y,indy)

    do i=1,n
       rnkx(indx(i))=i
       rnky(indy(i))=i
    end do

    spearman=0.0_gp
    do i=1,n
       spearman=spearman+real(rnkx(i)-rnky(i),gp)**2
    end do

    spearman=1.0_gp-spearman*6.0_gp/(real(n,gp)*(real(n,gp)**2-1.0_gp))

  end function spearman

  subroutine heap_sort_index(num_items,weight,indx)

    !=========================================================================!
    ! This subroutine sorts the list of weights into descending order.        !
    ! The weights are unchanged, with index returning the result              !
    !                                                                         !
    ! This is a heap sort                                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   num_items (input) :: The number of items to sort                      !
    !   weight (in)       :: The weights of each item.                        !
    !   indx (inout)      :: The indices of the weights in decending order    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Chris Pickard 26th July 2012                                 !
    !=========================================================================!

    implicit none

    ! Arguments

    integer, intent(in) :: num_items
    real(kind=gp), dimension(num_items), intent(in) :: weight
    integer, dimension(num_items), intent(inout) :: indx

    ! Local variables

    integer :: i,ir,j,l,indxa ! Loop counters
    real(kind=gp) :: wta
    real(kind=gp), dimension(num_items) :: wgtemp

    if(num_items.lt.2) return

    wgtemp=weight
!!$    do i=1,num_items
!!$       indx(i)=i
!!$    end do

    l=num_items/2+1
    ir=num_items

    do
       if(l.gt.1) then
          l=l-1
          wta=wgtemp(l)
          indxa=indx(l)
       else
          wta=wgtemp(ir)
          indxa=indx(ir)
          wgtemp(ir)=wgtemp(1)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1) then
             wgtemp(1)=wta
             indx(1)=indxa
             return
          end if
       end if
       i=l
       j=l+l
20     if(j.le.ir) then
          if(j.lt.ir) then
             if(wgtemp(j).lt.wgtemp(j+1)) j=j+1
          end if
          if(wta.lt.wgtemp(j)) then
             wgtemp(i)=wgtemp(j)
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          goto 20
       end if
       wgtemp(i)=wta
       indx(i)=indxa
    end do

  end subroutine heap_sort_index

  subroutine init_pseudorandom(sstring)

    integer, parameter                 :: maxr = HUGE(0) - 1
    integer(kind=8)                    :: count,countrate,countmax,pid
    integer                            :: i,j,iters,size
    real(kind=gp)                      :: rr
    character(len=*), optional         :: sstring

    pid=getpid()

    ! ** Seed the random number generator

    call random_seed(size=size)

    if(allocated(seed)) deallocate(seed)

    allocate(seed(size))

    if(present(sstring).and.(len_trim(sstring).gt.0)) then
       seed=0
       read(sstring,*,end=99,err=99) seed
       call random_seed(put=seed)
       return
99     stop 'init_pseudorandom : error reading seed'
    end if

    call random_seed(get=seed)

    do j=1,size
       call system_clock(count,countrate,countmax)
       iters = int(mod(count,pid),4)
       do i=1,iters+1
          call random_number(rr)
       end do
       seed(j) = int(2*(rr-0.5_gp)*maxr)
    end do

    call random_seed(put=seed)

  end subroutine init_pseudorandom

  function outer(a,b) result(c)

    real(kind=gp), intent(in) :: a(:),b(:)

    real(kind=gp) :: c(size(a),size(b))

    c=spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))

  end function outer

  
end program flock
