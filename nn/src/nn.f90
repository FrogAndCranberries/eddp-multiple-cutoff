!==================================================================================!
!                                      nn                                          !
!==================================================================================!
!                                                                                  !
! This file is part of the nn package.                                             !
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
! This module is a fortran implemention of neural networks                         !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
module nn

  implicit none

  private

  integer, parameter :: sp=4 ! Single = 4 Double = 8 Quad = 16
  integer, parameter :: dp=8
  integer, parameter :: qp=16
  
#ifdef SINGLE
  integer, parameter :: gp=sp
#else
  integer, parameter :: gp=dp
#endif
  
  integer, parameter :: stdin=5
  integer, parameter :: stdout=6
  integer, parameter :: stderr=0

  integer, parameter :: unit_nn=66
  integer, parameter :: unit_stop=55
  integer, parameter :: unit_track=44
  
  public :: nn_create
  public :: nn_initialise
  public :: nn_summary
  public :: nn_data
  public :: nn_normalise
  public :: nn_train_tpsd
  public :: nn_train_lm  
  public :: nn_train_galm  
  public :: nn_train_sgd
  public :: nn_train_adamw
  public :: nn_cost
  public :: nn_gradient
  public :: nn_test_gradient
  public :: nn_forward
  public :: nn_backward
  public :: nn_save
  public :: nn_load
  public :: nn_solve
  public :: nn_shuffle

  type, public :: network
     integer                                        :: depth
     integer                                        :: inputs
     integer                                        :: outputs
     integer                                        :: nwgts
     integer                                        :: nnds
     integer,           allocatable, dimension(:)   :: layers
     real(kind=gp),     allocatable, dimension(:)   :: weights
     real(kind=gp),     allocatable, dimension(:)   :: pre
     real(kind=gp),     allocatable, dimension(:)   :: act
     real(kind=gp),     allocatable, dimension(:)   :: out
     real(kind=gp),     allocatable, dimension(:,:) :: del
     real(kind=gp),     allocatable, dimension(:)   :: gra
     real(kind=gp),     allocatable, dimension(:)   :: in_mean
     real(kind=gp),     allocatable, dimension(:)   :: in_sd
     real(kind=gp),     allocatable, dimension(:)   :: out_mean
     real(kind=gp),     allocatable, dimension(:)   :: out_sd
     character(len=10), allocatable, dimension(:)   :: fnc
     character(len=2000)                            :: cmt
  end type network

  type, public :: dataset
     integer :: num, idim, odim
     real(kind=gp), allocatable, dimension(:,:) :: in
     real(kind=gp), allocatable, dimension(:,:) :: out
  end type dataset
  
contains

  subroutine nn_train_tpsd(net,train,test,thresh,ncmin,ncmax,gamma,track)

    type(network),           intent(inout) :: net
    type(dataset),           intent(in)    :: train
    type(dataset),           intent(in)    :: test
    real(kind=gp),           intent(in)    :: thresh
    integer,       optional, intent(in)    :: ncmin
    integer,       optional, intent(in)    :: ncmax
    real(kind=gp), optional, intent(in)    :: gamma
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: g(net%nwgts),gprime(net%nwgts),g0(net%nwgts),w0(net%nwgts)
    real(kind=gp) :: sgg,step(net%nwgts),cst,lgrun

    integer :: nc,ncmin_use,ncmax_use

    logical :: track_use,opened,stopfile

    open(unit=unit_stop,file="STOP",err=100)
100 close(unit_stop,status='delete')
    
    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file='nn.track',form='formatted',status='unknown')
    end if
    
    if(present(ncmax)) then
       ncmax_use=ncmax
    else
       ncmax_use=huge(1)
    end if

    if(present(ncmin)) then
       ncmin_use=ncmin
    else
       ncmin_use=1
    end if
    
    step=1e-4_gp
    sgg=huge(1.0_gp)
    lgrun=0.0_gp
    nc=0
    do while(((lgrun.gt.log10(thresh)).and.(nc.lt.ncmax_use)).or.(nc.lt.ncmin_use))
       nc=nc+1

       inquire(file="STOP",exist=stopfile)
       
       if(stopfile) exit
       
!!$       if(mod(nc-1,100).eq.0) then
!!$          nw=0
!!$          do nl=2,net%depth
!!$             write (55,*) nl,&
!!$                  sum(abs(net%weights(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl))))/((net%layers(nl-1)+1)*net%layers(nl))
!!$             !write (55,*) nl,norm2(net%weights(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl)))/(net%layers(nl-1)+1)*net%layers(nl)
!!$             nw=nw+(net%layers(nl-1)+1)*net%layers(nl)
!!$          end do
!!$          write (55,*)
!!$       endif
       
       g=nn_gradient(net,train,cst)
       if(present(gamma)) then
          g=g+net%weights*gamma
       end if

       ! * Precondition

       gprime=g*net%pre

       ! * Calculate step size

       if(nc.gt.1) &
          step=max(abs(dot_product(net%weights-w0,gprime-g0)/(dot_product(gprime-g0,gprime-g0)+tiny(1.0_gp))),epsilon(1.0_gp))
       
       ! * Store this step

       g0=gprime
       w0=net%weights

       ! * Compute estimate of change in cost

       sgg=dot_product(step*gprime,g)

       lgrun=lgrun+(log10(abs(sgg))-lgrun)/min(nc,10)
       
       ! * Update

       net%weights=net%weights-step*gprime

       if(track_use.and.(mod(nc-1,100).eq.0)) then
          write (unit_track,*) log10(real(nc,gp)),log10(cst),log10(nn_cost(net,test))
          flush(unit_track)
       end if

!!$       if(mod(nc-1,100).eq.0) then
!!$          nw=0
!!$          do nl=2,net%depth
!!$             write (66,*) nl,log10(norm2(g(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl))))
!!$             nw=nw+(net%layers(nl-1)+1)*net%layers(nl)
!!$          end do
!!$          write (66,*)
!!$       endif
       
    end do

    if(track_use) write (unit_track,*)

    write (stderr,*) 'Number of steps:   ',nc

  end subroutine nn_train_tpsd

  subroutine nn_train_lm(net,train,test,thresh,track)

    type(network),           intent(inout) :: net
    type(dataset),           intent(in)    :: train
    type(dataset),           intent(in)    :: test
    real(kind=gp),           intent(in)    :: thresh
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: step(net%nwgts),cst,cst0,cstt,lambda

    integer :: nc

    logical :: track_use,opened,stopfile

    open(unit=unit_stop,file="STOP",err=100)
100 close(unit_stop,status='delete')

    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file='nn.track',form='formatted',status='unknown')
    end if

    step=1.0_gp
    lambda=1e1_gp
    cst=huge(1.0_gp)
    cst0=0.0_gp

    nc=0
    do while((abs(cst-cst0).gt.thresh).or.(nc.lt.2))
       
       inquire(file="STOP",exist=stopfile)

       if(stopfile) exit

       ! * LM step

       step=nn_lmstep(net,train,lambda,cst0)
       
       cstt=nn_cost(net,test)
       
       ! * Update
       
       net%weights=net%weights+step

       cst=nn_cost(net,train)
       
       if(cst.gt.cst0) then
          net%weights=net%weights-step
          lambda=lambda*1.5_gp
       else
          lambda=lambda/2.0_gp
          nc=nc+1
          if(track_use) then
             write (unit_track,*) log10(real(nc,gp)),log10(cst0),log10(cstt),log10(norm2(step)),log10(lambda),log10(abs(cst-cst0))
             flush(unit_track)
          end if
       end if

    end do

    if(track_use) write (unit_track,*)

    write (stderr,*) 'Number of steps:   ',nc

  end subroutine nn_train_lm

  subroutine nn_train_galm(net,train,test,thresh,track)

    type(network),           intent(inout) :: net
    type(dataset),           intent(in)    :: train
    type(dataset),           intent(in)    :: test
    real(kind=gp),           intent(in)    :: thresh
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: step(net%nwgts),cst,cst0,cstt,lambda

    integer :: nc

    logical :: track_use,opened,stopfile

    open(unit=unit_stop,file="STOP",err=100)
100 close(unit_stop,status='delete')

    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file='nn.track',form='formatted',status='unknown')
    end if

    step=1.0_gp
    lambda=1e1_gp
    cst=huge(1.0_gp)
    cst0=0.0_gp

    nc=0
    do while((abs(cst-cst0).gt.thresh).or.(nc.lt.2))
       
       inquire(file="STOP",exist=stopfile)

       if(stopfile) exit

       ! * GALM step

       step=nn_galmstep(net,train,lambda,cst0)
       
       cstt=nn_cost(net,test)
       
       ! * Update
       
       net%weights=net%weights+step

       cst=nn_cost(net,train)
       
       if((cst.gt.cst0)) then
          net%weights=net%weights-step
          lambda=lambda*1.5_gp
       else
          lambda=lambda/2.0_gp
          nc=nc+1
          if(track_use) then
             write (unit_track,*) log10(real(nc,gp)),log10(cst0),log10(cstt),log10(norm2(step)),log10(lambda),log10(abs(cst-cst0))
             flush(unit_track)
          end if
       end if

    end do

    if(track_use) write (unit_track,*)

    write (stderr,*) 'Number of steps:   ',nc

  end subroutine nn_train_galm

  subroutine nn_train_sgd(net,train,test,nbatch,nsteps,step,gamma,track)

    type(network),           intent(inout) :: net
    type(dataset),           intent(in)    :: train
    type(dataset),           intent(in)    :: test
    integer,                 intent(in)    :: nbatch
    integer,                 intent(in)    :: nsteps
    real(kind=gp),           intent(in)    :: step
    real(kind=gp), optional, intent(in)    :: gamma
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: g(net%nwgts),gprime(net%nwgts),wt(net%nwgts)
    real(kind=gp) :: cst,step_use,alpha

    integer :: n,nc,batch(train%num)

    logical :: track_use,opened,stopfile

    open(unit=unit_stop,file="STOP",err=101)
101 close(unit_stop,status='delete')
    
    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file='nn.track',form='formatted',status='unknown')
    end if

    step_use=step

    alpha=0.0_gp
    wt=0.0_gp
    do nc=1,nsteps

       inquire(file="STOP",exist=stopfile)
       
       if(stopfile) exit
       
       do n=1,train%num
          batch(n)=n
       end do

       call nn_shuffle(batch)

       do n=1,train%num,nbatch

          g=nn_gradient(net,train,cst,batch(n:min(n+nbatch-1,train%num)))
          
          if(present(gamma)) then
             g=g+net%weights*gamma
          end if

          ! * Precondition

          gprime=g*net%pre

          wt=alpha*wt+(1.0_gp-alpha)*step_use*gprime
                    

          ! * Update

!!$          net%weights=net%weights-step_use*gprime

          net%weights=net%weights-wt
          
       end do

       if(track_use.and.(mod(nc-1,100).eq.0)) then
          write (unit_track,*) log10(real(nc,gp)),log10(nn_cost(net,train)),log10(nn_cost(net,test))
          flush(unit_track)
       end if

       step_use=step_use-step/real(nsteps,gp)
          
    end do

    if(track_use) write (unit_track,*)

  end subroutine nn_train_sgd

  subroutine nn_train_adamw(net,train,test,nbatch,nsteps,step,gamma,track)

    type(network),           intent(inout) :: net
    type(dataset),           intent(in)    :: train
    type(dataset),           intent(in)    :: test
    integer,                 intent(in)    :: nbatch
    integer,                 intent(in)    :: nsteps
    real(kind=gp),           intent(in)    :: step
    real(kind=gp),           intent(in)    :: gamma
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: g(net%nwgts),mvec(net%nwgts),vvec(net%nwgts),mvech(net%nwgts),vvech(net%nwgts)

    real(kind=gp) :: cst,step_use,alpha,beta1,beta2,eps

    integer :: n,nc,nt,batch(train%num)

    logical :: track_use,opened,stopfile

    open(unit=unit_stop,file="STOP",err=101)
101 close(unit_stop,status='delete')
    
    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file='nn.track',form='formatted',status='unknown')
    end if

    alpha=0.001_gp
    beta1=0.9_gp
    beta2=0.999_gp
    eps=1e-8_gp

    mvec=0.0_gp
    vvec=0.0_gp
    
    step_use=step

    nt=0
    
    do nc=1,nsteps

       inquire(file="STOP",exist=stopfile)
       
       if(stopfile) exit
       
       do n=1,train%num
          batch(n)=n
       end do

       call nn_shuffle(batch)

       do n=1,train%num,nbatch

          nt=nt+1
          
          g=nn_gradient(net,train,cst,batch(n:min(n+nbatch-1,train%num)))
          
          ! * Update

          mvec=beta1*mvec+(1-beta1)*g
          vvec=beta2*vvec+(1-beta2)*g*g
          
          mvech=mvec/(1-beta1**nt)
          vvech=vvec/(1-beta2**nt)
          
          net%weights=net%weights-step_use*(alpha*mvech/(sqrt(vvech)+eps)+gamma*net%weights)   
          
       end do

       if(track_use.and.(mod(nc-1,100).eq.0)) then
          write (unit_track,*) log10(real(nc,gp)),log10(nn_cost(net,train)),log10(nn_cost(net,test))
          flush(unit_track)
       end if

       step_use=step_use-step/real(nsteps,gp)
          
    end do

    if(track_use) write (unit_track,*)

  end subroutine nn_train_adamw

  subroutine nn_create(net,layers,actfn,outfn)

    type(network), intent(inout) :: net
    
    integer, dimension(:), intent(in) :: layers

    character(len=*), optional, intent(in) :: actfn
    character(len=*), optional, intent(in) :: outfn

    integer :: nl

    net%depth=size(layers)
    
    if(allocated(net%layers)) deallocate(net%layers)
    allocate(net%layers(net%depth))
    if(allocated(net%fnc)) deallocate(net%fnc)
    allocate(net%fnc(net%depth))
    net%layers=layers
    net%inputs=net%layers(1)
    net%outputs=net%layers(net%depth)
    if(present(actfn)) then
       net%fnc=actfn
    else
       net%fnc='tanh'
    end if
    if(present(outfn)) then
       net%fnc(net%depth)=outfn
    else
       net%fnc(net%depth)='linear'
    end if

    if(allocated(net%act)) deallocate(net%act)
    if(allocated(net%out)) deallocate(net%out)
    if(allocated(net%del)) deallocate(net%del)
    if(allocated(net%gra)) deallocate(net%gra)
    if(allocated(net%pre)) deallocate(net%pre)
    if(allocated(net%weights)) deallocate(net%weights)

    net%nnds=net%layers(1)+1 ! ** Add one for offset
    net%nwgts=0
    do nl=2,net%depth
       net%nnds=net%nnds+net%layers(nl)+1
       net%nwgts=net%nwgts+(net%layers(nl-1)+1)*net%layers(nl)
    end do
    
    allocate(net%act(net%nnds))
    net%act=0.0_gp

    allocate(net%out(net%nnds))
    net%out=0.0_gp
    
    allocate(net%del(net%nnds,net%outputs))
    net%del=0.0_gp

    allocate(net%gra(net%nnds))
    net%gra=0.0_gp
    
    allocate(net%weights(net%nwgts))
    net%weights=0.0_gp

    allocate(net%pre(net%nwgts))
    net%pre=1.0_gp

    ! ** Random initialisation
    
!!$    call random_number(net%weights)
!!$
!!$    net%weights=(net%weights-0.5_gp)*2*sqrt(6.0_gp)
!!$
!!$    nw=0
!!$    do nl=2,net%depth
!!$       
!!$       net%weights(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl))=net%weights(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl))/&
!!$            real((net%layers(nl-1)+1*0)+(net%layers(nl)),gp)**(1.0_gp/2.0_gp)
!!$       
!!$       nw=nw+(net%layers(nl-1)+1)*net%layers(nl)
!!$       
!!$    end do

    call nn_initialise(net)
    
    net%cmt='NULL'

    ! ** Input and output scaling

    allocate(net%in_mean(net%inputs))
    allocate(net%in_sd(net%inputs))
    allocate(net%out_mean(net%outputs))
    allocate(net%out_sd(net%outputs))

    net%in_mean=0.0_gp
    net%in_sd=1.0_gp

    net%out_mean=0.0_gp
    net%out_sd=1.0_gp
    
  end subroutine nn_create

  subroutine nn_initialise(net)

    type(network), intent(inout) :: net

    integer :: nw,nl
    
    call random_number(net%weights)

    net%weights=(net%weights-0.5_gp)*2*sqrt(6.0_gp)

    nw=0
    do nl=2,net%depth
       
       net%weights(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl))=net%weights(nw+1:nw+(net%layers(nl-1)+1)*net%layers(nl))/&
            real((net%layers(nl-1)+1*0)+(net%layers(nl)),gp)**(1.0_gp/2.0_gp)
       
       nw=nw+(net%layers(nl-1)+1)*net%layers(nl)
       
    end do
    
  end subroutine nn_initialise

  subroutine nn_summary(net)

    type(network), intent(in) :: net

    write (stderr,'(a,a,a)') ' "',trim(net%cmt(1:80)),' ..."'
    write (stderr,*)
    write (stderr,*) 'number of inputs:  ',net%inputs
    write (stderr,*) 'number of outputs: ',net%outputs
    write (stderr,*) 'number of layers:  ',net%depth
    write (stderr,*) 'number of weights: ',net%nwgts
    write (stderr,*)
    write (stderr,'(a,10a10)') ' network architecture:        ',  net%fnc
    write (stderr,'(a,10i10)   ') '                      ',  net%layers
    write (stderr,*)

  end subroutine nn_summary

  subroutine nn_data(net,in,out,dat)

    type(network), intent(in) :: net

    real(kind=gp), dimension(:,:), intent(in) :: in
    real(kind=gp), dimension(:,:), intent(in) :: out

    type(dataset), intent(out) :: dat

    if(size(in,2).ne.size(out,2))  stop 'nn_data:    num in =/= out'
    if(size(in,1).ne.net%inputs)   stop 'nn_data:   size in =/= net%inputs'
    if(size(out,1).ne.net%outputs) stop 'nn_data:  size out =/= net%outputs'
    
    dat%num=size(in,2)

    if(allocated(dat%in)) deallocate(dat%in)
    if(allocated(dat%out)) deallocate(dat%out)

    allocate(dat%in(net%inputs,dat%num),dat%out(net%outputs,dat%num))

    dat%in=in
    dat%out=out
    
  end subroutine nn_data

  subroutine nn_normalise(net,dat,outonly)

    type(network), intent(inout)    :: net    
    type(dataset), intent(inout)    :: dat
    logical, optional, intent(in)   :: outonly

    logical :: outonly_use,scaleonly

    integer :: n

    if(.not.allocated(net%out_mean)) stop 'nn_normalise: net not initialised'
    if(.not.allocated(dat%in)) stop 'nn_normalise: dat not initialised'

    if(present(outonly)) then
       outonly_use=outonly
    else
       outonly_use=.false.
    end if

    scaleonly=.not.(all(net%in_mean.eq.0.0_gp).and.all(net%in_sd.eq.1.0_gp)&
         .and.all(net%out_mean.eq.0.0_gp).and.all(net%out_sd.eq.1.0_gp))
    
    if(.not.scaleonly) then

       ! Outputs
       
       net%out_mean=sum(dat%out,2)/real(dat%num,gp)

       net%out_sd=epsilon(1.0_gp)
       do n=1,dat%num
          net%out_sd(:)=net%out_sd(:)+(dat%out(:,n)-net%out_mean(:))**2
       end do
       
       net%out_sd=sqrt(net%out_sd/real(dat%num,gp))

       ! Inputs

       if(.not.outonly_use) then

          net%in_mean=sum(dat%in,2)/real(dat%num,gp)

          net%in_sd=epsilon(1.0_gp)
          
          do n=1,dat%num
             net%in_sd(:)=net%in_sd(:)+(dat%in(:,n)-net%in_mean(:))**2
          end do
          
          net%in_sd=sqrt(net%in_sd/real(dat%num,gp))
         
       else

          net%in_mean=0.0_gp
          
          net%in_sd=1.0_gp

       end if
       
    end if

    ! ** Scaling
    
    do n=1,dat%num
       dat%out(:,n)=(dat%out(:,n)-net%out_mean(:))/net%out_sd(:)
    end do
    
    do n=1,dat%num
       dat%in(:,n)=(dat%in(:,n)-net%in_mean(:))/net%in_sd(:)
    end do

  end subroutine nn_normalise
  
  function nn_cost(net,dat) result(cost)

    type(network), intent(inout) :: net
    type(dataset), intent(in)    :: dat

    real(kind=gp) :: cost

    integer :: n
    
    cost=0.0_gp
    do n=1,dat%num
       cost=cost+0.5_gp*norm2(dat%out(:,n)-nn_forward(net,dat%in(:,n)))**2
    end do
    cost=cost/real(dat%num,gp)
    
  end function nn_cost

  subroutine nn_test_gradient(net,dat)

    type(network), intent(inout) :: net
    type(dataset), intent(in)    :: dat


    real(kind=gp), dimension(net%nwgts) :: w0,g,g0

    real(kind=gp) :: step=1e-4_gp,c,cp,cm
    
    integer :: nw

    w0=net%weights
    
    do nw=1,net%nwgts

       net%weights=w0
       net%weights(nw)=net%weights(nw)+step

       cp=nn_cost(net,dat)

       net%weights=w0
       net%weights(nw)=net%weights(nw)-step

       cm=nn_cost(net,dat)

       g(nw)=(cp-cm)/2.0_gp/step
       
    end do

    net%weights=w0

    write (stderr,'(10f10.5)') g
    g0=nn_gradient(net,dat,c)
    write (stderr,'(10f10.5)') g0
    write (stderr,*) norm2(g0-g)

  end subroutine nn_test_gradient
  
  function nn_gradient(net,dat,cost,batch) result(gradient)

    type(network),                   intent(inout) :: net
    type(dataset),                   intent(in)    :: dat
    real(kind=gp),                   intent(out)   :: cost
    integer, optional, dimension(:), intent(in)    :: batch

    type(network) :: net_use
    
    real(kind=gp), dimension(net%nwgts)            :: gradient
    
    real(kind=gp), allocatable, save, dimension(:)   :: deltay,grad
    real(kind=gp), allocatable, save, dimension(:,:) :: back
    
    integer, allocatable, dimension(:) :: batch_use
    
    integer :: n

    if(present(batch)) then
       allocate(batch_use(size(batch)))
       batch_use=batch
    else
       allocate(batch_use(dat%num))
       do n=1,dat%num
          batch_use(n)=n
       end do
    end if

    if(.not.allocated(deltay)) &
         allocate(deltay(net%outputs),grad(net%nwgts),back(net%nwgts,net%outputs))
    
    cost=0.0_gp
    grad=0.0_gp

    !$omp parallel do private(deltay,back,net_use) reduction(+:cost,grad) schedule(dynamic)
    do n=1,size(batch_use)

       net_use=net
       
       ! ** Forward
       
       deltay=dat%out(:,batch_use(n))-nn_forward(net_use,dat%in(:,batch_use(n)))
       
       cost=cost+0.5_gp*norm2(deltay)**2
       
       ! ** Backward

       back=nn_backward(net_use)
       
       if(gp.eq.sp) then
          call sgemv('N',net_use%nwgts,net_use%outputs,1.0_gp,back(1,1),net_use%nwgts,deltay(1),1,1.0_gp,grad(1),1)
       else
          call dgemv('N',net_use%nwgts,net_use%outputs,1.0_gp,back(1,1),net_use%nwgts,deltay(1),1,1.0_gp,grad(1),1)          
       end if
       
    end do
    !$omp end parallel do

    
    cost=cost/real(size(batch_use),gp)
    
    gradient=grad/real(size(batch_use),gp)
    
  end function nn_gradient

  function nn_lmstep(net,dat,lambda,cost,batch) result(step)

    type(network),                   intent(inout) :: net
    type(dataset),                   intent(in)    :: dat
    real(kind=gp),                   intent(in)    :: lambda
    real(kind=gp),                   intent(out)   :: cost
    integer, optional, dimension(:), intent(in)    :: batch


    real(kind=gp), dimension(net%nwgts)            :: step
    real(kind=gp), dimension(net%nwgts)            :: bvec

    type(network) :: net_use
    
    real(kind=gp), allocatable, dimension(:)      :: deltay
    real(kind=gp), allocatable, dimension(:,:)    :: jj,jac
    
    integer, allocatable, dimension(:) :: batch_use
    
    integer :: n
    
    if(present(batch)) then
       allocate(batch_use(size(batch)),jac(net%nwgts,net%outputs))
       allocate(deltay(net%outputs),jj(net%nwgts,net%nwgts))
       batch_use=batch
    else
       allocate(batch_use(dat%num),jac(net%nwgts,net%outputs),deltay(net%outputs),jj(net%nwgts,net%nwgts))
       do n=1,dat%num
          batch_use(n)=n
       end do
    end if
    
    cost=0.0_gp

    jj=0.0_gp

    bvec=0.0_gp
    
    !$omp parallel do private(net_use,jac,deltay) reduction(+:cost,jj,bvec) schedule(dynamic)
    do n=1,size(batch_use)
       
       net_use=net
       
       ! ** Forward
       
       deltay(:)=dat%out(:,batch_use(n))-nn_forward(net_use,dat%in(:,batch_use(n)))

       cost=cost+0.5_gp*norm2(deltay(:))**2
       
       ! ** Backward

       jac(:,:)=nn_backward(net_use)
   
       ! ** Accumulate jj, through outer product

       if(gp.eq.sp) then
          call ssyrk('U','N',net%nwgts,net_use%outputs,1.0_gp,jac(1,1),net%nwgts,1.0_gp,jj(1,1),net%nwgts)
       else
          call dsyrk('U','N',net%nwgts,net_use%outputs,1.0_gp,jac(1,1),net%nwgts,1.0_gp,jj(1,1),net%nwgts)
       end if
       
       ! ** Accumulate bvec
       
       bvec=bvec+matmul(jac(:,:),deltay(:))
       
    end do
    !$omp end parallel do
    
    cost=cost/real(size(batch_use),gp) ! ** CHECK THIS NORMALISATION
    
    do n=1,net%nwgts
       jj(n,n)=jj(n,n)+lambda*max(1.0_gp,jj(n,n))
    end do
    
    step=-nn_solve(jj,bvec)
    
  end function nn_lmstep

 function nn_galmstep(net,dat,lambda,cost,batch) result(step)

    type(network),                   intent(inout) :: net
    type(dataset),                   intent(in)    :: dat
    real(kind=gp),                   intent(in)    :: lambda
    real(kind=gp),                   intent(out)   :: cost
    integer, optional, dimension(:), intent(in)    :: batch


    real(kind=gp), dimension(net%nwgts)            :: step
    real(kind=gp), dimension(net%nwgts)            :: bvec

    type(network) :: net_use,net_plus,net_zero,net_minus

    real(kind=gp), dimension(net%nwgts) :: accel
    
    real(kind=gp), allocatable, dimension(:,:)      :: deltay,jj,fvv
    real(kind=gp), allocatable, dimension(:,:,:)    :: jac
    
    integer, allocatable, dimension(:) :: batch_use
    
    real(kind=gp) :: hstep=0.1_gp
    
    integer :: n
    
    if(present(batch)) then
       allocate(batch_use(size(batch)),jac(net%nwgts,net%outputs,size(batch)))
       allocate(deltay(net%outputs,size(batch)),jj(net%nwgts,net%nwgts),fvv(net%outputs,size(batch)))
       batch_use=batch
    else
       allocate(batch_use(dat%num),jac(net%nwgts,net%outputs,dat%num),deltay(net%outputs,dat%num))
       allocate(jj(net%nwgts,net%nwgts),fvv(net%outputs,dat%num))
       do n=1,dat%num
          batch_use(n)=n
       end do
    end if

    cost=0.0_gp

    jj=0.0_gp

    !$omp parallel do private(net_use) reduction(+:cost,jj) schedule(dynamic)
    do n=1,size(batch_use)

       net_use=net
       
       ! ** Forward
       
       deltay(:,n)=dat%out(:,batch_use(n))-nn_forward(net_use,dat%in(:,batch_use(n)))

       cost=cost+0.5_gp*norm2(deltay(:,n))**2
       
       ! ** Backward

       jac(:,:,n)=nn_backward(net_use)
   
       ! ** Accumulate jj, through outer product

       if(gp.eq.sp) then
          call ssyrk('U','N',net%nwgts,net_use%outputs,1.0_gp,jac(1,1,n),net%nwgts,1.0_gp,jj(1,1),net%nwgts)
       else
          call dsyrk('U','N',net%nwgts,net_use%outputs,1.0_gp,jac(1,1,n),net%nwgts,1.0_gp,jj(1,1),net%nwgts)   
       end if
       
    end do
    !$omp end parallel do

    cost=cost/real(size(batch_use),gp) ! ** CHECK THIS NORMALISATION
    
    do n=1,net%nwgts
       jj(n,n)=jj(n,n)+lambda*max(1.0_gp,jj(n,n))
    end do

    bvec=matmul(reshape(jac,(/net%nwgts,net%outputs*dat%num/)),reshape(deltay,(/net%outputs*dat%num/)))

    step=-nn_solve(jj,bvec)

    ! ** Compute geodesic acceleration

    fvv=0.0_gp
    
    !$omp parallel do private(net_plus,net_zero,net_minus) schedule(dynamic)
    do n=1,size(batch_use)

       net_plus=net
       net_plus%weights=net_plus%weights+hstep*step

       net_zero=net
       
       net_minus=net
       net_minus%weights=net_minus%weights-hstep*step

       ! ** Forward
       
       fvv(:,n)=nn_forward(net_plus,dat%in(:,batch_use(n)))
       fvv(:,n)=fvv(:,n)-2*nn_forward(net_zero,dat%in(:,batch_use(n)))
       fvv(:,n)=fvv(:,n)+nn_forward(net_minus,dat%in(:,batch_use(n)))
       
    end do
    !$omp end parallel do

    fvv=fvv/hstep**2

    bvec=matmul(reshape(jac,(/net%nwgts,net%outputs*dat%num/)),reshape(fvv,(/net%outputs*dat%num/)))

    accel=-nn_solve(jj,bvec)

    if(norm2(accel)/norm2(step).lt.0.5_gp) step=step-accel/2
    
  end function nn_galmstep

  function nn_forward(net,x,scale) result(y)

    type(network), intent(inout)                     :: net
    real(kind=gp), dimension(net%inputs), intent(in) :: x
    real(kind=gp), dimension(net%outputs)            :: y
    logical, optional                                :: scale
    
    integer :: nnl,mnl,mnlpo,nl,nnds,nwgts

    logical :: scale_use

    if(present(scale)) then
       scale_use=scale
    else
       scale_use=.false.
    end if
    
    ! ** Fill the input nodes

    if(scale_use) then
       net%out(1:net%inputs)=(x(1:net%inputs)-net%in_mean(1:net%inputs))/net%in_sd(1:net%inputs)
    else
       net%out(1:net%inputs)=x(1:net%inputs)
    end if
    net%out(net%inputs+1)=1.0_gp

    nnds=net%inputs+1

    ! ** Work though the layers of the network
    
    nwgts=0
    
    do nl=2,net%depth

       mnl=net%layers(nl-1)
       nnl=net%layers(nl)
       mnlpo=mnl+1
       
       if(gp.eq.sp) then
          call sgemv('N',nnl,mnlpo,1.0_gp,net%weights(nwgts+1),nnl,net%out(nnds-mnl),1,0.0_gp,net%act(nnds+1),1)
       else
          call dgemv('N',nnl,mnlpo,1.0_gp,net%weights(nwgts+1),nnl,net%out(nnds-mnl),1,0.0_gp,net%act(nnds+1),1)
       end if
       
       ! * Activation
       
       net%out(nnds+1:nnds+nnl)=nn_activation(net%act(nnds+1:nnds+nnl),net%fnc(nl))
       
       net%gra(nnds+1:nnds+nnl)=nn_activationp(net%act(nnds+1:nnds+nnl),net%fnc(nl))
       
       ! * Add dummy node for offset
       
       net%out(nnds+nnl+1)=1.0_gp

       net%gra(nnds+nnl+1)=0.0_gp

       ! * Increment weight and node counters
       
       nnds=nnds+nnl+1

       nwgts=nwgts+nnl*(mnl+1)
       
    end do
    
    ! ** Set the output of the network
    
    y(1:net%outputs)=net%out(nnds-net%layers(net%depth):nnds-1)

    if(scale_use) y(1:net%outputs)=y(1:net%outputs)*net%out_sd(1:net%outputs)+net%out_mean(1:net%outputs)
    
  end function nn_forward

  function nn_backward(net) result(dydw)

    type(network), intent(inout) :: net

    real(kind=gp), dimension(net%nwgts,net%outputs) :: dydw
    
    integer :: nnl,nnlpo,pnl,mnl,nl,nnds,nwgts,n,no
    
    nnds=net%nnds
    nwgts=net%nwgts
    
    ! * start off

    dydw=0
    
    nnl=net%layers(net%depth)
    mnl=net%layers(net%depth-1)


    net%del=0
   
    
    do no=1,net%outputs

       do n=1,net%outputs
          if(n.eq.no) net%del(nnds-nnl-1+n,no)=1.0_gp
       end do
       
       net%del(nnds-nnl:nnds,no)=-net%gra(nnds-nnl:nnds)*net%del(nnds-nnl:nnds,no)
       
    end do

    do no=1,net%outputs

       ! * outer product

       if(gp.eq.sp) then
          call sger(nnl,mnl+1,1.0_gp,net%del(nnds-nnl,no),1,net%out(nnds-mnl-nnl-1),1,dydw(nwgts-(mnl+1)*nnl+1,no),nnl) !! Use DGEMM to get out of no loop 
       else
          call dger(nnl,mnl+1,1.0_gp,net%del(nnds-nnl,no),1,net%out(nnds-mnl-nnl-1),1,dydw(nwgts-(mnl+1)*nnl+1,no),nnl) !! Use DGEMM to get out of no loop
       end if
       
    end do
    
    net%pre(:nwgts)=1
    
    nnds=nnds-nnl-1
    nwgts=nwgts-(mnl+1)*nnl
    
    ! * work backwards
    
    do nl=net%depth-1,2,-1

       pnl=net%layers(nl+1)
       nnl=net%layers(nl)
       mnl=net%layers(nl-1)
       nnlpo=nnl+1

       do no=1,net%outputs
          
          if(gp.eq.sp) then
             call sgemv('T',pnl,nnlpo,1.0_gp,net%weights(nwgts+1),pnl,net%del(nnds+1,no),1,0.0_gp,net%del(nnds-nnl,no),1) !! Use DGEMM ??
          else
             call dgemv('T',pnl,nnlpo,1.0_gp,net%weights(nwgts+1),pnl,net%del(nnds+1,no),1,0.0_gp,net%del(nnds-nnl,no),1) !! Use DGEMM ??
          end if
          
          net%del(nnds-nnl:nnds,no)=net%gra(nnds-nnl:nnds)*net%del(nnds-nnl:nnds,no)
          
          ! * outer product

          if(gp.eq.sp) then
             call sger(nnl,mnl+1,1.0_gp,net%del(nnds-nnl,no),1,net%out(nnds-mnl-nnl-1),1,dydw(nwgts-(mnl+1)*nnl+1,no),nnl) !! Use DGEMM to get out of no loop
          else
             call dger(nnl,mnl+1,1.0_gp,net%del(nnds-nnl,no),1,net%out(nnds-mnl-nnl-1),1,dydw(nwgts-(mnl+1)*nnl+1,no),nnl) !! Use DGEMM to get out of no loop
          end if
          
       end do
       
!!$       net%pre(:nwgts)=net%pre(:nwgts)*1.0_gp
    
       nnds=nnds-nnl-1
       nwgts=nwgts-(mnl+1)*nnl

    end do
    
  end function nn_backward

  subroutine nn_save(net,name)
    
    type(network), intent(in) :: net

    character(len=*) :: name

    open(unit=unit_nn,file=name,form='formatted',status='unknown')
    
    write (unit_nn,'(5(1x,i0))')    net%depth,net%inputs,net%outputs,net%nwgts,net%nnds
    write (unit_nn,'(*(1x,i0))')    net%layers
    write (unit_nn,'(*(1x,f0.14))') net%weights
    write (unit_nn,'(*(1x,f0.14))') net%pre
    write (unit_nn,'(*(1x,f0.14))') net%act
    write (unit_nn,'(*(1x,f0.14))') net%out
    write (unit_nn,'(*(1x,f0.14))') net%del
    write (unit_nn,'(*(1x,f0.14))') net%gra
    write (unit_nn,'(*(1x,a))')     net%fnc
    write (unit_nn,'(*(1x,f0.14))') net%in_mean
    write (unit_nn,'(*(1x,f0.14))') net%in_sd
    write (unit_nn,'(*(1x,f0.14))') net%out_mean
    write (unit_nn,'(*(1x,f0.14))') net%out_sd
    write (unit_nn,'(a)')           net%cmt

    close(unit_nn)
    
  end subroutine nn_save

  subroutine nn_load(net,name)

    type(network), intent(out) :: net
    
    character(len=*) :: name

    open(unit=unit_nn,file=name,form='formatted',status='old',err=55)
    
    read (unit_nn,*) net%depth,net%inputs,net%outputs,net%nwgts,net%nnds

    if(allocated(net%layers)) deallocate(net%layers)
    allocate(net%layers(net%depth))
    if(allocated(net%fnc)) deallocate(net%fnc)
    allocate(net%fnc(net%depth))
    if(allocated(net%act)) deallocate(net%act)
    allocate(net%act(net%nnds))
    if(allocated(net%out)) deallocate(net%out)
    allocate(net%out(net%nnds))
    if(allocated(net%del)) deallocate(net%del)
    allocate(net%del(net%nnds,net%outputs))
    if(allocated(net%gra)) deallocate(net%gra)
    allocate(net%gra(net%nnds))
    if(allocated(net%pre)) deallocate(net%pre)
    allocate(net%pre(net%nwgts))
    if(allocated(net%weights)) deallocate(net%weights)
    allocate(net%weights(net%nwgts))
    if(allocated(net%in_mean)) deallocate(net%in_mean)
    allocate(net%in_mean(net%inputs))
    if(allocated(net%in_sd)) deallocate(net%in_sd)
    allocate(net%in_sd(net%inputs))
    if(allocated(net%out_mean)) deallocate(net%out_mean)
    allocate(net%out_mean(net%outputs))
    if(allocated(net%out_sd)) deallocate(net%out_sd)
    allocate(net%out_sd(net%outputs))
    
    read (unit_nn,*) net%layers
    read (unit_nn,*) net%weights
    read (unit_nn,*) net%pre
    read (unit_nn,*) net%act
    read (unit_nn,*) net%out
    read (unit_nn,*) net%del
    read (unit_nn,*) net%gra
    read (unit_nn,*) net%fnc
    read (unit_nn,*) net%in_mean
    read (unit_nn,*) net%in_sd
    read (unit_nn,*) net%out_mean
    read (unit_nn,*) net%out_sd
    read (unit_nn,'(a)') net%cmt

    net%cmt=adjustl(net%cmt)

    close(unit_nn)

    return

55  stop 'nn_load: loadfile cannot be opened'
    
  end subroutine nn_load
  
  elemental real(kind=gp) function nn_activation(x,f) result(y)

    real(kind=gp),    intent(in) :: x
    character(len=*), intent(in) :: f

    select case(f)
    case('sigmoid')
       y=1.0_gp/(1.0_gp+exp(-x))
    case('tanh')
       y=tanh(x)
       y=2*y+1e-6_gp*x
    case('relu')
       y=max(0.0_gp,x)
    case('qrelu')
       y=max(0.0_gp,x)**2
    case('lrelu')
       y=max(0.1_gp*x,x)
    case('linear')
       y=x
    case('quadratic')
       y=x**2/2
    case('swish','silu')
       y=x/(1.0_gp+exp(-x))
    case('lj')
       y=1/x**12-1/x**6
    case('softplus')
       y=log(1+exp(x))
    case default
       y=huge(1.0_gp)
    end select
    
  end function nn_activation

  elemental real(kind=gp) function nn_activationp(x,f) result(y)

    real(kind=gp),    intent(in) :: x
    character(len=*), intent(in) :: f

    select case(f)
    case('sigmoid')
       y=exp(-x)/(1.0_gp+exp(-x))**2
    case('tanh')
       y=1.0_gp-tanh(x)**2
       y=2*y+1e-6_gp
    case('relu')
       if(x.gt.0.0_gp) then
          y=1.0_gp
       else
          y=0.0_gp
       end if
    case('qrelu')
       if(x.gt.0.0_gp) then
          y=2*x
       else
          y=0.0_gp
       end if
    case('lrelu')
       if(x.gt.0.0_gp) then
          y=1.0_gp
       else
          y=0.1_gp
       end if
    case('linear')
       y=1.0_gp
    case('quadratic')
       y=x
    case('swish','silu')
       y=x/(1.0_gp+exp(-x))+(1-x/(1.0_gp+exp(-x)))/(1+exp(-x))
    case('lj')
       y=-12/x**13+6/x**7
    case('softplus')
       y=1/(1+exp(-x))
    case default
       y=huge(1.0_gp)
    end select
    
  end function nn_activationp

  subroutine nn_shuffle(list)

    integer, dimension(:), intent(out) :: list

    real(kind=gp) :: rnd
    
    integer :: i,j,itmp

    do i=1,size(list)
       list(i)=i
    end do

    do i=size(list),2,-1
       call random_number(rnd)
       j=int(rnd*i)+1
       itmp=list(j)
       list(j)=list(i)
       list(i)=itmp
    end do

  end subroutine nn_shuffle

  function nn_solve(A,B) result(BB)

    real(kind=gp), dimension(:,:), intent(in) :: A
    real(kind=gp), dimension(:),   intent(in) :: B

!!$    real(kind=gp), dimension(size(A,1),size(A,2)) :: AA
    real(kind=gp), save, allocatable, dimension(:,:) :: AA
    real(kind=gp), dimension(size(B,1)) :: BB
    
    integer,       dimension(size(A,1)) :: ipiv

    real(kind=gp), save, allocatable, dimension(:) :: work 
!!$    real(kind=gp), dimension(size(B)**2/2) :: work 

    integer :: n, info, lwork

    external dsysv,ssysv

    if(.not.allocated(work)) allocate(work(size(B)**2/2))
    
    AA=A
    BB=B

    info=0
    lwork=-1
    
    n=size(A,1)
    
    if(gp.eq.sp) then
       call ssysv('U',n,1,AA,n,ipiv,BB,n,work,lwork,info)       
    else
       call dsysv('U',n,1,AA,n,ipiv,BB,n,work,lwork,info)
    end if
    
    lwork = min(size(work),int(work(1)))
    if(info.ne.0) stop 'nn_solve: init failed'
    
    if(gp.eq.sp) then
       call ssysv('U',n,1,AA,n,ipiv,BB,n,work,lwork,info)
    else
       call dsysv('U',n,1,AA,n,ipiv,BB,n,work,lwork,info)
    end if
    
    if(info.ne.0) stop 'nn_solve: solve failed'

  end function nn_solve
   
end module nn
