!==================================================================================!
!                                    forge                                         !
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
! This program generates ddps by training neural networks                          !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
program forge

  use constants
  use omp_lib
  use data
  use nn
#ifdef IFORT
  use ifport, only : getpid
#endif

  implicit none

  integer, parameter :: unit_track=44
  integer, parameter :: unit_delta=46
  integer, parameter :: unit_agr=47
  integer, parameter :: unit_stop=55
  integer, parameter :: unit_seed=56

  integer, parameter :: max_layers=100

  real(kind=gp) :: eta,lambda_min,pval,kT,wmin

  integer :: nchunks,nsteps,nbatch,nnodes(max_layers),nlayers,nomp,nestop,npot,maxpot

  integer, allocatable, dimension(:) :: seed, netarch

  type(network) :: net

  character(len=240) :: trainfile,validfile,testfile,netfile,seedname

  type(dataset) :: training,validation,testing

  type(metadata) :: metatrain,metavalid,metatest

  real(kind=gp) :: thresh

  character(len=10) :: act,opt
  character(len=20) :: ctemp

  character(len=1000) :: randomseed

  logical :: load=.false.,trk=.false.,quiet=.false.,earlystop=.true.,pdf=.false.,plot=.true.,filestop=.false.
  logical :: loadseed

  real(kind=gp) :: total_time,start_time,end_time,load_time=0.0_gp,train_time=0.0_gp,test_time=0.0_gp
  real(kind=gp) :: jac_time=0.0_gp,bvec_time=0.0_gp,solve_time=0.0_gp,cost_time=0.0_gp,dsyrk_time=0.0_gp

  ! ************************************************************  

  call cpu_time(start_time)
  
  inquire (file="randomseed",exist=loadseed)

  if(loadseed) then
     open (unit=unit_seed,file='randomseed',form='formatted',status='old')
     read(unit_seed,'(a)') randomseed
     call init_pseudorandom(trim(randomseed))
  else
     call init_pseudorandom()
  end if

  call get_arguments()

  call omp_set_num_threads(nomp) 

  call banner()

  ! ************************************************************

  ! ** Load data

  if(.not.quiet) write (stderr,'(a)') ' loading training data from: '//trim(trainfile)
  call data_read(trainfile,training,metatrain)
  if(.not.quiet) write (stderr,'(a,i10)') ' structures found:   ',metatrain%nstruct
  if(.not.quiet) write (stderr,'(a,i10)') ' centres found:      ',training%num
  if(.not.quiet) write (stderr,*)

  if(.not.quiet) write (stderr,'(a)') ' loading validation data from: '//trim(validfile)
  call data_read(validfile,validation,metavalid)
  if(.not.quiet) write (stderr,'(a,i10)') ' structures found:   ',metavalid%nstruct
  if(.not.quiet) write (stderr,'(a,i10)') ' centres found:      ',validation%num
  if(.not.quiet) write (stderr,*)

  if(.not.data_compatible(training,metatrain,validation,metavalid)) stop 'validation data not consistent with training data'

  if(.not.quiet) write (stderr,'(a)') ' loading testing data from: '//trim(testfile)
  call data_read(testfile,testing,metatest)
  if(.not.quiet) write (stderr,'(a,i10)') ' structures found:   ',metatest%nstruct
  if(.not.quiet) write (stderr,'(a,i10)') ' centres found:      ',testing%num
  if(.not.quiet) write (stderr,*)

  if(.not.data_compatible(training,metatrain,testing,metatest)) stop 'testing data not consistent with training data'

  call cpu_time(load_time)

  load_time=load_time-start_time
  
  ! ** Create network

  if(load) then
     if(.not.quiet) write (stderr,'(a,/)') ' loading data derived potential from '//trim(netfile)
     call nn_load(net,netfile)
  else

     if(.not.quiet) write (stderr,'(a,/)') ' creating new data derived potential'

     allocate(netarch(2+nlayers))

     netarch(1)=training%idim
     netarch(2:1+nlayers)=nnodes(1:nlayers)
     netarch(2+nlayers)=training%odim

     call nn_create(net,netarch,actfn=act)

     ! ** Normalise based on the training set

     call nn_normalise(net,training,outonly=.false.)
     call nn_normalise(net,validation,outonly=.false.)
     call nn_normalise(net,testing,outonly=.false.)

     ! ** Add the comment to document the features and composition space

     net%cmt=trim(metatrain%cmnt)

  end if

  if(.not.quiet) call nn_summary(net)

  if(.false.) call forge_test_gradient(net,testing,metatest)

  ! ** Train network

  do npot=1,maxpot

     call nn_initialise(net)
     
     select case(opt)
     case('lin')
        call forge_train_lin()
     case('lm')    
        call forge_train_lm(thresh,track=trk)
     case('tpsd')
        call forge_train_tpsd(thresh,track=trk)
     case('sgdm')
        if(nbatch.le.0) nbatch=metatrain%nstruct/10
        call forge_train_sgdm(nbatch,nsteps,eta,track=trk)
     case('adamw')
        if(nbatch.le.0) nbatch=metatrain%nstruct/10
        call forge_train_adamw(nbatch,nsteps,eta,track=trk)
     case default
        stop 'optimiser not known'
     end select

     ! ** Save parameters

     write(ctemp,*) npot
     
     if(.not.quiet) write (stderr,*) 'saving data derived potential to '//trim(seedname)//'-'//trim(adjustl(ctemp))//'.ddp'

     call nn_save(net,trim(seedname)//'-'//trim(adjustl(ctemp))//'.ddp')

     ! ************************************************************

     ! ** Test against training data

     call forge_test(' training:   ',net,training,metatrain)
     call forge_test(' validation: ',net,validation,metavalid)
     call forge_test(' testing:    ',net,testing,metatest)

     write (stderr,*)

  end do

  call cpu_time(end_time)

  total_time=end_time-start_time
  
  write (stderr,*)
  write (stderr,'(a)') '------------------------------------'
  write (stderr,'(a,f13.4,a,f7.2,a)') 'total   : ',total_time,' sec',100.0_gp,' %'

  write (stderr,'(a,f13.4,a,f7.2,a)') ' -load    : ',load_time,' sec',100*(load_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -train   : ',train_time,' sec',100*(train_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') '  -jac    : ',jac_time,' sec',100*(jac_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') '   -dsyrk : ',dsyrk_time,' sec',100*(dsyrk_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') '  -cost   : ',cost_time,' sec',100*(cost_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') '  -bvec   : ',bvec_time,' sec',100*(bvec_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') '  -solve  : ',solve_time,' sec',100*(solve_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -test    : ',test_time,' sec',100*(test_time)/total_time,' %'
  
contains

  subroutine forge_train_lin()

    real(kind=gp) :: M(net%nwgts,net%nwgts),X(net%nwgts,metatrain%nstruct)
    real(kind=gp) :: Y(metatrain%nstruct),W(metatrain%nstruct)
    real(kind=gp) :: beta(net%nwgts)

    real(kind=gp) :: Emin,Emean

    integer :: nw,ns

    if((nlayers.ne.0).or.(act.ne.'linear')) stop 'forge_train_lin: network architecture not linear'

    if(any(metatrain%ncentres.gt.1)) stop 'forge_train_lin: mean features only'

    X(1:net%nwgts-1,:)=training%in
    X(net%nwgts,:)=1.0_gp

    Y(:)=training%out(1,:)

    Emin=minval(Y)

    Emean=sum(Y)/size(Y)

    W=exp(-(Y-Emin)/kT)

    do ns=1,metatrain%nstruct
       X(:,ns)=sqrt(W(ns))*X(:,ns)
    end do

    Y=sqrt(W)*Y

!!$    M=matmul(X,transpose(X))

   if(gp.eq.sp) then
      call sgemm('N','T',net%nwgts,net%nwgts,metatrain%nstruct,1.0_gp,X,net%nwgts,X,net%nwgts,1.0_gp,M,net%nwgts)
   else
      call dgemm('N','T',net%nwgts,net%nwgts,metatrain%nstruct,1.0_gp,X,net%nwgts,X,net%nwgts,1.0_gp,M,net%nwgts)
    end if
    
    do nw=1,net%nwgts
       M(nw,nw)=M(nw,nw)+lambda_min
    end do

    beta=nn_solve(M,matmul(X,Y))

    net%weights=beta

  end subroutine forge_train_lin

  subroutine forge_train_tpsd(thresh,track)

    real(kind=gp),     intent(in) :: thresh
    logical, optional, intent(in) :: track

    integer :: nc

    real(kind=gp) :: cst,g(net%nwgts),step
    real(kind=gp) :: g0(net%nwgts),w0(net%nwgts)

    logical :: stopfile,track_use,opened

    if(filestop) then
       open(unit=unit_stop,file="STOP",err=100)
100    close(unit_stop,status='delete')
    end if

    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file=trim(seedname)//'.track',form='formatted',status='unknown')
    end if

    step=1e-5_gp
    g=1.0_gp
    nc=0
    do while(norm2(g).gt.thresh)
       nc=nc+1

       if(filestop) then
          inquire(file="STOP",exist=stopfile)
          if(stopfile) return
       end if

       g=forge_gradient(net,training,metatrain,cst)

       if(nc.gt.1) step=abs(dot_product(net%weights-w0,g-g0)/dot_product(g-g0,g-g0))

       g0=g
       w0=net%weights

       net%weights=net%weights-step*g

       if(mod(nc-1,100).eq.0) then
          write (unit_track,*) log10(real(nc,gp)),log10(cst*2)/2,log10(forge_cost(net,validation,metavalid)*2)/2
          flush(unit_track)
       end if

    end do

  end subroutine forge_train_tpsd

  subroutine forge_train_lm(thresh,track)

    real(kind=gp),     intent(in) :: thresh
    logical, optional, intent(in) :: track

    integer :: ncount,ncb

    real(kind=gp) :: lambda,cst,cst0,cstv,step(net%nwgts),best_weights(net%nwgts),best_cstv,emin

    real(kind=gp) :: start_time,end_time
    
    logical :: stopfile,track_use,opened

    call cpu_time(start_time)
    
    if(filestop) then
       open(unit=unit_stop,file="STOP",err=100)
100    close(unit_stop,status='delete')
    end if

    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file=trim(seedname)//'.track',form='formatted',status='unknown')
    end if

    emin=minval(training%out)

    step=1.0_gp
    lambda=1e5_gp
    best_cstv=huge(1.0_gp)
    ncount=0
    ncb=0
    do while((norm2(step).gt.thresh).or.(ncount.lt.10))

       if(filestop) then
          inquire(file="STOP",exist=stopfile)
          if(stopfile) exit
       end if

       if((ncb.gt.nestop).and.earlystop) exit

       if(ncount.gt.nsteps) exit

       step=forge_lmstep(net,training,metatrain,lambda,cst0,emin)

       cstv=forge_cost(net,validation,metavalid,emin) !!!! Pass in weight ...

       net%weights=net%weights+step

       cst=forge_cost(net,training,metatrain,emin)

       if(cst.gt.cst0) then
          net%weights=net%weights-step
          lambda=lambda*2.0_gp
       else
          lambda=max(lambda/5.0_gp,lambda_min)
          ncount=ncount+1
          ncb=ncb+1
          if(track_use) then
             write (unit_track,'(6f10.5)') &
                  log10(real(ncount,gp)),log10(cst0*2)/2,log10(cstv*2)/2,log10(norm2(step)),log10(lambda),log10(norm2(net%weights))
             flush(unit_track)
          end if
       end if

       if(cstv.lt.best_cstv) then
          best_cstv=cstv
          best_weights=net%weights
          ncb=0
       end if

    end do

    if(earlystop) net%weights=best_weights

    call cpu_time(end_time)

    train_time=train_time+end_time-start_time
    
  end subroutine forge_train_lm

  subroutine forge_train_sgdm(nbatch,nsteps,step,track)

    integer,                 intent(in)    :: nbatch
    integer,                 intent(in)    :: nsteps
    real(kind=gp),           intent(in)    :: step
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: g(net%nwgts),v(net%nwgts)
    real(kind=gp) :: cst,step_use,beta,alpha

    integer :: n,nc,batch(metatrain%nstruct)

    logical :: track_use,opened,stopfile

    if(filestop) then
       open(unit=unit_stop,file="STOP",err=101)
101    close(unit_stop,status='delete')
    end if

    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file=trim(seedname)//'.track',form='formatted',status='unknown')
    end if

    step_use=step
    v=0.0_gp
    alpha=0.00001_gp
    beta=0.999_gp

    do nc=1,nsteps

       if(filestop) then
          inquire(file="STOP",exist=stopfile)
          if(stopfile) return
       end if

       do n=1,metatrain%nstruct
          batch(n)=n
       end do

       call nn_shuffle(batch)

       do n=1,metatrain%nstruct,nbatch

          g=forge_gradient(net,training,metatrain,cst,batch(n:min(n+nbatch-1,metatrain%nstruct)))          

          ! * Update

          v=beta*v+(1-beta)*g

          net%weights=net%weights-step_use*alpha*v

       end do

       if(track_use.and.(mod(nc-1,100).eq.0)) then
          write (unit_track,*) log10(real(nc,gp)),log10(forge_cost(net,training,metatrain)*2)/2,&
               log10(forge_cost(net,validation,metavalid)*2)/2
          flush(unit_track)
       end if

       step_use=step_use-step/real(nsteps,gp)

    end do

    if(track_use) write (unit_track,*)

  end subroutine forge_train_sgdm

  subroutine forge_train_adamw(nbatch,nsteps,step,track)

    integer,                 intent(in)    :: nbatch
    integer,                 intent(in)    :: nsteps
    real(kind=gp),           intent(in)    :: step
    logical,       optional, intent(in)    :: track

    real(kind=gp) :: g(net%nwgts),mvec(net%nwgts),vvec(net%nwgts),mvech(net%nwgts),vvech(net%nwgts)
    real(kind=gp) :: cst

    real(kind=gp) :: alpha,beta1,beta2,eps,lambda,step_use

    integer :: n,nc,nt,batch(metatrain%nstruct)

    logical :: track_use,opened,stopfile

    if(filestop) then
       open(unit=unit_stop,file="STOP",err=101)
101    close(unit_stop,status='delete')
    end if

    if(present(track)) then
       track_use=track
    else
       track_use=.false.
    end if

    if(track_use) then
       inquire( unit=unit_track, opened=opened ) 
       if(.not.opened) open(unit=unit_track,file=trim(seedname)//'.track',form='formatted',status='unknown')
    end if


    alpha=0.001_gp
    beta1=0.9_gp
    beta2=0.999_gp
    eps=1e-8_gp
    lambda=lambda_min

    mvec=0.0_gp
    vvec=0.0_gp

    step_use=step

    nt=0

    do nc=1,nsteps

       if(filestop) then
          inquire(file="STOP",exist=stopfile)
          if(stopfile) return
       end if

       do n=1,metatrain%nstruct
          batch(n)=n
       end do

       call nn_shuffle(batch)

       do n=1,metatrain%nstruct,nbatch

          nt=nt+1

          g=forge_gradient(net,training,metatrain,cst,batch(n:min(n+nbatch-1,metatrain%nstruct)))

          ! * Update

          mvec=beta1*mvec+(1-beta1)*g
          vvec=beta2*vvec+(1-beta2)*g*g

          mvech=mvec/(1-beta1**nt)
          vvech=vvec/(1-beta2**nt)

          net%weights=net%weights-step_use*(alpha*mvech/(sqrt(vvech)+eps)+lambda*net%weights)   

       end do

       if(track_use.and.(mod(nc-1,100).eq.0)) then
          write (unit_track,*) log10(real(nc,gp)),log10(forge_cost(net,training,metatrain)*2)/2,&
               log10(forge_cost(net,validation,metavalid)*2)/2
          flush(unit_track)
       end if

       step_use=step_use-step/real(nsteps,gp)

    end do

    if(track_use) write (unit_track,*)

  end subroutine forge_train_adamw

  function forge_cost(net,data,meta,emin) result(cost)

    type(network),  intent(inout) :: net
    type(dataset),  intent(in)   :: data
    type(metadata), intent(in)   :: meta
    real(kind=gp),  optional, intent(in) :: emin

    real(kind=gp) :: cost,dcost

    type(network), save :: net_use

    real(kind=gp) :: din(net%inputs)

    real(kind=gp) :: es(1),edata,start_time,end_time

    integer :: ns,nc

    call cpu_time(start_time)
    
    cost=0.0_gp

    !$omp parallel do private(es,din,nc,net_use,dcost,edata) reduction(+:cost) schedule(dynamic)
    do ns=1,meta%nstruct

       net_use=net

       es(1:1)=0.0_gp

       do nc=1,meta%ncentres(ns) 

          din(:)=data%in(:,sum(meta%ncentres(1:ns-1))+nc)

          es(1:1)=es(1:1)+nn_forward(net_use,din)

       end do

       edata=data%out(1,sum(meta%ncentres(1:ns)))*meta%ncentres(ns) 

       dcost=(es(1)-edata)**2/2/real(meta%nstruct,gp)*max(1e-4_gp,abs(es(1)-edata))**(pval-2)

       if(present(emin)) dcost=dcost*max(wmin,exp(-(edata-emin)/kT))

       cost=cost+dcost

    end do
    !$omp end parallel do

    call cpu_time(end_time)

    cost_time=cost_time+end_time-start_time
    
  end function forge_cost

  function forge_gradient(net,data,meta,cost,batch) result(grad)

    type(network),                   intent(inout) :: net
    type(dataset),                   intent(in)    :: data
    type(metadata),                  intent(in)    :: meta
    real(kind=gp),                   intent(out)   :: cost
    integer, optional, dimension(:), intent(in)    :: batch

    real(kind=gp) :: grad(net%nwgts)
    real(kind=gp) :: din(net%inputs)

    real(kind=gp), allocatable, save, dimension(:)   :: g
    real(kind=gp), allocatable, save, dimension(:,:) :: back

    type(network), save :: net_use

    real(kind=gp) :: es(1)

    integer :: n,ns,nc

    integer, allocatable, dimension(:) :: batch_use

    if(present(batch)) then
       allocate(batch_use(size(batch)))
       batch_use=batch
    else
       allocate(batch_use(meta%nstruct))
       do n=1,meta%nstruct
          batch_use(n)=n
       end do
    end if

    if(.not.allocated(g)) &
         allocate(g(net%nwgts),back(net%nwgts,1))

    cost=0.0_gp
    g=0.0_gp

    !$omp parallel do private(es,din,nc,back,net_use) reduction(+:cost,g) schedule(dynamic)
    do ns=1,size(batch_use)

       net_use=net

       es(1:1)=0.0_gp

       back=0.0_gp

       do nc=1,meta%ncentres(ns) 

          din(:)=data%in(:,sum(meta%ncentres(1:batch_use(ns)-1))+nc)

          es(1:1)=es(1:1)+nn_forward(net_use,din)

          back=back+nn_backward(net_use)

       end do

       g=g-back(:,1)*(es(1)-data%out(1,sum(meta%ncentres(1:batch_use(ns))))*meta%ncentres(ns))

       cost=cost+(es(1)-data%out(1,sum(meta%ncentres(1:batch_use(ns))))*meta%ncentres(ns))**2/2

    end do
    !$omp end parallel do

    grad=g/real(size(batch_use),gp)

    cost=cost/real(size(batch_use),gp)

  end function forge_gradient

  function forge_lmstep(net,data,meta,lambda,cost,emin) result(step)

    type(network),  intent(inout) :: net
    type(dataset),  intent(in)    :: data
    type(metadata), intent(in)    :: meta
    real(kind=gp),  intent(in)    :: lambda
    real(kind=gp),  intent(out)   :: cost
    real(kind=gp), optional, intent(in) :: emin

    real(kind=gp) :: step(net%nwgts)
    real(kind=gp) :: din(net%inputs)

    type(network), save :: net_use

    real(kind=gp) :: es(1),edata

    real(kind=gp), allocatable, dimension(:), save :: bvec,deltay,W

    real(kind=gp), allocatable, dimension(:,:), save :: back,jac,jj

    real(kind=gp) :: start_time,end_time,dsyin_time,dsyout_time
    
    integer :: ns,nc,nw,npos

    integer, allocatable, dimension(:) :: chunksize

    call cpu_time(start_time)
    
    if(.not.allocated(deltay)) & ! ** These need to be allocated and saved for openmp,
                                !    apparently there is no automatic static allocation
         allocate(deltay(meta%nstruct),jac(net%nwgts,meta%nstruct),jj(net%nwgts,net%nwgts),&
         back(net%nwgts,1),bvec(net%nwgts),W(meta%nstruct))

    cost=0.0_gp

    !$omp parallel do private(es,din,nc,back,net_use,edata) reduction(+:cost) schedule(dynamic)
    do ns=1,meta%nstruct

       net_use=net

       es(1:1)=0.0_gp

       back=0.0_gp

       do nc=1,meta%ncentres(ns) 

          din(:)=data%in(:,sum(meta%ncentres(1:ns-1))+nc)
          
          es(1:1)=es(1:1)+nn_forward(net_use,din)
          
          back=back+nn_backward(net_use)
          
       end do

       edata=data%out(1,sum(meta%ncentres(1:ns)))*meta%ncentres(ns)

       if(present(emin)) then
          W(ns)=max(wmin,exp(-(edata-emin)/kT))*max(1e-4_gp,abs(edata-es(1)))**(pval-2)       
       else
          W(ns)=max(1e-4_gp,abs(edata-es(1)))**(pval-2)
       end if

       jac(:,ns)=back(:,1)*sqrt(W(ns)) 

       deltay(ns)=(edata-es(1))*sqrt(W(ns)) 

       cost=cost+(edata-es(1))**2/2*W(ns)

    end do
    !$omp end parallel do

    cost=cost/real(meta%nstruct,gp)

    allocate(chunksize(nchunks))

    chunksize=meta%nstruct/nchunks

    do nc=1,meta%nstruct-sum(chunksize)
       chunksize(nc)=chunksize(nc)+1
    end do

    call cpu_time(dsyin_time)
    
    jj=0.0_gp

    !$omp parallel do private(npos) reduction(+:jj) schedule(dynamic)
    do nc=1,nchunks

       npos=sum(chunksize(1:nc-1))+1

       if(gp.eq.sp) then
          call ssyrk('U','N',net%nwgts,chunksize(nc),1.0_gp,jac(1,npos),net%nwgts,1.0_gp,jj(1,1),net%nwgts)
       else
          call dsyrk('U','N',net%nwgts,chunksize(nc),1.0_gp,jac(1,npos),net%nwgts,1.0_gp,jj(1,1),net%nwgts)
       end if
       
    end do
    !$omp end parallel do

    call cpu_time(dsyout_time)
    
    dsyrk_time=dsyrk_time+dsyout_time-dsyin_time

    do nw=1,net%nwgts
       jj(nw,nw)=jj(nw,nw)+lambda*max(1.0_gp,jj(nw,nw))
    end do
    
    call cpu_time(end_time)

    jac_time=jac_time+end_time-start_time
    
    bvec=matmul(jac,deltay) !! ** This still needs parallelising

    call cpu_time(start_time)

    bvec_time=bvec_time+start_time-end_time
    
    step=-nn_solve(jj,bvec)

    call cpu_time(end_time)

    solve_time=solve_time+end_time-start_time
    
  end function forge_lmstep

  subroutine forge_test_gradient(net,dat,met)

    type(network),  intent(inout) :: net
    type(dataset),  intent(in)    :: dat
    type(metadata), intent(in)    :: met


    real(kind=gp), dimension(net%nwgts) :: w0,g,g0

    real(kind=gp) :: step=1e-4_gp,c,cp,cm

    real(kind=gp) :: start_time,end_time
    
    integer :: nw

    w0=net%weights

    do nw=1,net%nwgts

       net%weights=w0
       net%weights(nw)=net%weights(nw)+step

       cp=forge_cost(net,dat,met)

       net%weights=w0
       net%weights(nw)=net%weights(nw)-step

       cm=forge_cost(net,dat,met)

       g(nw)=-(cp-cm)/2.0_gp/step

    end do

    net%weights=w0

    write (stderr,'(10f10.5)') g
    g0=forge_gradient(net,dat,met,c)
    write (stderr,'(10f10.5)') g0
    write (stderr,*) norm2(g0-g)

  end subroutine forge_test_gradient

  subroutine forge_test(ctask,net,dat,met)

    character(len=12), intent(in)  :: ctask
    type(network),   intent(inout) :: net
    type(dataset),   intent(in)    :: dat
    type(metadata),  intent(in)    :: met

    real(kind=gp) :: e_predict(met%nstruct),es(1),mae,rmse,emax,abserr,tick,din(net%inputs),dmax

    real(kind=gp) :: start_time,end_time
    
    integer :: ns,nc,nv,stat

    character(len=40)  :: ctemp
    character(len=240) :: ctemp2
    character(len=240) :: cmax

    call cpu_time(start_time)
    
    if(plot) then

       write(ctemp,*) npot
       
       open(unit=unit_delta,file=trim(seedname)//'-'//trim(adjustl(ctemp))//'-'//trim(met%setname)&
            //'.delta',form='formatted',status='unknown')
       open(unit=unit_agr,file=trim(seedname)//'-'//trim(adjustl(ctemp))//'-'//trim(met%setname)&
            //'.agr',form='formatted',status='unknown')

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

          din(:)=dat%in(:,nv)
          es(1:1)=es(1:1)+nn_forward(net,din)*net%out_sd+net%out_mean(1)

       end do

       e_predict(ns)=es(1)

       if(plot) then
          write (unit_agr,*) met%enthalpy(ns)/met%natoms(ns),e_predict(ns)/met%ncentres(ns)
       end if

    end do

    if(plot) then

       write (unit_agr,*) '&'

       dmax=1.1_gp*maxval(abs(met%enthalpy/met%natoms-e_predict/met%ncentres))
       call delta_head(minval(met%enthalpy/met%natoms),maxval(met%enthalpy/met%natoms),&
            0.0_gp,dmax,tick,dmax/5.0_gp)

    end if

    emax=0.0_gp
    mae=0.0_gp
    rmse=0.0_gp
    do ns=1,met%nstruct
       abserr=abs(e_predict(ns)/met%ncentres(ns)-met%enthalpy(ns)/met%natoms(ns))
       rmse=rmse+abserr**2
       mae=mae+abserr
       if(plot) write (unit_delta,*) &
            met%enthalpy(ns)/met%natoms(ns),abs(e_predict(ns)/met%ncentres(ns)-met%enthalpy(ns)/met%natoms(ns))
       if(abserr.gt.emax) then
          emax=abserr
          cmax=met%label(ns)
       end if
    end do
    if(plot) write(unit_delta,*)

    rmse=sqrt(rmse/real(met%nstruct,gp))
    mae=mae/real(met%nstruct,gp)

    if(plot) then
       write (unit_delta,*) '&'
       write (unit_delta,*) minval(met%enthalpy/met%natoms),rmse
       write (unit_delta,*) maxval(met%enthalpy/met%natoms),rmse
       write (unit_delta,*) '&'
       write (unit_delta,*) minval(met%enthalpy/met%natoms),mae
       write (unit_delta,*) maxval(met%enthalpy/met%natoms),mae
    end if

    write (stdout,*)
    write (ctemp,'(a)') ctask//trim(met%setname)//' RMSE/MAE:  '
    write (stdout,'(a,2f9.2,a)',advance='no') ctemp,rmse*1000.0_gp,mae*1000.0_gp,' meV'
    write (stdout,'(a,f7.5)',advance='no')  ' Spearman : ',spearman(e_predict/met%ncentres,met%enthalpy/met%natoms)
    write (stdout,'(a,f10.2,a,a)') ' Max : ', emax*1000.0_gp,' meV for ',trim(cmax)

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
    write (stdout,'(a,2f9.2,a)',advance='no') ctemp,rmse,mae,'  eV'
    write (stdout,'(a,f7.5)',advance='no')  ' Spearman : ',spearman(e_predict/met%ncentres*met%natoms,met%enthalpy)
    write (stdout,'(a,f10.2,a,a,a)') ' Max : ', emax,'  eV for ',trim(cmax),' )'

    if(plot) then

       close(unit_delta)
       close(unit_agr)

       if(pdf) then

          write(ctemp,*) npot
          
          ctemp2="gracebat "//trim(seedname)//'-'//trim(adjustl(ctemp))//'-'//trim(met%setname)&
               //".agr -hdevice PDF -hardcopy -printfile "&
               //trim(seedname)//'-'//trim(adjustl(ctemp))//'-'//trim(met%setname)//".pdf"
          stat=0
          call system(ctemp2,stat)
          ctemp2='Problem executing external command :: '//trim(ctemp2)
          if (stat.ne.0) then
             write (stderr,'(a)') trim(ctemp2)
             stop
          end if
          if(.not.quiet) write (stderr,*) " saved file "//trim(seedname)//'-'//trim(adjustl(ctemp))&
               //'-'//trim(met%setname)//".pdf"

       end if

    end if

    call cpu_time(end_time)

    test_time=test_time+end_time-start_time

  end subroutine forge_test

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
    write(unit_agr,'(a)') '@    xaxis  label "Energy (eV)"'
    write(unit_agr,'(a)') '@    xaxis  label font 4'
    write(unit_agr,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_agr,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_agr,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_agr,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_agr,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_agr,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_agr,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_agr,'(a)') '@    yaxis  label "Predicted Energy (eV)"'
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
    write(unit_delta,'(a)') '@    xaxis  label "Energy (eV)"'
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
    write(unit_delta,'(a)') '@    s1 line linestyle 1'
    write(unit_delta,'(a)') '@    s1 line linewidth 1.0'
    write(unit_delta,'(a)') '@    s1 line color 2'
    write(unit_delta,'(a)') '@    s1 line pattern 1'

    write(unit_delta,'(a)') '@    s2 line type 1'
    write(unit_delta,'(a)') '@    s2 line linestyle 1'
    write(unit_delta,'(a)') '@    s2 line linewidth 1.0'
    write(unit_delta,'(a)') '@    s2 line color 3'
    write(unit_delta,'(a)') '@    s2 line pattern 1'


  end subroutine delta_head

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

    nnodes=0
    nlayers=0
    if(any(cargs.eq."-nn")) then               ! * Number of nodes in each layer
       do i=1,num_args-1
          if(cargs(i).eq."-nn") exit
       end do
       do nlayers=1,max_layers
          if(i+nlayers.gt.num_args) exit
          read(cargs(i+nlayers),*,err=102,end=102) nnodes(nlayers)
       end do
102    nlayers=nlayers-1
    end if

    trainfile='training'
    if(any(cargs.eq."-tr")) then               ! * File name for training data
       do i=1,num_args-1
          if(cargs(i).eq."-tr") exit
       end do
       read(cargs(i+1),'(a)') trainfile
    end if

    validfile='validation'
    if(any(cargs.eq."-va")) then               ! * File name for validation data
       do i=1,num_args-1
          if(cargs(i).eq."-va") exit
       end do
       read(cargs(i+1),'(a)') validfile
    end if

    testfile='testing'
    if(any(cargs.eq."-te")) then               ! * File name for testing data
       do i=1,num_args-1
          if(cargs(i).eq."-te") exit
       end do
       read(cargs(i+1),'(a)') testfile
    end if

    thresh=1e-10_gp
    if(any(cargs.eq."-th")) then               ! * Convergence hreshold
       do i=1,num_args-1
          if(cargs(i).eq."-th") exit
       end do
       read(cargs(i+1),*) thresh
    end if

    nchunks=nomp
    if(any(cargs.eq."-nc")) then               ! * Chunk size
       do i=1,num_args-1
          if(cargs(i).eq."-nc") exit
       end do
       read(cargs(i+1),*) nchunks
    end if

    nsteps=1000
    if(any(cargs.eq."-n")) then                ! * Number of steps
       do i=1,num_args-1
          if(cargs(i).eq."-n") exit
       end do
       read(cargs(i+1),*) nsteps
    end if

    nbatch=0
    if(any(cargs.eq."-b")) then                ! * Batch size
       do i=1,num_args-1
          if(cargs(i).eq."-b") exit
       end do
       read(cargs(i+1),*) nbatch
    end if

    eta=1e-0_gp
    if(any(cargs.eq."-eta")) then              ! * Learning rate
       do i=1,num_args-1
          if(cargs(i).eq."-eta") exit
       end do
       read(cargs(i+1),*) eta
    end if

    if(nlayers.gt.0) then
       act='tanh'
    else
       act='linear'
    end if
    if(any(cargs.eq."-a")) then                ! * Activation function
       do i=1,num_args-1
          if(cargs(i).eq."-a") exit
       end do
       read(cargs(i+1),*) act
    end if

    opt='lm'
    if(any(cargs.eq."-o")) then                ! * Optimisation algorithm
       do i=1,num_args-1
          if(cargs(i).eq."-o") exit
       end do
       read(cargs(i+1),*) opt
    end if

    if(opt.eq.'lm') then
       lambda_min=epsilon(1.0_gp)
    else if(opt.eq.'irls') then
       lambda_min=1e-4_gp
    else
       lambda_min=1e-13_gp
    end if
    if(any(cargs.eq."-lmin")) then                ! * Minimum lambda
       do i=1,num_args-1
          if(cargs(i).eq."-lmin") exit
       end do
       read(cargs(i+1),*) lambda_min
    end if

    pval=1.25_gp
    if(any(cargs.eq."-p")) then                ! * P value for IRLS
       do i=1,num_args-1
          if(cargs(i).eq."-p") exit
       end do
       read(cargs(i+1),*) pval
    end if

    wmin=1e-3_gp
    if(any(cargs.eq."-w")) then                ! * Minimum weight for IRLS
       do i=1,num_args-1
          if(cargs(i).eq."-w") exit
       end do
       read(cargs(i+1),*) wmin
    end if

    kT=huge(1.0_gp)
    if(any(cargs.eq."-kT")) then               ! * Energy for Boltzmann weighting
       do i=1,num_args-1
          if(cargs(i).eq."-kT") exit
       end do
       read(cargs(i+1),*) kT
    end if

    nestop=20
    if(any(cargs.eq."-es")) then               ! * Number of steps for early stopping
       do i=1,num_args-1
          if(cargs(i).eq."-es") exit
       end do
       read(cargs(i+1),*) nestop
       earlystop=.true.
    end if

    netfile='forge.ddp'
    if(any(cargs.eq."-l")) then                ! * Load
       do i=1,num_args-1
          if(cargs(i).eq."-l") exit
       end do
       if(i.lt.num_args) then
          if(cargs(i+1)(1:1).ne.'-') read(cargs(i+1),*) netfile
       end if

       load=.true.
    end if

    seedname='forge'
    if(any(cargs.eq."-s")) then                ! * Seedname
       do i=1,num_args-1
          if(cargs(i).eq."-s") exit
       end do
       if(i.lt.num_args) then
          if(cargs(i+1)(1:1).ne.'-') read(cargs(i+1),'(a)') seedname
       end if
    end if

    maxpot=1
    if(any(cargs.eq."-numpot")) then           ! * Number of potentials to generate
       do i=1,num_args-1
          if(cargs(i).eq."-numpot") exit
       end do
       read(cargs(i+1),*) maxpot
    end if

    if(any(cargs.eq."-ne")) earlystop=.false.

    if(any(cargs.eq."-np")) plot=.false.

    if(any(cargs.eq."-fs")) filestop=.true.

    if(any(cargs.eq."-t")) trk=.true.

    if(any(cargs.eq."-pdf")) pdf=.true.

    if(any(cargs.eq."-q")) quiet=.true.



    deallocate(cargs)

    return

100 write (stderr,*) 'Usage: forge [-ompnp] [-nn] [-tr] [-va] [-te] [-th] [-nc] [-n] [-b] [-eta] [-lmin]'//&
         ' [-es] [-ne] [-np] [-fs] [-t] [-a] [-o] [-p] [-w] [-kT] [-l] [-s] [-q] [-pdf] [-numpot] [-h]'
    write (stderr,*) '-ompnp   : Number of omp threads (not mpi)'
    write (stderr,*) '-nn I    : Number of nodes in each layer'
    write (stderr,*) '-tr C    : File name for training data'
    write (stderr,*) '-va C    : File name for validation data'
    write (stderr,*) '-te C    : File name for testing data'
    write (stderr,*) '-th F    : Convergence threshold'
    write (stderr,*) '-nc I    : Chunk size'
    write (stderr,*) '-n I     : Number of steps'
    write (stderr,*) '-b I     : Batch size'
    write (stderr,*) '-eta F   : Learning rate'
    write (stderr,*) '-lmin F  : Minimum lambda'
    write (stderr,*) '-es      : Early stopping, number of steps'
    write (stderr,*) '-ne      : No early stopping'
    write (stderr,*) '-np      : No plotting'
    write (stderr,*) '-fs      : Use stop file'
    write (stderr,*) '-t       : Track the training'
    write (stderr,*) '-a       : Activation function'
    write (stderr,*) '-o       : Optimisation algorithm'
    write (stderr,*) '-p       : p value for IRLS'
    write (stderr,*) '-w       : Minimum weight for IRLS'
    write (stderr,*) '-kT      : Energy for Boltzmann weighting (eV)'
    write (stderr,*) '-l [C]   : Load data derived potential'
    write (stderr,*) '-s [C]   : Seedname'
    write (stderr,*) '-q       : Quiet - minimal output'
    write (stderr,*) '-pdf     : Output plots to pdf (requires grace)'
    write (stderr,*) '-numpot  : Number of potentials to generate'
    write (stderr,*) '-h       : Print this help message'

    stop

  end subroutine get_arguments

  subroutine banner()

    if(quiet) return

    write (stderr,'(a)') "                                                 "
    write (stderr,'(a)') "       ▄████  ████▄ █▄▄▄▄   ▄▀  ▄███▄            "
    write (stderr,'(a)') "       █▀   ▀ █   █ █  ▄▀ ▄▀    █▀   ▀           "
    write (stderr,'(a)') "       █▀▀    █   █ █▀▀▌  █ ▀▄  ██▄▄             "
    write (stderr,'(a)') "       █      ▀████ █  █  █   █ █▄   ▄▀          "
    write (stderr,'(a)') "        █             █    ███  ▀███▀            "
    write (stderr,'(a)') "         ▀           ▀                           "
    write (stderr,'(a)') "                                                 "                        
    write (stderr,'(a)') " Author: Chris J. Pickard, Cambridge, 2020-22 (c)"
    write (stderr,*)
    if(nlayers.gt.0) write (stderr,'(a,*((i0,1x)))') ' hidden nodes   : ',nnodes(1:nlayers)
    if(nlayers.gt.0) write (stderr,'(a,*((i0,1x)))')  ' hidden layers  : ',nlayers
    write (stderr,'(a,*((a,1x)))') ' activation     : ',act
    write (stderr,'(a,*((a,1x)))') ' optimisation   : ',opt
    write (stderr,'(a,*((g0.4,1x)))') ' irls p value   : ',pval
    write (stderr,'(a,*((i0,1x)))') ' max steps      : ',nsteps
    if(earlystop) then
       write (stderr,'(a,*((a,1x)))') ' early stopping : yes'
    else
       write (stderr,'(a,*((a,1x)))') ' early stopping : no'
    end if
    write (stderr,'(a,*((i0,1x)))')  ' nchunks        : ',nchunks
    if((opt.eq.'sgdm').or.(opt.eq.'adamw')) write (stderr,'(a,*((i0,1x)))')  ' batch size     : ',nbatch
    if((opt.eq.'sgdm').or.(opt.eq.'adamw')) write (stderr,'(a,*((f0.3,1x)))')  ' learning rate  : ', eta
    write (stderr,'(a,*((i0,1x)))')  ' num threads    : ', omp_get_max_threads()
    write (stderr,'(a,*((i0,1x)))')  ' max threads    : ', omp_get_num_procs()
    write (stderr,'(a,*((i0,1x)))')  ' num potentials : ', maxpot   
    write (stderr,'(a,*((i0,1x)))')  ' random seed    : ', seed
    write (stderr,*)

  end subroutine banner

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

end program forge
