!==================================================================================!
!                                      ddp                                         !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This module reads, knows, and evaluates the data derived potential               !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module ddp_new

   use constants
   use cell
   use nn
   use spglib_f08

   implicit none

   private

   public :: read_ddp
   public :: eval_ddp
   public :: rcore
   public :: acore
   public :: rmin,rcut
   public :: ball_time,feat_time,netw_time

   integer,parameter :: unit_ddp=22
   integer,parameter :: unit_eddp=23

   ! ** Data defining the neural network potential

   real(kind=gp) :: rcut(5),rcut_sq(5),rcore,acore,rmin

   real(kind=gp),allocatable,dimension(:)      :: pow
   real(kind=gp),allocatable,dimension(:)      :: ddpweight

   character(len=120),allocatable,dimension(:) :: ddpfile
   character(len=6),allocatable,dimension(:)   :: spec

   integer       :: nddp,ncomp,ntwo,nthree,nfour,nfive,npow,nfeat

   integer,allocatable,dimension(:,:)   :: npair
   integer,allocatable,dimension(:,:,:) :: ntriple

   logical       :: twobody=.false.
   logical       :: threebody=.false.
   logical       :: full=.false.
   logical       :: ensemble=.false.

   type(network),allocatable,dimension(:) :: ddpot

   real(kind=gp) :: ball_time=0.0_gp
   real(kind=gp) :: feat_time=0.0_gp
   real(kind=gp) :: netw_time=0.0_gp

   !---------------------------------------------------!

   ! ** Private variables

   ! Ball data

   integer :: nbmax,nbmax_save=0
   integer,allocatable :: nball(:),iball(:,:)
   character(len=6),allocatable :: ballnames(:,:)
   real(kind=gp),allocatable :: lvec(:,:),ll(:),ball12(:,:,:),mod_ball12(:,:)

   ! Feature data

   real(kind=gp),allocatable :: features(:,:),forcefeat(:,:,:,:),stressfeat(:,:,:,:)
   real(kind=gp),allocatable :: fcore(:,:)
   real(kind=gp) :: score(3,3),ecore

contains

   subroutine read_ddp()

      integer  :: nc,n1,n2,n3,n,m,mm,nc1,nc2
      character(len=120) :: ctemp
      character(len=120) :: comp

      inquire(file=trim(seedname)//".eddp",exist=ensemble)

      if(ensemble) then

         open(unit=unit_eddp,file=trim(seedname)//".eddp",status="old",form="formatted",err=999)

         nddp=0

         do
            read(unit_eddp,*,end=90,err=999) ctemp
            nddp=nddp+1
         enddo

90       rewind(unit_eddp)

         allocate(ddpfile(nddp),ddpweight(nddp))

         do n=1,nddp
            read(unit_eddp,'(a)',err=999) ctemp
            ctemp=adjustl(ctemp)
            ddpfile(n)=ctemp(1:index(ctemp," "))
            read(ctemp(index(ctemp," ")+1:),*) ddpweight(n)
         enddo

         ddpweight=ddpweight/sum(ddpweight)

         close(unit_eddp)

      else

         nddp=1

         allocate(ddpfile(nddp),ddpweight(nddp))

         ddpfile(1)=trim(seedname)//".ddp"
         ddpweight(1)=1.0_gp

      endif

      allocate(ddpot(nddp))

      do n=1,nddp
         call nn_load(ddpot(n),ddpfile(n))
      enddo

      if(.not.all(ddpot(1)%cmt.eq.ddpot(:)%cmt)) stop 'ddp: inconsistent ensemble potential'

      if(.not.quiet) then
         write(stderr,*)
         if(ensemble) then
            write(stderr,'(a,i0,a)') ' ensemble data derived potential summary (',nddp,'):'
         else
            write(stderr,*) 'data derived potential summary: '
         endif
         write(stderr,*)
         call nn_summary(ddpot(1))
         write(stderr,*)
      endif

      read(ddpot(1)%cmt,*) ctemp,comp,ctemp,rcut,ctemp,ntwo,nthree,nfour,nfive

      twobody=ntwo.gt.0
      threebody=nthree.gt.0

      if(nfour.gt.0) stop 'ddp: four body not implemented'
      if(nfive.gt.0) stop 'ddp: five body not implemented'

      if(threebody) then
         if(nthree.eq.ntwo**2) then
            full=.true.
         else
            stop 'diagonal not implemented'
         endif
      endif

      rcut_sq=rcut**2

      if(rcore.lt.0.0_gp) rcore=rcut(1)

      ncomp=1
      do nc=1,len_trim(comp)
         if(comp(nc:nc).eq.'-') then
            ncomp=ncomp+1
            comp(nc:nc)=' '
         endif
      enddo

      allocate(spec(ncomp),npair(ncomp,ncomp),ntriple(ncomp,ncomp,ncomp))

      read(comp,*) spec

      nfeat=ncomp+ntwo*ncomp**2+nthree*ncomp**2*(ncomp+1)/2

      npow=ntwo

      allocate(pow(npow))

      read(ddpot(1)%cmt,*) ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,ctemp,pow

      ! ** Set up indexing

      n=0
      do n1=1,ncomp
         do n2=1,ncomp
            n=n+1
            npair(n2,n1)=n
         enddo
      enddo

      do n1=1,ncomp
         do n2=1,ncomp
            do n3=1,ncomp

               nc1=min(n2,n3)
               nc2=max(n2,n3)

               mm=0
               do n=1,ncomp
                  do m=n,ncomp
                     mm=mm+1
                     if((n.eq.nc1).and.(m.eq.nc2)) exit
                  enddo
                  if((n.eq.nc1).and.(m.eq.nc2)) exit
               enddo

               ntriple(n3,n2,n1)=(n1-1)*ncomp*(ncomp+1)/2+mm

            enddo
         enddo
      enddo

      return

999   stop 'read_ddp: problem opening/reading eddp file'

   endsubroutine read_ddp

   subroutine eval_ddp(str,e,f,s,variance)

      type(structure),intent(in)  :: str
      real(kind=gp),intent(out) :: e
      real(kind=gp),optional,intent(out) :: f(3,str%num_ions)
      real(kind=gp),optional,intent(out) :: s(3,3)
      real(kind=gp),optional,intent(out) :: variance

      integer       :: ion_i,ion_j,i,j,n,np

      real(kind=gp) :: sloc(3,3),de,eddp(nddp)

      real(kind=gp),save,allocatable :: dndf(:),dndfw(:)

      logical :: fcalc,scalc,vcalc

      real(kind=gp) :: start_time,end_time

      type(network),save,allocatable,dimension(:,:) :: ddpot_use

      fcalc=present(f)

      scalc=present(s)

      if(.not.fcalc.and.scalc) stop 'eval_ddp: forces required for stresses'

      vcalc=present(variance)

      if(allocated(ddpot_use)) then
         if(size(ddpot_use,2).ne.str%num_ions) deallocate(ddpot_use)
      endif

      if(.not.allocated(ddpot_use)) then
         allocate(ddpot_use(nddp,str%num_ions))
         !$omp parallel do schedule(dynamic)
         do n=1,str%num_ions
            ddpot_use(:,n)=ddpot(:)
         enddo
         !$omp end parallel do
      endif

      ! * --------------------------------

      ! ** Compute the balls

      call ball_ddp(str)

      ! ** Compute the features

      call feature_ddp(str,fcalc,scalc)

      ! ** Predict the forces, stress and energy

      call cpu_time(start_time)

      if(.not.allocated(dndf).and.fcalc) allocate(dndf(nfeat),dndfw(nfeat))

      eddp=0.0_gp
      e=ecore
      if(fcalc) f=fcore
      if(scalc) sloc=score

      !$omp parallel do private(dndf,dndfw,de) reduction(+:e,f,sloc,eddp) schedule(dynamic)
      do ion_i=1,str%num_ions

         if(fcalc) dndfw=0.0_gp

         do np=1,nddp

            if(fcalc) then
               de=predict_ddp(ddpot_use(np,ion_i),features(:,ion_i),dndf)
            else
               de=predict_ddp(ddpot_use(np,ion_i),features(:,ion_i))
            endif

            e=e+de*ddpweight(np)

            if(vcalc) eddp(np)=eddp(np)+de

            if(fcalc) dndfw=dndfw+dndf*ddpweight(np)

         enddo

         if(fcalc) then
            do ion_j=1,nball(ion_i)
               f(:,iball(ion_j,ion_i))=f(:,iball(ion_j,ion_i))+matmul(dndfw(:),forcefeat(:,:,ion_j,ion_i))
            enddo
         endif

         if(scalc) then
            do i=1,3
               do j=1,3
                  sloc(j,i)=sloc(j,i)+dot_product(dndfw(:),stressfeat(:,j,i,ion_i))
               enddo
            enddo
         endif

      enddo
      !$omp end parallel do

      if(scalc) s=sloc

      if(vcalc) variance=sum((eddp-sum(eddp)/size(eddp))**2)/size(eddp)

      call cpu_time(end_time)

      netw_time=netw_time+end_time-start_time

      ! ----------------------------------------------------------------------

      call symmetrise_ddp(str,f,s)

   endsubroutine eval_ddp

   subroutine ball_ddp(str)

      type(structure),intent(in)  :: str

      real(kind=gp) :: ion_positions_use(3,str%num_ions)
      real(kind=dp) :: car_use_dp(3,3)
      real(kind=gp) :: car_use(3,3),rec_use(3,3)

      real(kind=gp) :: h(3),r12(3),rij(3),mod_r12_sq,mod_rij

      real(kind=gp) :: start_time,end_time

      integer :: ion_i,ion_j,na,nb,nc,nn,nn0,nna,nnb,nnc,nnmx
      integer :: info

      call cpu_time(start_time)

      ! * Niggli reduce lattice 
      
      car_use_dp=transpose(str%lattice_car)
      
      info=spg_niggli_reduce(car_use_dp,1e-4_dp)
      
      car_use=transpose(car_use_dp)
      
      rec_use=cell_rec(car_use)

      ! * Map to cell for ball generate

      ion_positions_use=str%ion_positions
      
      do ion_i=1,str%num_ions
         ion_positions_use(:,ion_i)=matmul(rec_use,ion_positions_use(:,ion_i))
         ion_positions_use(:,ion_i)=ion_positions_use(:,ion_i)-real(floor(ion_positions_use(:,ion_i),8),gp)
         ion_positions_use(:,ion_i)=matmul(car_use,ion_positions_use(:,ion_i))
      enddo

      ! Determine the supercell required

      h=car2habc(car_use)

      nna=int(2*rcut(1)/h(1))+1
      nnb=int(2*rcut(1)/h(2))+1
      nnc=int(2*rcut(1)/h(3))+1

      if(cluster) then
         nna=0; nnb=0; nnc=0
      endif

      nnmx=(2*nna+1)*(2*nnb+1)*(2*nnc+1)

      if(any((/nna,nnb,nnc/).gt.10)) then
         write(stderr,'(3f10.5)') car_use
         write(stderr,'(a,*(1x,i0))') 'nna,nnb,nnc=',nna,nnb,nnc
         stop 'eval_ddp: nna,nnb,nnc too big'
      endif

      if(allocated(ll)) then
         if((nnmx.ne.size(ll)).or.(str%num_ions.ne.size(nball))) deallocate(lvec,ll,ball12,mod_ball12,nball,iball,ballnames)
      endif

      if(.not.(allocated(lvec))) allocate(lvec(3,nnmx),ll(nnmx))

      nn=0
      nn0=0
      do na=-nna,nna
         do nb=-nnb,nnb
            do nc=-nnc,nnc
               nn=nn+1
               lvec(1:3,nn)=na*car_use(1:3,1)+nb*car_use(1:3,2)+nc*car_use(1:3,3)
               ll(nn)=norm2(lvec(1:3,nn))
               if(((na==0).and.(nb==0).and.(nc==0))) nn0=nn
            enddo
         enddo
      enddo

      ! ** Compute ball

      if(.not.(allocated(ball12))) then !! ** N^2 memory usage here
         allocate(ball12(3,str%num_ions*nnmx,str%num_ions),mod_ball12(str%num_ions*nnmx,str%num_ions),nball(str%num_ions), &
                  iball(str%num_ions*nnmx,str%num_ions),ballnames(str%num_ions*nnmx,str%num_ions))
      endif

      nball=0
      !$omp parallel do private(rij,mod_rij,r12,mod_r12_sq) schedule(dynamic)
      do ion_i=1,str%num_ions

         do ion_j=1,str%num_ions

            rij(:)=ion_positions_use(:,ion_i)-ion_positions_use(:,ion_j)

            mod_rij=norm2(rij)

            do nn=1,nnmx

               r12(:)=rij(:)-lvec(:,nn)

               mod_r12_sq=r12(1)**2+r12(2)**2+r12(3)**2
               
               ! Check separation is not too big

               if(mod_r12_sq.lt.rcut_sq(1)) then

                  nball(ion_i)=nball(ion_i)+1

                  iball(nball(ion_i),ion_i)=ion_j

                  ball12(:,nball(ion_i),ion_i)=r12(:)

                  mod_ball12(nball(ion_i),ion_i)=sqrt(mod_r12_sq)

                  if(.not.((nn.eq.nn0).and.(ion_i.eq.ion_j))) then
                     if(mod_ball12(nball(ion_i),ion_i).lt.rmin) stop 'eval_ddp: close contact detected'
                  endif

                  ballnames(nball(ion_i),ion_i)=str%ion_names(ion_j)

               endif

            enddo !nn

         enddo  !ion_j

      enddo !ion_i
      !$omp end parallel do

      nbmax=maxval(nball)

      call cpu_time(end_time)

      ball_time=ball_time+end_time-start_time

   endsubroutine ball_ddp

   subroutine feature_ddp(str,fcalc,scalc)

      type(structure),intent(in)  :: str
      logical,intent(in) ::fcalc
      logical,intent(in) ::scalc

      integer       :: ion_i,ion_j,ion_k,ion_ib
      integer       :: nci(1),ncj(1),nck(1),nvec(1)
      integer       :: i,j,npos,n1,n2,n

      real(kind=gp) :: r12(3),r13(3),r23(3),rr12(3,3),rr13(3,3),rr23(3,3)

      real(kind=gp) :: mod_r12,mod_r13,mod_r23,mod_r23_sq

      real(kind=gp),save,allocatable :: fact1(:),fact2(:),fc12p(:),fc13p(:),fc23p(:),f12(:),f13(:),f23(:)
      real(kind=gp),save,allocatable :: fact12(:,:),fact13(:,:),fact23(:,:),sfact12(:,:,:),sfact13(:,:,:),sfact23(:,:,:)
      real(kind=gp),save,allocatable :: sfact(:,:,:),ffact1(:,:),ffact2(:,:),ffact3(:,:)
      
      real(kind=gp) :: start_time,end_time

      if(.not.allocated(ball12)) stop 'feature_ddp: called before balls generated'

      call cpu_time(start_time)

      ! ** Allocate feature vectors

      if((nbmax.ne.nbmax_save).and.allocated(forcefeat)) deallocate(forcefeat)

      if(allocated(features)) then
         if(size(features,2).ne.str%num_ions) deallocate(features)
      endif

      if(allocated(fcore)) then
         if(size(fcore,2).ne.str%num_ions) deallocate(forcefeat,fcore)
      endif

      if(allocated(stressfeat)) then
         if(size(stressfeat,4).ne.str%num_ions) deallocate(stressfeat)
      endif

      ! ** Compute features

      if(.not.allocated(features)) allocate(features(nfeat,str%num_ions))

      if(.not.allocated(forcefeat).and.fcalc) then
         allocate(forcefeat(nfeat,3,nbmax,str%num_ions))
         nbmax_save=nbmax
      endif

      if(.not.allocated(fcore).and.fcalc) allocate(fcore(3,str%num_ions))
      if(.not.allocated(fact12).and.fcalc) allocate(fact12(npow,3),fact13(npow,3),fact23(npow,3))
      if(.not.allocated(fact1).and.fcalc) allocate(fact1(npow),fact2(npow*npow))
      if(.not.allocated(fc12p).and.fcalc) allocate(fc12p(npow),fc13p(npow),fc23p(npow),f12(npow),f13(npow),f23(npow))
      if(.not.allocated(ffact1).and.fcalc) allocate(ffact1(npow*npow,3),ffact2(npow*npow,3),ffact3(npow*npow,3))

      if(.not.allocated(stressfeat).and.scalc) allocate(stressfeat(nfeat,3,3,str%num_ions))
      if(.not.allocated(sfact12).and.scalc) allocate(sfact12(npow,3,3),sfact13(npow,3,3),sfact23(npow,3,3),sfact(npow*npow,3,3))


      ecore=0.0_gp
      if(scalc) score=0.0_gp

      !$omp parallel do private(nci,ncj,nck,npos,r12,r13,r23,mod_r12,mod_r13,mod_r23,mod_r23_sq,n,nvec) &
      !$omp private(fact12,fact13,fact23,sfact12,sfact13,sfact23,sfact,rr12,rr13,rr23,f12,f13,f23,ion_ib) &
      !$omp private(fc12p,fc13p,fc23p,fact1,fact2,ffact1,ffact2,ffact3) reduction(+:ecore,score) schedule(dynamic)
      do ion_i=1,str%num_ions

         features(:,ion_i)=0.0_gp
         if(fcalc) fcore(:,ion_i)=0.0_gp
         if(fcalc) forcefeat(:,:,:,ion_i)=0.0_gp
         if(scalc) stressfeat(:,:,:,ion_i)=0.0_gp

         nvec=findloc(iball(1:nball(ion_i),ion_i),ion_i)

         ion_ib=nvec(1)

         nci=findloc(spec,str%ion_names(ion_i))

         features(nci(1),ion_i)=1.0_gp

         if(twobody) then

            do ion_j=1,nball(ion_i)

               mod_r12=mod_ball12(ion_j,ion_i)

               if(mod_r12.le.tiny(1.0_gp)) cycle

               r12(:)=ball12(:,ion_j,ion_i)

               rr12=outer3(r12,r12)

               ! ** Add hard core repulsion

               if(mod_r12.lt.rcore) then

                  ecore=ecore+acore*fncore(mod_r12,rcore)

                  if(fcalc) fcore(:,ion_i)=fcore(:,ion_i)-2*acore*fncorep(mod_r12,rcore)*r12(:)/mod_r12

                  if(scalc) then
                     do i=1,3 !! REMOVE LOOPS
                        do j=1,3
                           score(j,i)=score(j,i)-acore*fncorep(mod_r12,rcore)*rr12(j,i)/mod_r12
                        enddo
                     enddo
                  endif

               endif

               ncj=findloc(spec,ballnames(ion_j,ion_i))

               npos=ncomp+(npair(ncj(1),nci(1))-1)*ntwo

               fc12p(:)=fcpow(mod_r12,rcut(2),pow)

               f12(:)=fcpfx(mod_r12,rcut(2))*pow(:)

               do i=1,3

                  if(fcalc) fact12(:,i)=f12(:)*r12(i)

                  if(scalc) then
                     do j=1,3
                        sfact12(:,j,i)=fact12(:,i)*r12(j)
                     enddo
                  endif

               enddo

               features(npos+1:npos+ntwo,ion_i)=features(npos+1:npos+ntwo,ion_i)+fc12p(:)

               ! ** Compute forces and stresses

               do i=1,3

                  if(fcalc) then
                     if(ion_ib.gt.0) forcefeat(npos+1:npos+ntwo,i,ion_ib,ion_i)=forcefeat(npos+1:npos+ntwo,i,ion_ib,ion_i)- &
                                                                                 fc12p(:)*fact12(:,i)

                     forcefeat(npos+1:npos+ntwo,i,ion_j,ion_i)=forcefeat(npos+1:npos+ntwo,i,ion_j,ion_i)+ &
                                                                fc12p(:)*fact12(:,i)
                  endif

                  if(scalc) then
                     do j=1,3
                        stressfeat(npos+1:npos+ntwo,j,i,ion_i)=stressfeat(npos+1:npos+ntwo,j,i,ion_i)-fc12p(:)*sfact12(:,j,i)
                     enddo
                  endif

               enddo

               if(threebody) then

                  do ion_k=ion_j+1,nball(ion_i)

                     if(ion_j.eq.ion_k) cycle

                     mod_r13=mod_ball12(ion_k,ion_i)

                     if(mod_r13.le.tiny(1.0_gp)) cycle

                     r13=ball12(:,ion_k,ion_i)

                     r23=r13-r12

                     mod_r23_sq=r23(1)**2+r23(2)**2+r23(3)**2

                     if(mod_r23_sq.gt.rcut_sq(3)) cycle

                     mod_r23=sqrt(mod_r23_sq)

                     rr13=outer3(r13,r13)

                     rr23=outer3(r23,r23)

                     nck=findloc(spec,ballnames(ion_k,ion_i))

                     npos=ncomp+ntwo*ncomp**2+(ntriple(nck(1),ncj(1),nci(1))-1)*nthree

                     f13(:)=fcpfx(mod_r13,rcut(3))*pow(:)
                     f23(:)=fcpfx(mod_r23,rcut(3))*pow(:)

                     fc13p(:)=fcpow(mod_r13,rcut(3),pow)
                     fc23p(:)=fcpow(mod_r23,rcut(3),pow)

                     do i=1,3

                        if(fcalc) then
                           fact13(:,i)=f13(:)*r13(i)
                           fact23(:,i)=f23(:)*r23(i)
                        endif

                        if(scalc) then
                           do j=1,3
                              sfact13(:,j,i)=fact13(:,i)*r13(j)
                              sfact23(:,j,i)=fact23(:,i)*r23(j)
                           enddo
                        endif

                     enddo

                     if(fcalc) then

                        do i=1,3

                           n=0
                           do n1=1,npow
                              do n2=1,npow
                                 n=n+1
                                 ffact1(n,i)=-fact12(n1,i)-fact13(n1,i)
                                 ffact2(n,i)=fact12(n1,i)-fact23(n2,i)
                                 ffact3(n,i)=fact13(n1,i)+fact23(n2,i)
                              end do
                           end do

                        end do

                     end if

                     if(scalc) then
                        do i=1,3
                           do j=1,3

                              n=0
                              do n1=1,npow
                                 do n2=1,npow
                                    n=n+1
                                    sfact(n,j,i)=sfact12(n1,j,i)+sfact13(n1,j,i)+sfact23(n2,j,i)
                                 end do
                              end do

                           end do
                        end do
                     end if

                     fact1(:)=fc12p(:)*fc13p(:)/npow

                     fact2=0.0_gp

                     if(gp.eq.sp) then
                        call sger(npow,npow,1.0_gp,fc23p(1),1,fact1(1),1,fact2(1),npow)
                     else
                        call dger(npow,npow,1.0_gp,fc23p(1),1,fact1(1),1,fact2(1),npow)
                     endif

                     ! ** Features

                     features(npos+1:npos+npow*npow,ion_i)=features(npos+1:npos+npow*npow,ion_i)+fact2(:)

                     ! ** Forces

                     if(fcalc) then
                        do i=1,3
                           forcefeat(npos+1:npos+npow*npow,i,ion_ib,ion_i)=forcefeat(npos+1:npos+npow*npow,i,ion_ib,ion_i)+&
                              fact2(:)*ffact1(:,i)
                            forcefeat(npos+1:npos+npow*npow,i,ion_j,ion_i)=forcefeat(npos+1:npos+npow*npow,i,ion_j,ion_i)+&
                              fact2(:)*ffact2(:,i)
                            forcefeat(npos+1:npos+npow*npow,i,ion_k,ion_i)=forcefeat(npos+1:npos+npow*npow,i,ion_k,ion_i)+&
                              fact2(:)*ffact3(:,i)
                        end do
                     end if

                     ! ** Stress

                     if(scalc) then
                        do i=1,3
                           do j=1,3
                              stressfeat(npos+1:npos+npow*npow,j,i,ion_i)=stressfeat(npos+1:npos+npow*npow,j,i,ion_i)-&
                                 fact2(:)*sfact(:,j,i)
                           enddo
                        enddo
                     endif

                  enddo !ion_k

               endif

            enddo  !ion_j

         endif

      enddo !ion_i
      !$omp end parallel do

      call cpu_time(end_time)

      feat_time=feat_time+end_time-start_time

   endsubroutine feature_ddp

   function predict_ddp(pot,x,gpdx)

      type(network),intent(inout) :: pot
      real(kind=gp),dimension(:),intent(in)    :: x
      real(kind=gp),optional,dimension(:),intent(out)   :: gpdx

      real(kind=gp) :: predict_ddp

      real(kind=gp) :: pred(1),back(pot%nwgts,1)

      ! ** Predict output

      pred=nn_forward(pot,x,scale=.true.)

      predict_ddp=pred(1)

      if(.not.present(gpdx)) return

      ! ** Compute gradient with respect to inputs

      back=-nn_backward(pot)

      if(gp.eq.sp) then
         call sgemv('T',pot%layers(2),nfeat,1.0_gp,pot%weights(1),pot%layers(2),back(nfeat*pot%layers(2)+1,1),1,0.0_gp,gpdx,1)
      else
         call dgemv('T',pot%layers(2),nfeat,1.0_gp,pot%weights(1),pot%layers(2),back(nfeat*pot%layers(2)+1,1),1,0.0_gp,gpdx,1)
      endif

      ! Scale and shift the results

      gpdx(:)=gpdx(:)*pot%out_sd(1)/pot%in_sd(:)

   endfunction predict_ddp

   subroutine symmetrise_ddp(str,f,s)

      type(structure),intent(in)  :: str
      real(kind=gp),optional,intent(inout) :: f(3,str%num_ions)
      real(kind=gp),optional,intent(inout) :: s(3,3)

      real(kind=gp) :: work(3,str%num_ions),swork(3,3),fav(3)

      integer :: ni,ns

      logical :: fcalc,scalc

      fcalc=present(f)

      scalc=present(s)

      if(fcalc) then

         ! * Balance forces

         fav(1)=sum(f(1,:))/real(str%num_ions,gp)
         fav(2)=sum(f(2,:))/real(str%num_ions,gp)
         fav(3)=sum(f(3,:))/real(str%num_ions,gp)

         f(1,:)=f(1,:)-fav(1)
         f(2,:)=f(2,:)-fav(2)
         f(3,:)=f(3,:)-fav(3)

      endif

      if(str%num_symm.gt.1) then

         if(fcalc) then

            ! * Symmetrise forces

            work=f
            f=0.0_gp

            do ni=1,str%num_ions
               do ns=1,str%num_symm
                  f(:,ni)=f(:,ni)+matmul(str%symm_ops(1:3,1:3,ns),work(:,str%ion_equiv(ns,ni)))
               enddo
               f(:,ni)=f(:,ni)/real(str%num_symm,gp)
            enddo

         endif

         if(scalc) then

            ! * Symmetrise stresses

            swork=s/real(str%num_symm,gp)
            s=0.0_gp

            do ns=1,str%num_symm
               s=s+matmul(str%symm_ops(1:3,1:3,ns),matmul(swork,transpose(str%symm_ops(1:3,1:3,ns))))
            enddo

         endif

      endif

      ! ** Constrain ions

      do ni=1,str%num_ions
         if(str%ion_cons(ni)) f(:,ni)=0.0_gp
      enddo

   endsubroutine symmetrise_ddp

   ! -----------------------------------------------------------

!    function fast_distance_2(b,a) result(d2)
!
!       real(kind=gp), intent(in) :: b(:,:)
!       real(kind=gp), intent(in) :: a(:,:)
!
!       real(kind=gp) :: d2(size(b,2),size(a,2))
!
!       real(kind=gp) :: l2a(size(a,2)),l2b(size(b,2))
!
!       integer :: na,nb
!
!       do na=1,size(a,2)
!          l2a(na)=dot_product(a(:,na),a(:,na))
!       end do
!
!       do nb=1,size(b,2)
!          l2b(nb)=dot_product(b(:,nb),b(:,nb))
!       end do
!
!        do na=1,size(a,2)
!          do nb=1,size(b,2)
!             d2(nb,na)=l2b(nb)+l2a(na)
!          end do
!       end do
!
!       if(gp.eq.sp) then
!          call sgemm('T','N',size(b,2),size(a,2),3,-2.0_gp,b(1,1),3,a(1,1),3,1.0_gp,d2(1,1),size(b,2))
!       else
!          call dgemm('T','N',size(b,2),size(a,2),3,-2.0_gp,b(1,1),3,a(1,1),3,1.0_gp,d2(1,1),size(b,2))
!       end if
!
!       d2=abs(d2)
!
!    end function

   elemental function fc(x,c)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp) :: fc

      ! Cutoff function

      if((x.gt.tiny(1.0_gp)).and.(x.lt.c)) then
         fc=2*(1-x/c)
      else
         fc=0.0_gp
      endif

   endfunction fc

   function fcpow(x,c,p)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp),intent(in) :: p(:)

      real(kind=gp) :: fcpow(size(p))

      if((x.gt.tiny(1.0_gp)).and.(x.lt.c)) then
         fcpow=exp(p*log(2*(1-x/c)))
      else
         fcpow=0.0_gp
      endif

   endfunction fcpow

   elemental function fcp(x,c)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp) :: fcp

      ! Cutoff function

      if((x.gt.tiny(1.0_gp)).and.(x.lt.c)) then
         fcp=-2/c
      else
         fcp=0.0_gp
      endif

   endfunction fcp

   elemental function fcpfx(x,c)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp) :: fcpfx

      ! Cutoff function

      if((x.gt.tiny(1.0_gp)).and.(x.lt.c)) then
         fcpfx=-1/x/(c-x)
      else
         fcpfx=0.0_gp
      endif

   endfunction fcpfx

   elemental function fcpp(x,c)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp) :: fcpp

      ! Cutoff function

      if((x.gt.tiny(1.0_gp)).and.(x.lt.c)) then
         fcpp=(-2/c)/(2*(1-x/c))/x
      else
         fcpp=0.0_gp
      endif

   endfunction fcpp

   function fncore(x,c)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp) :: fncore

      ! Hard core function

      if(x.lt.c) then
!!$       fncore=((c/x)**6-1)**2
!!$       fncore=((c/x)**3-1)**4
         fncore=(c/x-1)**12
      else
         fncore=0.0_gp
      endif

   endfunction fncore

   function fncorep(x,c)
      real(kind=gp),intent(in) :: x
      real(kind=gp),intent(in) :: c
      real(kind=gp) :: fncorep

      ! Hard core function

      if(x.lt.c) then
!!$       fncorep=-12*((c/x)**6-1)*(c/x)**6/x
!!$       fncorep=-12*((c/x)**3-1)**3*(c/x)**3/x
         fncorep=-12*(c/x-1)**11*c/x**2
      else
         fncorep=0.0_gp
      endif

   endfunction fncorep

   function detab(string)

      implicit none

      character(len=*),intent(in) :: string
      character(len=400)           :: detab

      integer                      :: i

      detab=string
      do while(scan(detab,char(9)).gt.0)
         i=scan(detab,char(9))
         detab=trim(detab(1:i-1))//'      '//trim(detab(i+1:))
      enddo

   endfunction detab

   function outer3(a,b) result(c)

      real(kind=gp),intent(in) :: a(3),b(3)

      real(kind=gp) :: c(3,3)

      c(1,1)=a(1)*b(1)
      c(2,1)=a(1)*b(2)
      c(3,1)=a(1)*b(3)

      c(1,2)=a(2)*b(1)
      c(2,2)=a(2)*b(2)
      c(3,2)=a(2)*b(3)

      c(1,3)=a(3)*b(1)
      c(2,3)=a(3)*b(2)
      c(3,3)=a(3)*b(3)

   endfunction outer3

endmodule ddp_new
