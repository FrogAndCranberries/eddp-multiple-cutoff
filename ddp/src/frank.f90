!==================================================================================!
!                                    frank                                         !
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
! This program generates environment vectors                                       !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
program frank

  use constants
  use omp_lib

  implicit none

  integer, parameter :: max_words=60
  integer, parameter :: max_species=200
  integer, parameter :: max_lines=1000 ! ** STDIN can only be read once, so need to estimate a maximum input length

  character(len=max_words*2), dimension(max_lines)           :: buff
  character(len=max_words*2), dimension(max_words,max_lines) :: buff_words
  integer,                    dimension(max_lines)           :: num_words
  
  integer :: num_lines,num_ions,num_species=0,num_super,num_body,num_comp,nomp

  integer, allocatable, dimension(:) :: ion_index,ion_superindex,num_feat,num_comb

  integer, allocatable, dimension(:,:)       :: indx_2b
  integer, allocatable, dimension(:,:,:)     :: indx_3b
  integer, allocatable, dimension(:,:,:,:)   :: indx_4b
  integer, allocatable, dimension(:,:,:,:,:) :: indx_5b
  
  real(kind=gp)                                :: lattice_abc(6),lattice_car(3,3),lattice_rec(3,3)
  real(kind=gp)                                :: lattice_volume
  real(kind=gp), allocatable, dimension(:)     :: ion_occ,ion_spin
  real(kind=gp), allocatable, dimension(:,:)   :: ion_fractional,ion_absolute,ion_supercluster
  real(kind=gp), allocatable, dimension(:,:,:) :: features

  character(len=3)                            :: species_names(max_species),comp_names(max_species)
  character(len=3), allocatable, dimension(:) :: ion_names,ion_supernames

  ! * =======================================================================

  integer :: status,npow=10,deltap
  
  real(kind=gp) :: rmax=3.75_gp,pmin=2.0_gp,pmax=15.0_gp

  real(kind=gp), allocatable, dimension(:) :: pow

  logical :: cluster=.false.,lj=.false.,mean=.false.

  ! * =======================================================================

  call get_arguments()

  do 

     lattice_car=0.0_gp
     lattice_rec=0.0_gp
     lattice_abc=0.0_gp
     lattice_volume=0.0_gp
     species_names=''
     
     if(allocated(ion_index)) then
        deallocate(ion_index,ion_superindex,ion_occ,ion_spin,ion_fractional,ion_absolute)
        deallocate(ion_supercluster,ion_names,ion_supernames,num_feat,num_comb)
        deallocate(features)
        num_species=0
     end if
     
     call read_buff(status)

     if(status.ne.0) exit
     
     call read_res()
     
     call consolidate()
     
     call niggli()
     
     call generate_supercluster()
     
     call generate_features()

  end do

contains

  subroutine generate_features()

    integer :: ni,nj,nk,nl,nm,n,npows(4),ncentres

    integer :: n1,n2,n3,n4,m,mm,o,p,nv(num_body)

    integer :: boff

    integer :: nspec(num_super)

    real(kind=gp) :: length2(num_super),rmax2
    real(kind=gp) :: elj,alpha,beta,enth,vol

    real(kind=gp), allocatable, dimension(:,:) :: dij,dij2
    
    character(len=600) :: ctemp

    character(len=60) :: centhalpy, cvolume

    allocate(dij(num_super,num_super),dij2(num_super,num_super))
    
    if(.not.(allocated(pow))) then
       allocate(pow(npow))

       if(npow.gt.1) then
          beta=(pmax/pmin)**(1.0_gp/real(npow-1,gp))
       else
          beta=1.0_gp
       endif
       alpha=pmin

       do n=1,npow
          if(.true.) then
             pow(n)=alpha*beta**(n-1)
          else
             pow(n)=pmin+real(n-1,gp)/real(npow-1,gp)*(pmax-pmin)
          end if
          write (ctemp,'(f8.3)') pow(n)
          read (ctemp,*) pow(n)
       end do

    end if

    rmax2=rmax**2
    
    ! ** Fast calculation of distance matrix

    do ni=1,num_super
       length2(ni)=dot_product(ion_supercluster(:,ni),ion_supercluster(:,ni))
    end do

    !$omp parallel do
    do ni=1,num_super
       do nj=1,num_super
          dij2(nj,ni)=length2(ni)+length2(nj)
       end do
    end do
    !$omp end parallel do

    if(gp.eq.sp) then
       call sgemm('T','N',num_super,num_super,3,-2.0_gp,ion_supercluster,3,ion_supercluster,3,1.0_gp,dij2,num_super)
    else
       call dgemm('T','N',num_super,num_super,3,-2.0_gp,ion_supercluster,3,ion_supercluster,3,1.0_gp,dij2,num_super)
    end if
    
    !$omp parallel do
    do ni=1,num_super
       do nj=1,num_super
          dij(nj,ni)=sqrt(abs(dij2(nj,ni)))
       end do
    end do
    !$omp end parallel do

    if(num_body.gt.5) stop 'frank: only up to five body terms coded'

    allocate(num_feat(num_body),num_comb(num_body))

    ! Indexing composition

    do ni=1,num_super
       nspec(ni:ni)=findloc(comp_names,ion_supernames(ni))
    end do

    ! One body

    if(num_body.ge.1) then

       num_feat(1)=num_comp
       num_comb(1)=1

    end if

    ! Two body

    if(num_body.ge.2) then

       num_feat(2)=npow

       num_comb(2)=num_comp**2

       if(allocated(indx_2b)) deallocate(indx_2b)
       allocate(indx_2b(num_comp,num_comp))

       indx_2b=0

       do ni=1,num_comp

          do nj=1,num_comp

             indx_2b(nj,ni)=1+(ni-1)*num_comp+nj

          end do

       end do

    end if

    ! Three body

    if(num_body.ge.3) then

       if(deltap.gt.0) then
          num_feat(3)=npow**2
       else
          num_feat(3)=npow
       end if

       num_comb(3)=num_comp**2*(num_comp+1)/2

       if(allocated(indx_3b)) deallocate(indx_3b)
       allocate(indx_3b(num_comp,num_comp,num_comp))

       indx_3b=0

       do ni=1,num_comp

          boff=sum(num_comb(1:2))+(ni-1)*num_comp*(num_comp+1)/2

          do nj=1,num_comp

             do nk=1,num_comp

                ! Sort the pair

                nv(1:2)=sort((/nj,nk/))

                mm=0
                do n=1,num_comp
                   do m=n,num_comp
                      mm=mm+1
                      if((n.eq.nv(1)).and.(m.eq.nv(2))) exit
                   end do
                   if((n.eq.nv(1)).and.(m.eq.nv(2))) exit
                end do

                indx_3b(nk,nj,ni)=boff+mm 

             end do

          end do

       end do

    end if

    ! Four body

    if(num_body.ge.4) then

       if(deltap.gt.0) then
          num_feat(4)=npow**2
       else
          num_feat(4)=npow
       end if

       num_comb(4)=num_comp**2*(num_comp+1)*(num_comp+2)/6 

       if(allocated(indx_4b)) deallocate(indx_4b)
       allocate(indx_4b(num_comp,num_comp,num_comp,num_comp))

       indx_4b=0

       do ni=1,num_comp

          boff=sum(num_comb(1:3))+(ni-1)*num_comp*(num_comp+1)*(num_comp+2)/6

          do nj=1,num_comp

             do nk=1,num_comp

                do nl=1,num_comp

                   ! Sort the triple

                   nv(1:3)=sort((/nj,nk,nl/))

                   mm=0
                   do n=1,num_comp 
                      do m=n,num_comp
                         do o=m,num_comp
                            mm=mm+1
                            if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3))) exit
                         end do
                         if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3))) exit
                      end do
                      if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3))) exit
                   end do

                   indx_4b(nl,nk,nj,ni)=boff+mm

                end do
             end do
          end do
       end do

    end if

    ! Five body

    if(num_body.ge.5) then

       if(deltap.gt.0) then
          num_feat(5)=npow**2
       else
          num_feat(5)=npow
       end if

       num_comb(5)=num_comp**2*(num_comp+1)*(num_comp+2)*(num_comp+3)/24

       if(allocated(indx_5b)) deallocate(indx_5b)
       allocate(indx_5b(num_comp,num_comp,num_comp,num_comp,num_comp))

       indx_5b=0

       do ni=1,num_comp

          boff=sum(num_comb(1:4))+(ni-1)*num_comp*(num_comp+1)*(num_comp+2)*(num_comp+3)/24

          do nj=1,num_comp

             do nk=1,num_comp

                do nl=1,num_comp

                   do nm=1,num_comp

                      ! Sort the quadruple

                      nv(1:4)=sort((/nj,nk,nl,nm/))

                      mm=0
                      do n=1,num_comp 
                         do m=n,num_comp
                            do o=m,num_comp
                               do p=o,num_comp
                                  mm=mm+1
                                  if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3)).and.(p.eq.nv(4))) exit
                               end do
                               if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3)).and.(p.eq.nv(4))) exit
                            end do
                            if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3)).and.(p.eq.nv(4))) exit
                         end do
                         if((n.eq.nv(1)).and.(m.eq.nv(2)).and.(o.eq.nv(3)).and.(p.eq.nv(4))) exit
                      end do

                      indx_5b(nm,nl,nk,nj,ni)=boff+mm

                   end do
                end do
             end do
          end do
       end do

    end if
    
    allocate(features(maxval(num_feat),sum(num_comb(1:num_body)),num_ions))

    features=0.0_gp

    elj=0.0_gp

    !$omp parallel do
    do ni=1,num_ions

       do n=1,num_comp
          if(comp_names(n).eq.ion_names(ni)) then
             features(n,1,ni)=1.0_gp
             exit
          end if
       end do

       if(num_body.ge.2) then

          do nj=1,num_super
             if(ni.eq.nj) cycle

             if(dij2(nj,ni).lt.rmax2) then

                n1=indx_2b(nspec(nj),nspec(ni))

                if(n1.gt.0) features(1:num_feat(2),n1,ni)=features(1:num_feat(2),n1,ni)+add_twobody(pow,dij(nj,ni))

                if(lj) elj=elj+(2/dij2(nj,ni))**6-(2/dij2(nj,ni))**3

                if(num_body.ge.3) then

                   do nk=nj+1,num_super
                      if(ni.eq.nk) cycle

                      if((dij2(nk,ni).lt.rmax2).and.(dij2(nk,nj).lt.rmax2)) then

                         n2=indx_3b(nspec(nk),nspec(nj),nspec(ni)) 

                         if(n2.gt.0) features(1:num_feat(3),n2,ni)=features(1:num_feat(3),n2,ni)+&
                              add_threebody(num_feat(3),pow,dij(nj,ni),dij(nk,ni),dij(nj,nk))

                         if(num_body.ge.4) then

                            do nl=nk+1,num_super
                               if(ni.eq.nl) cycle

                               if((dij2(nl,ni).lt.rmax2).and.(dij2(nk,nl).lt.rmax2).and.(dij2(nj,nl).lt.rmax2)) then 

                                  n3=indx_4b(nspec(nl),nspec(nk),nspec(nj),nspec(ni))

                                  if(n3.gt.0) features(1:num_feat(4),n3,ni)=features(1:num_feat(4),n3,ni)+&
                                       add_fourbody(num_feat(4),pow,dij(nj,ni),dij(nk,ni),dij(nl,ni),&
                                       dij(nj,nk),dij(nj,nl),dij(nk,nl))

                                  if(num_body.ge.5) then

                                     do nm=nl+1,num_super
                                        if(ni.eq.nm) cycle

                                        if((dij2(nm,ni).lt.rmax2)&
                                             .and.(dij2(nk,nm).lt.rmax2).and.(dij2(nj,nm).lt.rmax2).and.(dij2(nl,nm).lt.rmax2)) then 

                                           n4=indx_5b(nspec(nm),nspec(nl),nspec(nk),nspec(nj),nspec(ni))

                                           if(n4.gt.0) features(1:num_feat(5),n4,ni)=features(1:num_feat(5),n4,ni)+&
                                                add_fivebody(num_feat(5),pow,dij(nj,ni),dij(nk,ni),dij(nl,ni),dij(nm,ni),&
                                                dij(nj,nk),dij(nj,nl),dij(nk,nl),dij(nm,nj),dij(nm,nk),dij(nm,nl))
                                           
                                        end if

                                     end do

                                  end if

                               end if

                            end do

                         end if

                      end if

                   end do

                end if

             end if

          end do

       end if

    end do
    !$omp end parallel do

    if(lj) write(buff_words(5,1),*) elj

    npows=0
    npows(1:num_body-1)=num_feat(2:num_body)

    write(ctemp,*) comp_names(1:num_comp)    

    ctemp=trim(comp_names(1))
    do n=2,num_comp
       ctemp=trim(ctemp)//'-'//trim(comp_names(n))
    end do

    if(.not.mean) then
       ncentres=num_ions
       cvolume=trim(buff_words(4,1))
       centhalpy=trim(buff_words(5,1))
    else
       ncentres=-num_ions
       read(buff_words(4,1),*) vol
       read(buff_words(5,1),*) enth

       write(cvolume,*) vol
       write(centhalpy,*) enth

    end if

    write (stdout,'(1x,10(a,1x),a,f7.3,a,1x,i0,a,4(1x,i0),a,*(1x,f0.3))') &
         '  structure:',trim(buff_words(2,1)),&
         '  composition:',trim(ctemp),&
         '  pressure:',trim(buff_words(3,1)),&
         '  volume:',trim(cvolume),&
         '  enthalpy:',trim(centhalpy),&
         '  rmax:',rmax,&
         '  centres:',ncentres,&
         '  length:',npows(1:4),&
         '  powers:',pow

    if(.not.mean) then

       do ni=1,num_ions
          if(num_body.ge.1) write (stdout,'(*(3x,f0.14))') features(1:num_feat(1),1,ni)
          if(num_body.ge.2) then
             do n=num_comb(1)+1,num_comb(1)+num_comb(2)
                write (stdout,'(*(3x,f0.14))') features(1:num_feat(2),n,ni)
             end do
          end if
          if(num_body.ge.3) then
             do n=num_comb(1)+num_comb(2)+1,num_comb(1)+num_comb(2)+num_comb(3)
                write (stdout,'(*(3x,f0.14))') features(1:num_feat(3),n,ni)
             end do
          end if
          if(num_body.ge.4) then
             do n=num_comb(1)+num_comb(2)+num_comb(3)+1,num_comb(1)+num_comb(2)+num_comb(3)+num_comb(4)
                write (stdout,'(*(3x,f0.14))') features(1:num_feat(4),n,ni)
             end do
          end if
          if(num_body.ge.5) then
             do n=num_comb(1)+num_comb(2)+num_comb(3)+num_comb(4)+1,num_comb(1)+num_comb(2)+num_comb(3)+num_comb(4)+num_comb(5)
                write (stdout,'(*(3x,f0.14))') features(1:num_feat(5),n,ni)
             end do
          end if
       end do

    else

       if(num_body.ge.1) write (stdout,*) sum(features(1:num_feat(1),1,:),2)/num_ions
       if(num_body.ge.2) then
          do n=num_comb(1)+1,num_comb(1)+num_comb(2)
             write (stdout,'(*(3x,f0.14))') sum(features(1:num_feat(2),n,:),2)/num_ions
          end do
       end if
       if(num_body.ge.3) then
          do n=num_comb(1)+num_comb(2)+1,num_comb(1)+num_comb(2)+num_comb(3)
             write (stdout,'(*(3x,f0.14))') sum(features(1:num_feat(3),n,:),2)/num_ions
          end do
       end if
       if(num_body.ge.4) then
          do n=num_comb(1)+num_comb(2)+num_comb(3)+1,num_comb(1)+num_comb(2)+num_comb(3)+num_comb(4)
             write (stdout,'(*(3x,f0.14))') sum(features(1:num_feat(4),n,:),2)/num_ions
          end do
       end if
       if(num_body.ge.5) then
          do n=num_comb(1)+num_comb(2)+num_comb(3)+num_comb(4)+1,num_comb(1)+num_comb(2)+num_comb(3)+num_comb(4)+num_comb(5)
             write (stdout,'(*(3x,f0.14))') sum(features(1:num_feat(5),n,:),2)/num_ions
          end do
       end if

    end if

  end subroutine generate_features

  function f(r)

    real(kind=gp), intent(in) :: r

    real(kind=gp) :: f

    if(r.gt.rmax) stop 'f(r): r>rmax'
    
    f=2*(1-r/rmax)
    
  end function f

  function add_twobody(p,r) result(fr)

    real(kind=gp), dimension(:), intent(in) :: p
    real(kind=gp), intent(in) :: r

    real(kind=gp), dimension(size(p)) :: fr

    fr=f(r)**p
    
  end function add_twobody

  function add_threebody(nf,p,r1,r2,r3) result(fp)

    integer,                     intent(in) :: nf
    real(kind=gp), dimension(:), intent(in) :: p
    real(kind=gp),               intent(in) :: r1,r2,r3

    real(kind=gp), dimension(nf) :: fp

    real(kind=gp), dimension(3) :: fr

    real(kind=gp) :: fa
    
    integer :: n,n1,n2

    fr(1)=f(r1)
    fr(2)=f(r2)
    fr(3)=f(r3)
    
    fa=product(fr)/norm2(fr) !!! CHECK
     
    if(deltap.gt.0) then

       n=0
       do n1=1,size(p)
          do n2=1,size(p)
             
                n=n+1

                fp(n)=fr(1)**p(n1)*fr(2)**p(n1)*fr(3)**p(n2)
               
          end do
       end do
       
       fp=fp/size(p)

    else

       fp=fa**p
       
    end if

  end function add_threebody

  function add_fourbody(nf,p,r1,r2,r3,r4,r5,r6) result(fp)

    integer,                     intent(in) :: nf
    real(kind=gp), dimension(:), intent(in) :: p
    real(kind=gp),               intent(in) :: r1,r2,r3,r4,r5,r6

    real(kind=gp), dimension(nf) :: fp

    real(kind=gp), dimension(6) :: fr

    real(kind=gp) :: fa
    
    integer :: n,n1,n2

    fr(1)=f(r1)
    fr(2)=f(r2)
    fr(3)=f(r3)
    fr(4)=f(r4)
    fr(5)=f(r5)
    fr(6)=f(r6)
    
    fa=product(fr)/norm2(fr) !!! CHECK
     
    if(deltap.gt.0) then

       n=0
       do n1=1,size(p)
          do n2=1,size(p)
             
                n=n+1

                fp(n)=fr(1)**p(n1)*fr(2)**p(n1)*fr(3)**p(n1)*fr(4)**p(n2)*fr(5)**p(n2)*fr(6)**p(n2)
               
          end do
       end do
       
       fp=fp/size(p)

    else

       fp=fa**p
       
    end if

  end function add_fourbody

  function add_fivebody(nf,p,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10) result(fp)

    integer,                     intent(in) :: nf
    real(kind=gp), dimension(:), intent(in) :: p
    real(kind=gp),               intent(in) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10

    real(kind=gp), dimension(nf) :: fp

    real(kind=gp), dimension(10) :: fr

    real(kind=gp) :: fa
    
    integer :: n,n1,n2

    fr(1)=f(r1)
    fr(2)=f(r2)
    fr(3)=f(r3)
    fr(4)=f(r4)
    fr(5)=f(r5)
    fr(6)=f(r6)
    fr(7)=f(r7)
    fr(8)=f(r8)
    fr(9)=f(r9)
    fr(10)=f(r10)
    
    fa=product(fr)/norm2(fr) !!! CHECK
     
    if(deltap.gt.0) then

       n=0
       do n1=1,size(p)
          do n2=1,size(p)
             
                n=n+1

                fp(n)=fr(1)**p(n1)*fr(2)**p(n1)*fr(3)**p(n1)*fr(4)**p(n1)*&
                     fr(5)**p(n2)*fr(6)**p(n2)*fr(7)**p(n2)*fr(8)**p(n2)*fr(9)**p(n2)*fr(10)**p(n2)
               
          end do
       end do
       
       fp=fp/size(p)

    else

       fp=fa**p
       
    end if

  end function add_fivebody
  
  subroutine get_arguments()

    character(len=80), allocatable,dimension(:) :: cargs

    integer :: i,n,iargc,num_args

    character(len=3) :: ctemp

    num_args = iargc()

    allocate(cargs(num_args))

    do i=1,num_args
       call getarg(i,cargs(i))
    end do

    nomp=1
    if(any(cargs.eq."-ompnp")) then               ! * Number of omp threads (note - not mpi)
       do i=1,num_args-1
          if(cargs(i).eq."-ompnp") exit
       end do
       read(cargs(i+1),*) nomp
       nomp=min(nomp,omp_get_num_procs())
    end if
    call omp_set_num_threads(nomp) 
    
    num_comp=0
    comp_names=''
    if(any(cargs.eq."-c")) then                  ! * Composition space
       do i=1,num_args-1
          if(cargs(i).eq."-c") exit
       end do

       do n=i,num_args-1
          read(cargs(n+1),'(a)') ctemp

          if(index(ctemp,'-').gt.0) exit

          num_comp=num_comp+1

          comp_names(num_comp)=ctemp
          
       end do

    end if

    if(any(cargs.eq."-r")) then                  ! * Maximum radius scaling
       do i=1,num_args-1
          if(cargs(i).eq."-r") exit
       end do
       read(cargs(i+1),*) rmax
    end if

    num_body=2
    if(any(cargs.eq."-nb")) then                 ! * nbody terms
       do i=1,num_args-1
          if(cargs(i).eq."-nb") exit
       end do
       read(cargs(i+1),*) num_body
    end if

    if(any(cargs.eq."-p")) then                  ! * Number of powers
       do i=1,num_args-1
          if(cargs(i).eq."-p") exit
       end do
       read(cargs(i+1),*) npow
    end if

    if(any(cargs.eq."-pmin")) then               ! * Minimum power
       do i=1,num_args-1
          if(cargs(i).eq."-pmin") exit
       end do
       read(cargs(i+1),*) pmin
    end if

    if(any(cargs.eq."-pmax")) then               ! * Maximum power
       do i=1,num_args-1
          if(cargs(i).eq."-pmax") exit
       end do
       read(cargs(i+1),*) pmax
    end if

    deltap=1
    if(any(cargs.eq."-gp")) then                 ! * Delta power
       do i=1,num_args-1
          if(cargs(i).eq."-gp") exit
       end do
       read(cargs(i+1),*) deltap
    end if

    if(any(cargs.eq."-lj")) lj=.true.

    if(any(cargs.eq."-cl")) cluster=.true.

    if(any(cargs.eq."-m")) mean=.true.
    
    if(any(cargs.eq."-h")) goto 100

    deallocate(cargs)

    return

100 write (*,*) 'Usage: frank [-ompnp] [-c] [-r] [-nb] [-p] [-pmin] [-pmax] [-gp] [-lj] [-cl] [-m] [-h] < [.res] '

    write (*,*) '-ompnp : Number of OMP threads'
    write (*,*) '-c     : Composition space'
    write (*,*) '-r     : Maximum radius'
    write (*,*) '-nb    : Number of body terms'    
    write (*,*) '-p     : Powers (number)'    
    write (*,*) '-pmin  : Minimum power (0.1)'    
    write (*,*) '-pmax  : Maximum power (10.0)'   
    write (*,*) '-gp    : Delta power (0)'   
    write (*,*) '-lj    : Lennard-Jones energy'    
    write (*,*) '-cl    : Cluster'    
    write (*,*) '-m     : Mean feature'    
    write (*,*) '-h     : Print this help message'

    stop

  end subroutine get_arguments

  subroutine read_buff(stat)

    integer, intent(out) :: stat

    integer :: n,nw,pos,i

    character(len=max_words*2) :: line

    stat=0
    n=0
    do 
       n=n+1
       if(n.gt.max_lines) stop 'frank: increase num_lines'
       read(5,'(a)',end=100,err=100) buff(n)
       if(buff(n).eq.'END') goto 99
    end do
100 stat=1
    return
99  continue

    num_lines=n

    do n=1,num_lines

       line=detab(buff(n))
       pos = 1 
       nw = 0 
       do 
          i = verify(line(pos:), ' ')  !-- Find next non-blank. 
          if (i == 0) exit             !-- No word found. 
          nw = nw + 1                  !-- Found something.
          pos = pos + i - 1            !-- Move to start of the word. 
          i = scan(line(pos:), ' ')    !-- Find next blank.
          buff_words(nw,n)=line(pos:pos+i-1)
          if (i == 0) exit             !-- No blank found.
          pos = pos + i - 1            !-- Move to the blank. 
       end do
       num_words(n)=nw

    end do

  end subroutine read_buff

  subroutine consolidate()

    integer :: n,ni
    
    if(all(lattice_car.eq.0.0_gp)) call abc2car(lattice_abc,lattice_car)

    if(all(lattice_abc.eq.0.0_gp)) call car2abc(lattice_car,lattice_abc)

    call car2rec(lattice_car,lattice_rec)

    if(all(ion_fractional.eq.0.0_gp)) then
       do ni=1,num_ions
          ion_fractional(:,ni)=matmul(lattice_rec,ion_absolute(:,ni))
       end do
    end if

    do ni=1,num_ions
       ion_fractional(:,ni)=ion_fractional(:,ni)-floor(ion_fractional(:,ni))
    end do

    do ni=1,num_ions
       ion_absolute(:,ni)=matmul(lattice_car,ion_fractional(:,ni))
    end do

    
    if(all(ion_index.eq.0)) then

       if(num_species.eq.0) then
          species_names=''
          do ni=1,num_ions
             if(.not.any(species_names==ion_names(ni))) then
                num_species=num_species+1
                species_names(num_species)=ion_names(ni)
             end if
          end do
       end if
       
       do ni=1,num_ions
          do n=1,num_species
             if(ion_names(ni).eq.species_names(n)) exit
          end do
          ion_index(ni)=n
       end do

       if(num_comp.eq.0) then
          num_comp=num_species
          comp_names=species_names
       end if
       
    end if
    
    lattice_volume=car_volume(lattice_car)

  end subroutine consolidate

  subroutine niggli()

    real(kind=gp) :: lattice_car0(3,3)
    integer       :: Cmat(3,3)

    integer :: ni

    interface
       subroutine niggli_reduce(Lm,Cm)
         use constants
         real(kind=gp), dimension(3,3), intent(inout) :: Lm
         integer,       dimension(3,3), intent(out)   :: Cm
       end subroutine niggli_reduce
    end interface

    lattice_car0=lattice_car
    call niggli_reduce(lattice_car,Cmat)
    call car2rec(lattice_car,lattice_rec)    
    call car2abc(lattice_car,lattice_abc)

    do ni=1,num_ions
       ion_fractional(:,ni)=matmul(lattice_rec,matmul(lattice_car0,ion_fractional(:,ni)))
       ion_absolute(:,ni)=matmul(lattice_car,ion_fractional(:,ni))
    end do

  end subroutine niggli

  ! * =======================================================================

  ! ** SHLX res (result) file format

  subroutine read_res()

    integer :: n,i,ns,ni

    if(.not.(buff_words(1,1).eq.'TITL')) stop 'frank: first line of res/shx data should start with TITL'
    if(.not.(buff_words(1,num_lines).eq.'END')) stop 'frank: last line of res/shx data should start be END'

    do n=1,num_lines
       if(buff_words(1,n).eq.'CELL') then
          do i=1,6
             read(buff_words(i+2,n),*) lattice_abc(i)
          end do
       end if
       if(buff_words(1,n).eq.'SFAC') then
          if(num_species.eq.0) then
             num_species=num_words(n)-1
             do ns=1,num_species
                read(buff_words(ns+1,n),*) species_names(ns)
             end do
          end if
          exit
       end if
    end do

    num_ions=num_lines-n-1

    call init_ions()

    ion_index=0
    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,n+ni))
       !read(buff_words(2,n+ni),*) ion_index(ni)
       read(buff_words(3:5,n+ni),*,err=99) ion_fractional(:,ni)
       read(buff_words(6,n+ni),*,err=99) ion_occ(ni)
       if(num_words(n+ni).gt.6) read(buff_words(7,n+ni),*,err=99) ion_spin(ni)
    end do

    return

99  write (stderr,*) buff_words(1:10,1)
    write (stderr,*) buff_words(3:5,n+ni)
    stop 'read_res: error reading res file'
    
  end subroutine read_res

  subroutine generate_supercluster()

    integer :: n,m,n1,n2,n3,nmax(3)

    real(kind=gp) :: vnm(3),vm(3),rnm2,rmax2,h(3),diag

    logical :: keep
    
    rmax2=rmax**2
    
    h=car2habc(lattice_car)
    
    diag=sqrt(dot_product(sum(lattice_car,2),sum(lattice_car,2)))

    nmax(1) = int((rmax+diag)/h(1))
    nmax(2) = int((rmax+diag)/h(2))
    nmax(3) = int((rmax+diag)/h(3))

    !nmax=10
    
    if(cluster) nmax=0

    num_super=0  
    
    do n1=-nmax(1),nmax(1)
       do n2=-nmax(2),nmax(2)
          do n3=-nmax(3),nmax(3)
             do m=1,num_ions

                vm=ion_absolute(:,m)+n1*lattice_car(:,1)+n2*lattice_car(:,2)+n3*lattice_car(:,3)

                keep=.false.

                do n=1,num_ions

                   vnm=vm-ion_absolute(:,n)
                   rnm2=dot_product(vnm,vnm)

                   if(rnm2.lt.rmax2)  then
                      keep=.true.
                      cycle
                   end if

                end do

                if(keep) num_super=num_super+1

             end do
          end do
       end do
    end do
    
    allocate(ion_supercluster(3,num_super),ion_supernames(num_super),ion_superindex(num_super))

    num_super=0

    ! ** Put the original atoms first
    
    do m=1,num_ions

       num_super=num_super+1

       ion_supercluster(:,num_super)=ion_absolute(:,m)
       ion_supernames(num_super)=ion_names(m)
       ion_superindex(num_super)=ion_index(m)
       
    end do
    
    do n1=-nmax(1),nmax(1)
       do n2=-nmax(2),nmax(2)
          do n3=-nmax(3),nmax(3)

             if((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0)) cycle
             
             do m=1,num_ions

                vm=ion_absolute(:,m)+n1*lattice_car(:,1)+n2*lattice_car(:,2)+n3*lattice_car(:,3)

                keep=.false.

                do n=1,num_ions

                   vnm=vm-ion_absolute(:,n)
                   rnm2=dot_product(vnm,vnm)

                   if(rnm2.lt.rmax2)  then
                      keep=.true.
                      cycle
                   end if

                end do

                if(keep) then
                   num_super=num_super+1
                   ion_supercluster(:,num_super)=vm
                   ion_supernames(num_super)=ion_names(m)
                   ion_superindex(num_super)=ion_index(m)
                end if

             end do
             
          end do
       end do
    end do

  end subroutine generate_supercluster

  !* =======================================================================

  subroutine init_ions()

    if(allocated(ion_names)) stop 'frank: multiple ion definitions'

    allocate(ion_names(num_ions),ion_index(num_ions),ion_occ(num_ions),ion_spin(num_ions))
    allocate(ion_fractional(3,num_ions),ion_absolute(3,num_ions))

    ion_names=''
    ion_index=0
    ion_occ=1.0_gp
    ion_spin=0.0_gp
    ion_fractional=0.0_gp
    ion_absolute=0.0_gp

  end subroutine init_ions

  subroutine car2abc(car,abc)

    real(kind=gp), intent(in)  :: car(3,3)
    real(kind=gp), intent(out) :: abc(6)

    abc(1) = sqrt(car(1,1)**2+car(2,1)**2+car(3,1)**2)
    abc(2) = sqrt(car(1,2)**2+car(2,2)**2+car(3,2)**2)
    abc(3) = sqrt(car(1,3)**2+car(2,3)**2+car(3,3)**2)
    abc(4) = acos(dot_product(car(:,2),car(:,3))/abc(2)/abc(3))/dgrd
    abc(5) = acos(dot_product(car(:,1),car(:,3))/abc(1)/abc(3))/dgrd
    abc(6) = acos(dot_product(car(:,1),car(:,2))/abc(1)/abc(2))/dgrd

  end subroutine car2abc

  subroutine abc2car(abc,car)

    real(kind=gp), intent(in) :: abc(6)
    real(kind=gp), intent(out):: car(3,3)

    car(:,1) = (/abc(1),0.0_gp,0.0_gp/)
    car(:,2) = (/abc(2)*cos(dgrd*abc(6)),abc(2)*sin(dgrd*abc(6)),0.0_gp/)
    car(1,3) = abc(3)*cos(dgrd*abc(5))
    car(2,3) = abc(3)*(cos(dgrd*abc(4))-cos(dgrd*abc(5))*cos(dgrd*abc(6)))/sin(dgrd*abc(6))
    car(3,3) = sqrt(abc(3)**2-car(1,3)**2-car(2,3)**2)

  end subroutine abc2car

  subroutine car2rec(car,rec)

    real(kind=gp), intent(in)  :: car(3,3)
    real(kind=gp), intent(out) :: rec(3,3)

    real(kind=gp) :: volume

    volume = car(1,1)*(car(2,2)*car(3,3)-car(3,2)*car(2,3))+&
         car(2,1)*(car(3,2)*car(1,3)-car(1,2)*car(3,3))+&
         car(3,1)*(car(1,2)*car(2,3)-car(2,2)*car(1,3))

    if(abs(volume).lt.epsilon(1.0_gp)) stop 'frank: zero volume cell detected'

    ! ** Calculate the reciprocal lattice 

    rec(1,1)=car(2,2)*car(3,3)-car(3,2)*car(2,3)
    rec(2,1)=car(2,3)*car(3,1)-car(3,3)*car(2,1)
    rec(3,1)=car(2,1)*car(3,2)-car(3,1)*car(2,2)
    rec(1,2)=car(3,2)*car(1,3)-car(1,2)*car(3,3)
    rec(2,2)=car(3,3)*car(1,1)-car(1,3)*car(3,1)
    rec(3,2)=car(3,1)*car(1,2)-car(1,1)*car(3,2)
    rec(1,3)=car(1,2)*car(2,3)-car(2,2)*car(1,3)
    rec(2,3)=car(1,3)*car(2,1)-car(2,3)*car(1,1)
    rec(3,3)=car(1,1)*car(2,2)-car(2,1)*car(1,2)

    rec(:,:)=rec(:,:)/volume

  end subroutine car2rec

  function car_volume(car)

    real(kind=gp), dimension(3,3), intent(in) :: car

    real(kind=gp) :: car_volume

    car_volume = car(1,1)*(car(2,2)*car(3,3)-car(3,2)*car(2,3))+&
         car(2,1)*(car(3,2)*car(1,3)-car(1,2)*car(3,3))+&
         car(3,1)*(car(1,2)*car(2,3)-car(2,2)*car(1,3))

    if(abs(car_volume).lt.epsilon(1.0_gp)) stop 'frank: zero volume cell detected'

  end function car_volume

  function car2habc(car)

    real(kind=gp), dimension(3,3), intent(in) :: car

    real(kind=gp), dimension(3) :: car2habc

    real(kind=gp), dimension(3) :: acb,acc,ccb
    
    real(kind=gp) :: vol

    vol=car_volume(car)
  
    ccb(1)=car(2,3)*car(3,2)-car(3,3)*car(2,2)
    ccb(2)=car(3,3)*car(1,2)-car(1,3)*car(3,2)
    ccb(3)=car(1,3)*car(2,2)-car(2,3)*car(1,2)

    car2habc(1)=vol/sqrt(dot_product(ccb,ccb))
    
    acc(1)=car(2,1)*car(3,3)-car(3,1)*car(2,3)
    acc(2)=car(3,1)*car(1,3)-car(1,1)*car(3,3)
    acc(3)=car(1,1)*car(2,3)-car(2,1)*car(1,3)

    car2habc(2)=vol/sqrt(dot_product(acc,acc))

    acb(1)=car(2,1)*car(3,2)-car(3,1)*car(2,2)
    acb(2)=car(3,1)*car(1,2)-car(1,1)*car(3,2)
    acb(3)=car(1,1)*car(2,2)-car(2,1)*car(1,2)

    car2habc(3)=vol/sqrt(dot_product(acb,acb))

  end function car2habc

  
  function detab(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=2*max_words)   :: detab

    integer                      :: i

    detab=string
    do while (scan(detab,char(9)).gt.0)
       i=scan(detab,char(9))
       detab=trim(detab(1:i-1))//'      '//trim(detab(i+1:))
    end do

  end function detab

  function sort(list)

    integer, dimension(:), intent(in) :: list

    integer, dimension(size(list)) :: sort
    
    integer :: i,ir,j,l ! Loop counters
    integer :: n,ita

    n=size(list)
    
    sort=list

    if(size(list).lt.2) return

    l=n/2+1
    ir=n

    do
       if(l.gt.1) then
          l=l-1
          ita=sort(l)
       else
          ita=sort(ir)
          sort(ir)=sort(1)
          ir=ir-1
          if(ir.eq.1) then
             sort(1)=ita
             return
          end if
       end if
       i=l
       j=l+l
20     if(j.le.ir) then
          if(j.lt.ir) then
             if(sort(j).lt.sort(j+1)) j=j+1
          end if
          if(ita.lt.sort(j)) then
             sort(i)=sort(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          goto 20
       end if
       sort(i)=ita
    end do
   
  end function sort
 
end program frank

