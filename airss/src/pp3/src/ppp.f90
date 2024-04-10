
!==================================================================================!
!                                      ppp                                         !
!==================================================================================!
!                                                                                  !
! This file is part of the AIRSS structure prediction package.                     !
!                                                                                  !
! AIRSS is free software; you can redistribute it and/or modify it under the terms !
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
! This module reads, knows, and evaluates the polynomial pair potential            !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module ppp

  use constants
  use cell

  implicit none

  private

  public :: read_ppp
  public :: init_ppp
  public :: eval_ppp

  integer, parameter :: unit_ppp=22
  character(len=80)  :: pppfile

  ! ** Data defining the potential

  real(kind=dp)                                   :: range
  
  real(kind=dp),    allocatable, dimension(:)     :: pow
  real(kind=dp),    allocatable, dimension(:,:,:) :: wgt
  real(kind=dp),    allocatable, dimension(:,:)   :: rcut, rcut_sq
  character(len=6), allocatable, dimension(:)     :: spec
 
  integer            :: nspec
  integer            :: npow
  
  !---------------------------------------------------!
  
  ! ** Private variables

  real(kind=dp) :: mod_r12, mod_r12_sq, inv_mod_r12_n, inv_mod_r12_m

contains

  subroutine read_ppp

    integer            :: ns,n,m,i,j
    character(len=200) :: ctemp,ctemp2

    nspec=num_spec

    ! * Allocate ppp arrays

    allocate(rcut(nspec,nspec))
    allocate(rcut_sq(nspec,nspec))

    ! ** Open the ppp file

    write(pppfile,'(a)') trim(seedname)//".ppp"

    open(unit=unit_ppp,file=pppfile,status="old",err=999)

    ! ** Check the species count, and read exponents

    read(unit_ppp,*,err=999) ctemp, range

    read(unit_ppp,'(a)',err=999) ctemp

    ns=words(ctemp)-1
    
    allocate(spec(ns))
    
    read(ctemp,*,err=999) ctemp2,spec

    read(unit_ppp,'(a)',err=999) ctemp

    npow=words(ctemp)-1

    allocate(pow(npow),wgt(npow,ns,ns))

    read(ctemp,*,err=999) ctemp2,pow

    do i=1,ns
       do j=i,ns
          read(unit_ppp,*,err=999) ctemp,ctemp2,wgt(:,i,j)
          wgt(:,j,i)=wgt(:,i,j)
       end do
    end do

    if(ns.lt.nspec) stop 'PPP is for the wrong number of species in read_ppp'
  
    close(unit_ppp)
    
    return

999 stop 'There is a problem reading the ppp information. Stopping.'

  end subroutine read_ppp

  subroutine init_ppp

    rcut = range
    
    rcut_sq=rcut**2

  end subroutine init_ppp

  subroutine eval_ppp(e,f,s)

    real(kind=dp),           intent(out) :: e
    real(kind=dp), optional, intent(out) :: f(3,num_ions)
    real(kind=dp), optional, intent(out) :: s(3,3)

    integer       :: ion_i,ion_j,ispeci,ispecj
    integer       :: na,nna,nb,nnb,nc,nnc,ni,ns,nn,nnmx,nn0,np
    
    real(kind=dp) :: r12(3),rij(3),f12(3)
    real(kind=dp) :: de,rcsq,rcmod,rijsq,rijmod,mod_f12,work(3,num_ions*num_symm),swork(3,3)
    real(kind=dp), allocatable :: lvec(:,:),ll(:)

    logical       :: sameion
    
    e=0.0_dp
    f12=0.0_dp
    if(present(f)) f=0.0_dp 
    if(present(s)) s=0.0_dp

    ! Determine the supercell required

    nna = (int(2*maxval(rcut)/lattice_abc(1))+1)/2+1
    nnb = (int(2*maxval(rcut)/lattice_abc(2))+1)/2+1
    nnc = (int(2*maxval(rcut)/lattice_abc(3))+1)/2+1
    
    if(cluster) then
       nna=0 ; nnb=0 ; nnc=0
    end if

    nnmx=(2*nna+1)*(2*nnb+1)*(2*nnc+1)
    allocate(lvec(3,nnmx),ll(nnmx))

    nn=0
    nn0=0
    do na=-nna,nna
       do nb=-nnb,nnb
          do nc=-nnc,nnc
             nn=nn+1
             lvec(1:3,nn)=na*lattice_car(1:3,1)+nb*lattice_car(1:3,2)+nc*lattice_car(1:3,3)
             ll(nn)=sqrt(lvec(1,nn)*lvec(1,nn)+lvec(2,nn)*lvec(2,nn)+lvec(3,nn)*lvec(3,nn))
             if(((na==0).and.(nb==0).and.(nc==0))) nn0=nn
          end do
       end do
    end do

    ! ** PPP calculation of energy, forces and stress

    do ion_i=1,num_ions

       ispeci=ion_species(ion_i)

       do ion_j=ion_i,num_ions

          ispecj=ion_species(ion_j)

          rij(:) = ion_positions(:,ion_i)-ion_positions(:,ion_j)

          rijsq = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
          rijmod=sqrt(rijsq)

          sameion=ion_i.eq.ion_j

          rcsq=rcut_sq(ispeci,ispecj)
          rcmod=sqrt(rcsq)

          do nn=1,nnmx

             if((nn.eq.nn0).and.sameion) cycle

             if(abs(rijmod-ll(nn)).gt.rcmod) cycle

             r12(1:3) = rij(1:3)-lvec(1:3,nn)

             mod_r12_sq=r12(1)*r12(1)+r12(2)*r12(2)+r12(3)*r12(3)

             ! Check separation is not too big

             if (mod_r12_sq <= rcsq) then

                mod_r12=sqrt(mod_r12_sq)

                ! Calculate PPP energy, forces and stress

!!$                inv_mod_r12_n=(sigma(ispeci,ispecj)/mod_r12)**PPP_n(ispeci,ispecj)*beta(ispeci,ispecj)
!!$                inv_mod_r12_m=(sigma(ispeci,ispecj)/mod_r12)**PPP_m(ispeci,ispecj)

!!$                de=2.0_dp*epsil(ispeci,ispecj)*(inv_mod_r12_m-inv_mod_r12_n)

                de=0.0_dp
                do np=1,npow
                   de=de+2.0_dp*wgt(np,ispeci,ispecj)/mod_r12**pow(np)
                end do
                
                if(present(f)) then

!!$                   mod_f12=2.0_dp*epsil(ispeci,ispecj)*&
!!$                        (real(PPP_m(ispeci,ispecj),dp)*inv_mod_r12_m/mod_r12-&
!!$                        real(PPP_n(ispeci,ispecj),dp)*inv_mod_r12_n/mod_r12)

                   mod_f12=0.0_dp
                   do np=1,npow
                      mod_f12=mod_f12+2.0_dp*wgt(np,ispeci,ispecj)*pow(np)/mod_r12**(pow(np)+1)
                   end do
                   
                   f12(:)=mod_f12*r12(:)/mod_r12
                   
                end if

                if(ion_j.ne.ion_i) then
                   if(present(f)) f12=f12*2.0_dp
                   de=de*2.0_dp
                end if

                ! Update energy in model

                e=e+de

                ! Update forces in model

                if(present(f)) then
                   f(:,ion_i)=f(:,ion_i)+f12(:)
                   if(ion_j.ne.ion_i) then                         
                      f(:,ion_j)=f(:,ion_j)-f12(:)
                   end if
                end if

                ! Update stress

                if(present(s)) then
                   s(1,1)=s(1,1)+f12(1)*r12(1)
                   s(2,2)=s(2,2)+f12(2)*r12(2)
                   s(3,3)=s(3,3)+f12(3)*r12(3)
                   s(2,3)=s(2,3)+(f12(2)*r12(3)+f12(3)*r12(2))/2.0_dp
                   s(3,1)=s(3,1)+(f12(3)*r12(1)+f12(1)*r12(3))/2.0_dp
                   s(1,2)=s(1,2)+(f12(1)*r12(2)+f12(2)*r12(1))/2.0_dp
                end if

             end if


          end do


       end do  !ion_j
    end do !ion_i

    deallocate(lvec)
    
    if(present(s)) then
       s(3,2) = s(2,3)
       s(1,3) = s(3,1)
       s(2,1) = s(1,2)
       s=s/volume       ! * Convert stress to absolute form
    end if

    ! * Balance forces

    if(present(f)) then
       f12(1) = sum(f(1,:))/real(num_ions,dp)
       f12(2) = sum(f(2,:))/real(num_ions,dp)
       f12(3) = sum(f(3,:))/real(num_ions,dp)
       
       f(1,:) = f(1,:)-f12(1)
       f(2,:) = f(2,:)-f12(2)
       f(3,:) = f(3,:)-f12(3)
    end if

    if(num_symm.gt.1) then
       
       if(present(f)) then
          
          ! * Symmetrise forces
          
          work=f
          f=0.0_dp
          
          do ni=1,num_ions
             do ns=1,num_symm
                f(:,ni) = f(:,ni) +  matmul(symm_ops(1:3,1:3,ns),work(:,ion_equiv(ns,ni))) 
             end do
             f(:,ni) = f(:,ni)/real(num_symm,dp)
          end do
          
       end if
       
       if(present(s)) then

          ! * Symmetrise stresses

          swork = s/real(num_symm,dp)
          s=0.0_dp
          
          do ns=1,num_symm
             s = s + matmul(symm_ops(1:3,1:3,ns),matmul(swork,transpose(symm_ops(1:3,1:3,ns))))
          end do

       end if

    end if

    ! ** Constrain ions

    do ni=1,num_ions
       if(ion_cons(ni)) f(:,ni)=0.0_dp
    end do
    

  end subroutine eval_ppp

  ! -----------------------------------------------------------

  function words(line) result(nw)

    character(len=*) :: line

    character(len=len(line)) :: lwork
    
    integer :: nw,n
    
    integer :: indx
    
    nw=0
    lwork=trim(adjustl(detab(line)))
    do while (len_trim(lwork).gt.0)
       nw=nw+1
       indx=index(lwork,' ')
       lwork=adjustl(lwork(indx+1:))
    end do
    
  end function words

  function detab(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: detab

    integer                      :: i

    detab=string
    do while (scan(detab,char(9)).gt.0)
       i=scan(detab,char(9))
       detab=trim(detab(1:i-1))//'      '//trim(detab(i+1:))
    end do

  end function detab
  
end module ppp
