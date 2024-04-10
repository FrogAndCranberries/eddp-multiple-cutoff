!==================================================================================!
!                                    data                                          !
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
! This module reads and manages the datasets                                       !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
module data

  use constants
  use nn

  implicit none

  private

  public :: data_read
  public :: data_compatible

  integer, parameter :: unit_data=45

  type, public :: metadata
     integer                                                :: ncomp,nstruct,nlength,ntlength,nflength,nfvlength
     integer,             public, allocatable, dimension(:) :: ncentres
     integer,             public, allocatable, dimension(:) :: natoms
     integer,             public, allocatable, dimension(:) :: struct
     real(kind=gp),       public, allocatable, dimension(:) :: enthalpy,pressure,volume,rmax
     character(len=240),  public, allocatable, dimension(:) :: label
     character(len=240),  public                            :: setname
     character(len=600),  public                            :: composition
     character(len=2000), public                            :: cmnt
  end type metadata

  real(kind=gp), public :: external_pressure
  
contains

  subroutine data_read(datafile,data,meta)

    character(len=*),  intent(in)  :: datafile
    type(dataset),     intent(out) :: data
    type(metadata),    intent(out) :: meta

    ! =======================================

    integer :: ns,nv,nc,ccheck,lcheck,tlcheck,flcheck,fvlcheck,num,n,ncomp,nstart

    character(len=2000) :: ctemp,cvec(19)

    open(unit=unit_data,file=datafile,form='formatted',status='old',err=101)

    read(unit_data,'(a)') ctemp
    meta%cmnt=trim(ctemp(index(ctemp,'composition:'):index(ctemp,'pressure:')-1))
    meta%cmnt=trim(meta%cmnt)//' '//trim(ctemp(index(ctemp,'rmax:'):index(ctemp,'centres:')-1))
    meta%cmnt=trim(meta%cmnt)//' '//trim(ctemp(index(ctemp,'length:'):))
    rewind(unit_data)
    
    meta%setname=datafile

    meta%nstruct=0
    meta%ncomp=0
    meta%nlength=0
    meta%ntlength=0
    meta%nflength=0
    meta%nfvlength=0

    num=0
    do
       read(unit_data,*,end=99,err=99) cvec
               
       if(cvec(1).eq.'structure:')  then
          
          meta%nstruct=meta%nstruct+1
          
          read(cvec(4),*) meta%composition
          
          ncomp=1
          do ns=1,len_trim(meta%composition)
             if(meta%composition(ns:ns).eq.'-') ncomp=ncomp+1
          end do
          
          read(cvec(14),*) ccheck
          if(ccheck.lt.0) ccheck=1
          read(cvec(16),*) lcheck 
          read(cvec(17),*) tlcheck
          read(cvec(18),*) flcheck
          read(cvec(19),*) fvlcheck
          
          if(ncomp.gt.meta%ncomp) meta%ncomp=ncomp
          if(lcheck.gt.meta%nlength) meta%nlength=lcheck
          if(tlcheck.gt.meta%ntlength) meta%ntlength=tlcheck
          if(flcheck.gt.meta%nflength) meta%nflength=flcheck
          if(fvlcheck.gt.meta%nfvlength) meta%nfvlength=fvlcheck
         
          do n=1,ccheck
             read(unit_data,*,end=99,err=99) ctemp
             if(lcheck.gt.0) then
                do ns=1,ncomp**2
                   read(unit_data,*,end=99,err=99) ctemp
                end do
             end if
             if(tlcheck.gt.0) then
                do ns=1,ncomp**2*(ncomp+1)/2
                   read(unit_data,*,end=99,err=99) ctemp
                end do
             end if
             if(flcheck.gt.0) then
                do ns=1,ncomp**2*(ncomp+1)*(ncomp+2)/6
                   read(unit_data,*,end=99,err=99) ctemp
                end do
             end if
             if(fvlcheck.gt.0) then
                do ns=1,ncomp**2*(ncomp+1)*(ncomp+2)*(ncomp+3)/24
                   read(unit_data,*,end=99,err=99) ctemp
                end do
             end if
             
          end do
          num=num+ccheck
       else
          stop 'data_read: problem reading data'
       end if
       
    end do

99  rewind(unit_data)

    data%num=num

    data%idim=meta%ncomp+meta%nlength*meta%ncomp**2
    data%idim=data%idim+meta%ntlength*meta%ncomp**2*(meta%ncomp+1)/2
    data%idim=data%idim+meta%nflength*meta%ncomp**2*(meta%ncomp+1)*(meta%ncomp+2)/6
    data%idim=data%idim+meta%nfvlength*meta%ncomp**2*(meta%ncomp+1)*(meta%ncomp+2)*(meta%ncomp+3)/24

    data%odim=1

    rewind(unit_data)
    
    allocate(data%in(data%idim,data%num),data%out(1,data%num))

    allocate(meta%struct(data%num),meta%ncentres(meta%nstruct),meta%label(meta%nstruct),meta%pressure(meta%nstruct),&
         meta%volume(meta%nstruct),meta%enthalpy(meta%nstruct),meta%rmax(meta%nstruct),meta%natoms(meta%nstruct))

    nv=0
    do ns=1,meta%nstruct

       read(unit_data,*,err=100,end=100) cvec

       read(cvec(2),*)  meta%label(ns)
       read(cvec(6),*)  meta%pressure(ns)
       read(cvec(8),*)  meta%volume(ns)
       read(cvec(10),*) meta%enthalpy(ns)
       meta%enthalpy(ns)=meta%enthalpy(ns)+external_pressure*meta%volume(ns)
       read(cvec(12),*) meta%rmax(ns)

       read(cvec(14),*) ccheck
       read(cvec(16),*) lcheck
       read(cvec(17),*) tlcheck   
       read(cvec(18),*) flcheck   
       read(cvec(19),*) fvlcheck   

       meta%ncentres(ns)=max(1,ccheck)
       meta%natoms(ns)=abs(ccheck)
       
       do nc=1,meta%ncentres(ns)
          nv=nv+1
          
          meta%struct(nv)=ns
          
          read(unit_data,*,err=100,end=100) data%in(1:meta%ncomp,nv)
          
          if(lcheck.gt.0) then
             do n=1,meta%ncomp**2
                
                nstart=meta%ncomp+(n-1)*meta%nlength

                read(unit_data,*,err=100,end=100) data%in(nstart+1:nstart+meta%nlength,nv)
                
             end do
          end if
          
          if(tlcheck.gt.0) then
             do n=1,meta%ncomp**2*(meta%ncomp+1)/2

                nstart=meta%ncomp+meta%ncomp**2*meta%nlength+(n-1)*meta%ntlength

                read(unit_data,*,err=100,end=100) data%in(nstart+1:nstart+meta%ntlength,nv)
                
             end do
          end if

          if(flcheck.gt.0) then
             do n=1,meta%ncomp**2*(meta%ncomp+1)*(meta%ncomp+2)/6
                
                nstart=meta%ncomp+meta%ncomp**2*meta%nlength+meta%ncomp**2*(meta%ncomp+1)/2*meta%ntlength+(n-1)*meta%nflength

                read(unit_data,*,err=100,end=100) data%in(nstart+1:nstart+meta%nflength,nv)
                
             end do
          end if

          if(fvlcheck.gt.0) then
             do n=1,meta%ncomp**2*(meta%ncomp+1)*(meta%ncomp+2)*(meta%ncomp+3)/24
                
                nstart=meta%ncomp+meta%ncomp**2*meta%nlength+meta%ncomp**2*(meta%ncomp+1)/2*meta%ntlength+&
                     meta%ncomp**2*(meta%ncomp+1)*(meta%ncomp+2)/6*meta%nflength+(n-1)*meta%nfvlength

                read(unit_data,*,err=100,end=100) data%in(nstart+1:nstart+meta%nfvlength,nv)
                
             end do
          end if

          data%out(1,nv)=meta%enthalpy(ns)/meta%natoms(ns)

       end do
       
    end do
    
    close(unit_data)

    return

100 stop 'data_read : error reading file'
101 stop 'data_read : file not found'

  end subroutine data_read

  function data_compatible(data1,meta1,data2,meta2) result(compat)

    type(dataset),     intent(in) :: data1
    type(metadata),    intent(in) :: meta1
    type(dataset),     intent(in) :: data2
    type(metadata),    intent(in) :: meta2
    
    logical :: compat
    
    compat=.true.

    if(data1%idim.ne.data2%idim)               compat=.false.
    if(meta1%composition.ne.meta2%composition) compat=.false.
    if(meta1%rmax(1).ne.meta2%rmax(1))         compat=.false.
    if(meta1%ncomp.ne.meta2%ncomp)             compat=.false.
    if(meta1%nlength.ne.meta2%nlength)         compat=.false.
    if(meta1%ntlength.ne.meta2%ntlength)       compat=.false.
    if(meta1%nflength.ne.meta2%nflength)       compat=.false.
    if(meta1%nfvlength.ne.meta2%nfvlength)     compat=.false.
    
  end function data_compatible
  ! -------------------------------------------------------------------

end module data
