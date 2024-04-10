!==================================================================================!
!                                    constants                                     !
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
! This module contains important constants                                         !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
module constants

  implicit none

  private

  integer, public, parameter :: sp=4  ! Single = 4 Double = 8 Quad = 16
  integer, public, parameter :: dp=8 
  integer, public, parameter :: qp=16 
  
#ifdef SINGLE
    integer, public, parameter :: gp=sp
#else
    integer, public, parameter :: gp=dp
#endif
  
  
  integer, public, parameter :: stdin=5
  integer, public, parameter :: stdout=6
  integer, public, parameter :: stderr=0

  real(kind=gp), public, parameter :: pi=3.141592653589793238462643383279502884197_gp
  real(kind=gp), public, parameter :: tpi=2.0_gp*pi
  real(kind=gp), public, parameter :: gr=(sqrt(5.0_gp)+1.0_gp)/2.0_gp
  real(kind=gp), public, parameter :: dgrd = pi/180.0_gp
  real(kind=gp), public, parameter :: evbyang3=160.2176487_gp
  real(kind=gp), public, parameter :: bohr2ang = 0.529177210903_gp
  real(kind=gp), public, parameter :: delta = 1e-13_gp
  
  real(kind=gp), public, parameter, dimension(3,3) :: &
       ident=reshape((/1.0_gp,0.0_gp,0.0_gp,0.0_gp,1.0_gp,0.0_gp,0.0_gp,0.0_gp,1.0_gp/),(/3,3/))
  
end module constants
