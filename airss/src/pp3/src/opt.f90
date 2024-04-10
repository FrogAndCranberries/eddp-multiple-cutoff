!==================================================================================!
!                                      opt                                         !
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
! This module performs the geometry optimisation                                   !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module opt

  use constants
  use cell
  use pp
  use ppp

  implicit none

  private

  public :: opt_tpsd

  logical, public :: optimise=.true.
  logical, public :: have_ppp=.false.
  integer, public :: maxsteps=9999

  !---------------------------------------------------!

  ! ** Private data

  real(kind=dp) :: step=1e-4_dp  ! ** Set the intial step size
  real(kind=dp) :: etol=1e-10_dp ! ** Set the energy/enthalpy tolerance

contains

  subroutine opt_tpsd(converged)

    !-------------------------------------------------------------!
    ! Two-point Step Size Gradient Methods - Barzilai and Borwein !
    ! IMA Journal of Numerical Analysis (1988) 8, 141-148         !
    !-------------------------------------------------------------!

    logical, intent(out) :: converged

    real(kind=dp) :: e,h,h0,dh
    real(kind=dp) :: ion_positions0(3,num_ions),lattice_car0(3,3)
    real(kind=dp) :: g(3,num_ions),g0(3,num_ions),sig(3,3),s(3,3),s0(3,3)
    real(kind=dp) :: xg,gg,glength

    integer :: i,j,n,steps

    dh = huge(1.0_dp)
    h  = 1.0_dp
    g  = 1.0_dp
    s  = 1.0_dp
    
    sig= 0.0_dp

    ion_positions0 = ion_positions
    lattice_car0   = lattice_car
    
    steps     = 0
    converged = .true.

    do while(abs(dh).gt.etol)
       
       steps=steps+1

       if(steps.gt.maxsteps) then 
          converged=.false.
          exit
       end if

       h0=h
       g0=g
       s0=s

       if(.not.fix_cell) then
          if(have_ppp) then
             call eval_ppp(e,g,sig)
          else
             call eval_pp(e,g,sig)
          end if
       else
          if(have_ppp) then
             call eval_ppp(e,g)
          else
             call eval_pp(e,g)
          end if
       end if

       s = matmul(sig-external_pressure,lattice_car)

       h = e + (external_pressure(1,1)+external_pressure(2,2)+external_pressure(3,3))/3.0_dp*volume
       
       if(.not.optimise) exit

       xg=0.0_dp
       gg=0.0_dp

       do n=1,num_ions
          do i=1,3
             xg = xg + (ion_positions(i,n)-ion_positions0(i,n))*(g(i,n)-g0(i,n))
             gg = gg + (g(i,n)-g0(i,n))*(g(i,n)-g0(i,n))
          end do
       end do
       do i=1,3
          do j=1,3
             xg = xg + (lattice_car(i,j)-lattice_car0(i,j))*(s(i,j)-s0(i,j))
             gg = gg + (s(i,j)-s0(i,j))*(s(i,j)-s0(i,j))
          end do
       end do

       if(abs(xg).gt.0.0_dp) then
          step = abs(xg/gg)
       else
          glength=0.0_dp
          do n=1,num_ions
             glength=glength+dot_product(g(:,n),g(:,n))
          end do
          glength=sqrt(glength)
          if(glength.gt.1e+8_dp) then
             step=1e-2_dp/sqrt(glength)
          else
             step=1e-8_dp
          end if
       end if
       
       ion_positions0 = ion_positions
       ion_positions  = ion_positions + step*g

       ! * Convert to fractional

       ion_positions0 = matmul(lattice_rec,ion_positions0)
       ion_positions  = matmul(lattice_rec,ion_positions)
       ion_positions0 = ion_positions0 - real(floor(ion_positions,8),dp)
       ion_positions  = ion_positions  - real(floor(ion_positions,8),dp)

       if(.not.fix_cell) then
          lattice_car0   = lattice_car
          lattice_car    = lattice_car + step*s
          call update_cell()
       end if

       ! * Convert back to absolute

       ion_positions0 = matmul(lattice_car,ion_positions0)
       ion_positions = matmul(lattice_car,ion_positions)

       dh=h-h0
       
       if(.not.quiet) then
          write (unit_conv,'(i6,2f65.25)') steps,h,log10(abs(dh))
          call flush(unit_conv)
       end if       
       
    end do

    write (*,'(a11,f20.10)') 'Volume:    ',  volume
    write (*,'(a11,f20.10)') 'Pressure:  ',  (sig(1,1)+sig(2,2)+sig(3,3))/3.0_dp
    write (*,'(a11,f50.10)') 'Energy:    ',  e
    write (*,'(a11,f50.10)') 'Enthalpy:  ',  h
    write (*,'(a11,6f10.5)') 'Stress:    ',  sig(1,1),sig(2,2),sig(3,3),sig(1,2),sig(1,3),sig(2,3)

  end subroutine opt_tpsd

end module opt
