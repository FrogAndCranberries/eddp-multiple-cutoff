!==================================================================================!
!                                    constants                                     !
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
! This module contains useful constants                                            !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module constants

  implicit none

  private

  integer, public, parameter :: sp=4 ! Single = 4 Double = 8 Quad = 16
  integer, public, parameter :: dp=8 
  integer, public, parameter :: qp=16
  
#ifdef SINGLE
  integer, public, parameter :: gp=sp
#else
  integer, public, parameter :: gp=dp
#endif
  
  integer, public, parameter :: ki=8
  
  integer, public, parameter :: stdin=5
  integer, public, parameter :: stdout=6
  integer, public, parameter :: stderr=0

  real(kind=gp), public, parameter :: pi=3.141592653589793238462643383279502884197_gp
  real(kind=gp), public, parameter :: tpi=2.0_gp*pi
  real(kind=gp), public, parameter :: gr=(sqrt(5.0_gp)+1.0_gp)/2.0_gp
  real(kind=gp), public, parameter :: dgrd = pi/180.0_gp
  real(kind=gp), public, parameter :: evbyang3=160.2176487_gp
  real(kind=gp), public, parameter :: ev2k=11604.51812_gp
  real(kind=gp), public, parameter :: bohr2ang = 0.529177210903_gp
  real(kind=gp), public, parameter :: ntu2ps = 0.010180509_gp ! sqrt((amu*ang^2)/eV)*1e12
  real(kind=gp), public, parameter :: delta = 1e-13_gp
  

  real(kind=gp), public, parameter, dimension(3,3) :: &
       ident=reshape((/1.0_gp,0.0_gp,0.0_gp,0.0_gp,1.0_gp,0.0_gp,0.0_gp,0.0_gp,1.0_gp/),(/3,3/))
  
  
  character(len=2), public, parameter, dimension(118) :: elements_alpha=(/&
       & 'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne'&
       &,'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca'&
       &,'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn'&
       &,'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr'&
       &,'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn'&
       &,'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd'&
       &,'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb'&
       &,'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg'&
       &,'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th'&
       &,'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'&
       &,'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'&
       &,'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/)

  
  integer, public, parameter, dimension(118) :: elements_valence=(/&
       &  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &, 0, 0, 0, 0, 0, 0, 0 ,0/)
  
  
  ! ** Radii taken from Fig 3 Clementi, E.; Raimond, D. L.; Reinhardt, W. P. (1967).
  !    "Atomic Screening Constants from SCF Functions. II. Atoms with 37 to 86 Electrons".
  !    Journal of Chemical Physics. 47 (4): 1300–1307.
  !    H=0.39 Missing elements set to 2.00 incl. La and Ce

  real(kind=gp), public, parameter, dimension(118) :: elements_radius=(/&
       &  0.39_gp,0.31_gp,1.67_gp,1.12_gp,0.87_gp,0.67_gp,0.56_gp,0.48_gp,0.42_gp,0.38_gp&
       &, 1.90_gp,1.45_gp,1.18_gp,1.11_gp,0.98_gp,0.88_gp,0.79_gp,0.71_gp,2.43_gp,1.94_gp&
       &, 1.84_gp,1.76_gp,1.71_gp,1.66_gp,1.61_gp,1.56_gp,1.52_gp,1.49_gp,1.45_gp,1.42_gp&
       &, 1.36_gp,1.25_gp,1.14_gp,1.03_gp,0.94_gp,0.88_gp,2.65_gp,2.19_gp,2.12_gp,2.06_gp&
       &, 1.98_gp,1.90_gp,1.83_gp,1.78_gp,1.73_gp,1.69_gp,1.65_gp,1.61_gp,1.56_gp,1.45_gp&
       &, 1.33_gp,1.23_gp,1.15_gp,1.08_gp,2.98_gp,2.53_gp,2.00_gp,2.00_gp,2.47_gp,2.06_gp&
       &, 2.05_gp,2.38_gp,2.31_gp,2.33_gp,2.25_gp,2.28_gp,2.26_gp,2.26_gp,2.22_gp,2.22_gp&
       &, 2.17_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp/)

  ! ** van der Waals raii taken from Table 12 of Mantina et al,
  !    "Consistent van der Waals Radii for the Whole Main Group",
  !    J Phys Chem A, 113(19), 5806-5812, 2009
  !    Missing elements set to 2.00

  real(kind=gp), public, parameter, dimension(118) :: elements_vdw=(/&
       &  1.10_gp,1.40_gp,1.81_gp,1.53_gp,1.92_gp,1.70_gp,1.55_gp,1.52_gp,1.47_gp,1.54_gp&
       &, 2.27_gp,1.73_gp,1.84_gp,2.10_gp,1.80_gp,1.80_gp,1.75_gp,1.88_gp,2.75_gp,2.31_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 1.87_gp,2.11_gp,1.85_gp,1.90_gp,1.83_gp,2.02_gp,3.03_gp,2.49_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.06_gp,2.06_gp,1.98_gp,2.16_gp,3.43_gp,2.68_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 1.96_gp,2.02_gp,2.07_gp,1.97_gp,2.02_gp,2.20_gp,3.48_gp,2.93_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp&
       &, 2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp,2.00_gp/)

  ! ** Masses of elements

  real(kind=gp), public, parameter, dimension(118) :: elements_mass=(/&
       & 1.00794_gp,4.00260_gp,6.941_gp,9.012187_gp,10.811_gp,12.0107_gp,14.00674_gp,15.9994_gp,&
       & 18.99840_gp,20.1797_gp,22.98977_gp,24.3050_gp,26.98154_gp,28.0855_gp,30.97376_gp,32.066_gp,&
       & 35.4527_gp,39.948_gp,39.0983_gp,40.078_gp,44.95591_gp,47.867_gp,50.9415_gp,51.9961_gp,&
       & 54.93805_gp,55.845_gp,58.93320_gp,58.6934_gp,63.546_gp,65.39_gp,69.723_gp,72.61_gp,74.92160_gp,&
       & 78.96_gp,79.904_gp,83.80_gp,85.4678_gp,87.62_gp,88.90585_gp,91.224_gp,92.90638_gp,95.94_gp,&
       & 98.0_gp,101.07_gp,102.90550_gp,106.42_gp,107.8682_gp,112.411_gp,114.818_gp,118.710_gp,121.760_gp,&
       & 127.60_gp,126.90447_gp,131.29_gp,132.90545_gp,137.327_gp,138.9055_gp,140.116_gp,140.90765_gp,&
       & 144.24_gp,145.0_gp,150.36_gp,151.964_gp,157.25_gp,158.92534_gp,162.50_gp,164.93032_gp,167.26_gp,&
       & 168.93421_gp,173.04_gp,174.967_gp,178.49_gp,180.9479_gp,183.84_gp,186.207_gp,190.23_gp,192.217_gp,&
       & 195.078_gp,196.96655_gp,200.59_gp,204.3833_gp,207.2_gp,208.98038_gp,209.0_gp,210.0_gp,222.0_gp,&
       & 223.0_gp,226.0_gp,227.0_gp,232.0381_gp,231.03588_gp,238.0289_gp,237.0_gp,244.0_gp,243.0_gp,247.0_gp,&
       & 247.0_gp,251.0_gp,252.0_gp,257.0_gp,258.0_gp,259.0_gp,262.0_gp,261.0_gp,262.0_gp,263.0_gp,264.0_gp,&
       & 265.0_gp,268.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp,0.0_gp/)
  
end module constants
