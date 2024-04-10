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
  
  integer, public, parameter :: ki=8
  
  integer, public, parameter :: stdin=5
  integer, public, parameter :: stdout=6
  integer, public, parameter :: stderr=0

  real(kind=dp), public, parameter :: pi=3.141592653589793238462643383279502884197_dp
  real(kind=dp), public, parameter :: tpi=2.0_dp*pi

end module constants
