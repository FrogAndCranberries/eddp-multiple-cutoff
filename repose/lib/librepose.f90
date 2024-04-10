module librepose
  use, intrinsic :: iso_c_binding
  use constants
  use cell
  use ddp

  implicit none

  contains

    subroutine c_set_seedname(c_seed) bind(c, name='c_set_seedname')
      implicit none
      character(len=1), dimension(80), intent(in)  ::  c_seed
      character(len=80)  ::  seed
      integer  ::  i
      ! convert between string types, might be a better way
      seed = ""
      do i=1,80
        if (c_seed(i) == c_null_char) then
          exit
        else
          seed(i:i) = c_seed(i)
        end if
      end do
      ! set seedname in cell module 
      seedname = trim(seed)
    end subroutine

    subroutine c_read_ddp() bind(c, name='c_read_ddp')
      implicit none
      call read_ddp()
    end subroutine

    function c_ddp_rcut() result(ddp_rcut) bind(c, name="c_ddp_rcut")
      implicit none
      real(c_double)  ::  ddp_rcut
      ddp_rcut = rcut
    end function

    subroutine c_eval_ddp(num_ions, ion_names, ion_positions, lattice_car, e, f, s) bind(c, name='c_eval_ddp')

      implicit none

      integer(c_int), value, intent(in)                   ::  num_ions
      character(len=1), dimension(6*num_ions), intent(in) ::  ion_names
      real(c_double), intent(in), dimension(3*num_ions)   ::  ion_positions
      real(c_double), intent(in), dimension(9)            ::  lattice_car
      real(c_double), intent(out)                         ::  e
      real(c_double), intent(out), dimension(3*num_ions)  ::  f
      real(c_double), intent(out), dimension(9)           ::  s

      integer  ::  i,j
      type(structure)  ::  str
      real(dp)  ::  e_fortran
      real(dp), allocatable, dimension(:,:)  ::  f_fortran
      real(dp), dimension(3,3)  ::  s_fortran

      allocate(str%ion_positions(3,num_ions))
      allocate(str%ion_names(num_ions))
      allocate(str%ion_cons(num_ions))
      allocate(f_fortran(3,num_ions))

      ! ----- build input -----
      str%num_ions = num_ions
      do i=1,num_ions  ! set ion_names
        ! convert between string types, might be a better way
        str%ion_names(i) = ""
        do j=1,6
          if (ion_names(6*(i-1)+j) == c_null_char) then
            exit
          else
            str%ion_names(i)(j:j) = ion_names(6*(i-1)+j)
          end if
        end do
      end do
      do j=1,num_ions  ! set ion_positions 
        do i=1,3
          str%ion_positions(i,j) = ion_positions(3*(j-1)+i)
        end do
      end do
      do i=1,3  ! set lattice_car
        do j=1,3
          str%lattice_car(j,i) = lattice_car(3*(i-1)+j)
        end do
      end do
      str%num_symm = 0  ! turn off symmetry
      str%ion_cons = .false.  ! don't constrain ions

      ! ----- evaluate the potential -----
      call eval_ddp(str,e_fortran,f_fortran,s_fortran)

      ! ----- extract output -----
      e = e_fortran
      do j=1,num_ions  ! set forces 
        do i=1,3
          f(3*(j-1)+i) = f_fortran(i,j)
        end do
      end do
      do i=1,3  ! set stress
        do j=1,3
          s(3*(i-1)+j) = s_fortran(j,i)
        end do
      end do

    end subroutine

end module
