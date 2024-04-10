!==================================================================================!
!                                      cell                                        !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This module read, knows and writes the unit cell                                 !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module cell

  use constants
  use omp_lib
  use rng
  use spglib_f08

  implicit none

  private

  public :: read_cell
  public :: write_cell
  public :: write_xyze
  public :: update_cell
  public :: findsymm_cell
  public :: shake_nearly_cubic_supercell
  public :: random_supercell
  public :: supercell
  public :: cell_volume
  public :: cell_rec
  public :: car2habc
  public :: deunder

  integer, parameter, public :: unit_xyze=31
  integer, parameter, public :: unit_xyze_run=32

  integer, public :: scell_attempts=1000000

  ! ** Public data defining the system

  type, public :: structure
     integer                                      :: num_ions
     integer                                      :: num_spec
     integer                                      :: num_symm
     integer                                      :: num_cells     
     integer                                      :: num_path
     integer                                      :: supercell_matrix(3,3)
     integer,       dimension(:), allocatable     :: ion_species
     integer,       dimension(:,:), allocatable   :: ion_equiv
     real(kind=gp)                                :: lattice_abc(6)
     real(kind=gp)                                :: lattice_car(3,3)
     real(kind=gp)                                :: lattice_rec(3,3)
     real(kind=gp)                                :: external_pressure(3,3)
     real(kind=gp)                                :: volume
     real(kind=gp)                                :: energy
     real(kind=gp)                                :: enthalpy
     real(kind=gp)                                :: stress(3,3)
     real(kind=gp)                                :: variance
     real(kind=gp), dimension(:,:,:), allocatable :: symm_ops
     real(kind=gp), dimension(:,:), allocatable   :: ion_positions
     real(kind=gp), dimension(:,:), allocatable   :: external_forces
     real(kind=gp), dimension(:,:), allocatable   :: subcell_vec
     real(kind=gp), dimension(:), allocatable     :: ion_rad
     real(kind=gp), dimension(:), allocatable     :: ion_vdw
     real(kind=gp), dimension(:), allocatable     :: ion_spins
     real(kind=gp), dimension(:), allocatable     :: ion_mass
     real(kind=gp), dimension(:,:), allocatable   :: kpoint_path
     character(len=1), dimension(:), allocatable  :: path_label
     logical,       dimension(:), allocatable     :: ion_cons
     character(len=11)                            :: spacegroup='n/a'
     character(len=6), dimension(:), allocatable  :: ion_names
     character(len=6), dimension(:), allocatable  :: ion_names_spec  
  end type structure

  real(kind=gp), public :: poveride=-9999999.99_gp,voveride=0.0_gp

  real(kind=gp), public :: symmtol=1e-3_gp

  character(len=80), public :: seedname
  
  integer, public :: track_every
  integer, public :: total_steps=0
  
  logical, public :: track=.false.
  logical, public :: fix_cell=.false.
  logical, public :: fix_vol=.false.
  logical, public :: cluster=.false.
  logical, public :: quiet=.false.
  logical, public :: symmgen=.false.

  !---------------------------------------------------!

  ! ** Private IO data

  integer, parameter                           :: unit_cell=20
  integer, parameter                           :: unit_out_cell=21
  character(len=80)                            :: cellfile

  ! ** Counters and temp variables

  integer :: i,j,m,n,ni,num_pot,num_hub,ns,na,nb,nc,nmax

  real(kind=gp) :: ang,scale,v(3),x,y,z

  character(len=120)                           :: ctemp,ctemp2 
  character(len=80)                            :: kpspacing
  character(len=80), dimension(:), allocatable :: spec_pot
  character(len=80), dimension(:), allocatable :: hubbard_u

contains

  subroutine read_cell(str)

    type(structure), intent(inout) :: str

    logical :: ionabs

    ! ** Open the cell file

    write(cellfile,'(a)') trim(seedname)//".cell"

    open(unit=unit_cell,file=cellfile,status="old",err=999)

    str%lattice_abc = 0.0_gp

    ! ** Read the lattice vectors

    do
       read(unit_cell,'(a)',end=100) ctemp

       if(ctemp(1:1).eq.'#') cycle

       if((index(ctemp,'%BLOCK LATTICE_ABC')>0).or.(index(ctemp,'%BLOCK lattice_abc')>0)) then

          read(unit_cell,*) str%lattice_abc(1:3)
          read(unit_cell,*) str%lattice_abc(4:6)

          ang = 1.0_gp-cos(dgrd*str%lattice_abc(4))**2-cos(dgrd*str%lattice_abc(5))**2-cos(dgrd*str%lattice_abc(6))**2+&
               2.0_gp*cos(dgrd*str%lattice_abc(4))*cos(dgrd*str%lattice_abc(5))*cos(dgrd*str%lattice_abc(6))

          str%volume = abs(str%lattice_abc(1)*str%lattice_abc(2)*str%lattice_abc(3)*sqrt(abs(ang)))

       end if
       if(trim(ctemp)=='FIX_ALL_CELL : true') fix_cell=.true.
       if(trim(ctemp)=='FIX_VOL : true') fix_vol=.true.
       if(trim(ctemp)=='SYMMETRY_GENERATE') symmgen=.true.
       if(trim(ctemp)=='#CLUSTER') then
          cluster=.true.
          fix_cell=.true.
       end if
    end do

100 rewind(unit_cell)

    do
       read(unit_cell,'(a)',end=101) ctemp
       if((index(ctemp,'%BLOCK LATTICE_CART')>0).or.(index(ctemp,'%BLOCK lattice_cart')>0)) then

          read(unit_cell,*) str%lattice_car(1:3,1)
          read(unit_cell,*) str%lattice_car(1:3,2)
          read(unit_cell,*) str%lattice_car(1:3,3)

          str%volume = cell_volume(str%lattice_car)

          read(unit_cell,*) ctemp

       end if

    end do

101 rewind(unit_cell)    


    ! ** Read the cell contents

    str%num_ions=0
    ionabs = .false.

    do
       read(unit_cell,'(a)',end=103) ctemp
       if((index(ctemp,'%BLOCK POSITIONS_FRAC')>0).or.(index(ctemp,'%BLOCK positions_frac')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK POSITIONS_FRAC')>0).or.(index(ctemp,'%ENDBLOCK positions_frac')>0)) exit
             str%num_ions=str%num_ions+1
          end do

          exit
       end if

    end do

103 rewind(unit_cell)

    if(str%num_ions==0) then

       do
          read(unit_cell,'(a)',end=999) ctemp
          if((index(ctemp,'%BLOCK POSITIONS_ABS')>0).or.(index(ctemp,'%BLOCK positions_abs')>0)) then
             do
                read(unit_cell,'(a)') ctemp
                if((index(ctemp,'%ENDBLOCK POSITIONS_ABS')>0).or.(index(ctemp,'%ENDBLOCK positions_abs')>0)) exit
                str%num_ions=str%num_ions+1
             end do

             exit
          end if

       end do

       ionabs = .true.
       rewind(unit_cell)

    end if

    if(allocated(str%ion_positions)) deallocate(str%ion_positions,str%ion_names,str%ion_species,&
         str%ion_equiv,str%ion_names_spec,str%external_forces,str%ion_spins)
    allocate(str%ion_positions(3,str%num_ions),str%external_forces(3,str%num_ions),str%ion_cons(str%num_ions))
    allocate(str%ion_names(str%num_ions),str%ion_names_spec(str%num_ions))
    allocate(str%ion_species(str%num_ions),str%ion_spins(str%num_ions))

    str%external_forces=0.0_gp
    str%ion_spins=0.0_gp
    do
       read(unit_cell,'(a)',end=110) ctemp
       if((index(ctemp,'%BLOCK POSITIONS_')>0).or.(index(ctemp,'%BLOCK positions_')>0)) then

          do i=1,str%num_ions
             read(unit_cell,'(a)') ctemp
             if(index(ctemp,'#').gt.0) then
                if(index(ctemp,'SPIN=').gt.0) then
                   read(ctemp(1:index(ctemp,'SPIN=')),*) str%ion_names(i),str%ion_positions(:,i)
                   read(ctemp(index(ctemp,'SPIN=')+5:),*) str%ion_spins(i)
                else
                   read(ctemp(1:index(ctemp,'#')),*) str%ion_names(i),str%ion_positions(:,i)
                end if
                if(index(ctemp,'F=').gt.0) then
                   read(ctemp(index(ctemp,'F=')+2:),*) str%external_forces(:,i)
                end if
             else
                if(index(ctemp,'SPIN=').gt.0) then
                   read(ctemp(1:index(ctemp,'SPIN=')),*) str%ion_names(i),str%ion_positions(:,i)
                   read(ctemp(index(ctemp,'SPIN=')+5:),*) str%ion_spins(i)
                else
                   read(ctemp,*) str%ion_names(i),str%ion_positions(:,i)
                end if
             end if
          end do

          exit
       end if

    end do

110 rewind(unit=unit_cell)

    ! ** Read pseudopotentials

    num_pot = 0

    do
       read(unit_cell,'(a)',end=106) ctemp
       if((index(ctemp,'%BLOCK SPECIES_POT')>0).or.(index(ctemp,'%BLOCK species_pot')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK SPECIES_POT')>0).or.(index(ctemp,'%ENDBLOCK species_pot')>0)) exit
             num_pot=num_pot+1
          end do
          exit
       end if

    end do

106 rewind(unit=unit_cell)    

    if(allocated(spec_pot)) deallocate(spec_pot)
    allocate(spec_pot(num_pot))

    do
       read(unit_cell,'(a)',end=107) ctemp
       if((index(ctemp,'%BLOCK SPECIES_POT')>0).or.(index(ctemp,'%BLOCK species_pot')>0)) then
          do i=1,num_pot
             read(unit_cell,'(a)',end=999) spec_pot(i)
          end do
          exit
       end if
    end do
107 rewind(unit=unit_cell)

    ! ** Read Hubbard U

    num_hub = 0

    do
       read(unit_cell,'(a)',end=113) ctemp
       if((index(ctemp,'%BLOCK HUBBARD_U')>0).or.(index(ctemp,'%BLOCK hubbard_u')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK HUBBARD_U')>0).or.(index(ctemp,'%ENDBLOCK hubbard_u')>0)) exit
             num_hub=num_hub+1
          end do
          exit
       end if

    end do

113 rewind(unit=unit_cell)   

    if(allocated(hubbard_u)) deallocate(hubbard_u)
    allocate(hubbard_u(num_hub))

    do
       read(unit_cell,'(a)',end=114) ctemp
       if((index(ctemp,'%BLOCK HUBBARD_U')>0).or.(index(ctemp,'%BLOCK hubbard_u')>0)) then
          do i=1,num_hub
             read(unit_cell,'(a)',end=999) hubbard_u(i)
          end do
          exit
       end if
    end do
114 rewind(unit=unit_cell)    

    ! ** Read fine phonon kpoint path

    str%num_path=0
    do
       read(unit_cell,'(a)',end=120) ctemp
       if((index(ctemp,'%BLOCK PHONON_FINE_KPOINT_PATH')>0).or.(index(ctemp,'%BLOCK phonon_fine_kpoint_path')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK PHONON_FINE_KPOINT_PATH')>0).or.(index(ctemp,'%ENDBLOCK phonon_fine_kpoint_path')>0)) exit
             str%num_path=str%num_path+1
          end do

          exit
       end if
    end do
120 rewind(unit=unit_cell)

    allocate(str%kpoint_path(3,str%num_path),str%path_label(str%num_path))
    
    n=0
    do
       read(unit_cell,'(a)',end=121) ctemp
       if((index(ctemp,'%BLOCK PHONON_FINE_KPOINT_PATH')>0).or.(index(ctemp,'%BLOCK phonon_fine_kpoint_path')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK PHONON_FINE_KPOINT_PATH')>0).or.(index(ctemp,'%ENDBLOCK phonon_fine_kpoint_path')>0)) exit
             n=n+1
             read(ctemp,*) str%kpoint_path(:,n),str%path_label(n)
          end do

          exit
       end if
    end do
121 rewind(unit=unit_cell)

!!$    write (stderr,*) str%num_path
!!$    write (stderr,'(3f10.5)') str%kpoint_path(:,:)
!!$    write (stderr,*) str%path_label
!!$    
!!$    stop
    
    ! ** Read other info

    kpspacing = ' '

    do
       read(unit_cell,'(a)',end=108) ctemp
       if(index(ctemp,'KPOINTS_MP_SPACING')>0) kpspacing=trim(ctemp)

    end do
108 rewind(unit=unit_cell)   

    str%external_pressure=0.0_gp
    do
       read(unit_cell,'(a)',end=109) ctemp
       if((index(ctemp,'%BLOCK EXTERNAL_PRESSURE')>0).or.(index(ctemp,'%BLOCK external_pressure')>0)) then
          read(unit_cell,*,end=999) str%external_pressure(1,1),str%external_pressure(1,2),str%external_pressure(1,3)
          read(unit_cell,*,end=999) str%external_pressure(2,2),str%external_pressure(2,3)
          read(unit_cell,*,end=999) str%external_pressure(3,3)
          str%external_pressure(2,1)=str%external_pressure(1,2)
          str%external_pressure(3,1)=str%external_pressure(1,3)
          str%external_pressure(3,2)=str%external_pressure(2,3)
          exit
       end if
    end do
109 rewind(unit=unit_cell)
    str%external_pressure=str%external_pressure/160.21766208_gp ! Convert to eV/ang^3

    str%supercell_matrix=0
    do
       read(unit_cell,'(a)',end=122) ctemp
       if((index(ctemp,'%BLOCK PHONON_SUPERCELL_MATRIX')>0).or.(index(ctemp,'%BLOCK phonon_supercell_matrix')>0)&
            .or.(index(ctemp,'%BLOCK SUPERCELL_MATRIX')>0).or.(index(ctemp,'%BLOCK supercell_matrix')>0)) then
          read(unit_cell,*,end=999) str%supercell_matrix(1,1),str%supercell_matrix(2,1),str%supercell_matrix(3,1)
          read(unit_cell,*,end=999) str%supercell_matrix(1,2),str%supercell_matrix(2,2),str%supercell_matrix(3,2)
          read(unit_cell,*,end=999) str%supercell_matrix(1,3),str%supercell_matrix(2,3),str%supercell_matrix(3,3)          
          exit
       end if
    end do
122 rewind(unit=unit_cell)
    
    str%num_symm=0
    allocate(str%symm_ops(3,4,48))
    str%symm_ops=0.0_gp
    do
       read(unit_cell,'(a)',end=111) ctemp
       if((index(ctemp,'%BLOCK SYMMETRY_OPS')>0).or.(index(ctemp,'%BLOCK symmetry_ops')>0)) then
          do
             str%num_symm=str%num_symm+1
             read(unit_cell,*,err=99,end=999) &
                  str%symm_ops(1,1,str%num_symm),str%symm_ops(2,1,str%num_symm),str%symm_ops(3,1,str%num_symm)
             read(unit_cell,*,err=99,end=999) &
                  str%symm_ops(1,2,str%num_symm),str%symm_ops(2,2,str%num_symm),str%symm_ops(3,2,str%num_symm)
             read(unit_cell,*,err=99,end=999) &
                  str%symm_ops(1,3,str%num_symm),str%symm_ops(2,3,str%num_symm),str%symm_ops(3,3,str%num_symm)
             read(unit_cell,*,err=99,end=999) &
                  str%symm_ops(1,4,str%num_symm),str%symm_ops(2,4,str%num_symm),str%symm_ops(3,4,str%num_symm)
          end do
99        str%num_symm=str%num_symm-1
          exit
       end if
    end do
111 rewind(unit=unit_cell)   


    ! ** Read the ionic constraints

    str%ion_cons=.false.
    do
       read(unit_cell,'(a)',end=112) ctemp
       if((index(ctemp,'%BLOCK IONIC_CONSTRAINTS')>0).or.(index(ctemp,'%BLOCK ionic_constraints')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK IONIC_CONSTRAINTS')>0).or.(index(ctemp,'%ENDBLOCK ionic_constraints')>0)) exit
             read(ctemp,*) i,ctemp2,j,x,y,z
             n=0
             do ni=1,str%num_ions
                if(ctemp2.eq.str%ion_names(ni)) n=n+1
                if(n.eq.j) exit
             end do
             str%ion_cons(ni)=.true.

          end do
          exit
       end if

    end do
112 rewind(unit=unit_cell)

    if(poveride.gt.-1000.00_gp) then
       str%external_pressure=0.0_gp
       str%external_pressure(1,1)=poveride
       str%external_pressure(2,2)=poveride
       str%external_pressure(3,3)=poveride
    end if

    if(.not.all(str%lattice_abc==0.0_gp)) then

       str%lattice_car=cell_car(str%lattice_abc)

    else

       str%lattice_abc=cell_abc(str%lattice_car)

    end if

    str%lattice_rec=cell_rec(str%lattice_car)

    ! ** Convert to absolute coordinates

    if(.not.ionabs) then
       do n=1,str%num_ions
          str%ion_positions(:,n) = str%ion_positions(:,n) - real(floor(str%ion_positions(:,n),8),gp)
          str%ion_positions(:,n) = str%ion_positions(1,n)*str%lattice_car(:,1) + str%ion_positions(2,n)*str%lattice_car(:,2) &
               + str%ion_positions(3,n)*str%lattice_car(:,3)
       end do
    else
       do n=1,str%num_ions
          str%ion_positions(:,n) = matmul(str%lattice_rec,str%ion_positions(:,n)) 
          str%ion_positions(:,n) = str%ion_positions(:,n) - real(floor(str%ion_positions(:,n),8),gp)
          str%ion_positions(:,n) = str%ion_positions(1,n)*str%lattice_car(:,1) + str%ion_positions(2,n)*str%lattice_car(:,2) &
               + str%ion_positions(3,n)*str%lattice_car(:,3)
       end do
    end if

    ! ** Scale to a new str%volume specified on command line

    if(voveride>0.0_gp) then

       scale = (voveride/str%volume)**(1.0_gp/3.0_gp)

       str%lattice_abc(1:3) = str%lattice_abc(1:3)*scale
       str%lattice_car      = str%lattice_car*scale
       str%lattice_rec      = str%lattice_rec/scale
       str%ion_positions    = str%ion_positions*scale
       str%volume           = str%volume*scale**3

    end if

    ! ** How many distinct "species"

    str%num_spec = 0
    do n=1,str%num_ions
       if(.not.any(str%ion_names(1:n-1)==str%ion_names(n))) then
          str%num_spec = str%num_spec+1
          str%ion_names_spec(str%num_spec) = str%ion_names(n)
       end if

    end do

    ! ** Set the str%ion_species array correctly

    do n=1,str%num_ions
       do m=1,str%num_spec
          if(str%ion_names(n)==str%ion_names_spec(m)) then
             str%ion_species(n) = m
             exit
          end if
       end do
    end do

    ! ** Set radii

    allocate(str%ion_rad(str%num_ions),str%ion_vdw(str%num_ions),str%ion_mass(str%num_ions))

    do ni=1,str%num_ions
       do i=1,118
          if(str%ion_names(ni).eq.elements_alpha(i)) then
             str%ion_rad(ni)=elements_radius(i)
             str%ion_vdw(ni)=elements_vdw(i)
             str%ion_mass(ni)=elements_mass(i)
             exit
          end if
       end do
    end do

    if(symmgen.and.(.not.any(str%ion_cons))) call findsymm_cell(str)

    ! ** Identify symmetry equivalent ions

    if(str%num_symm.gt.1) then

       allocate(str%ion_equiv(str%num_symm,str%num_ions))

       nmax=1

88     str%ion_equiv=-1

       do n=1,str%num_ions
          do m=1,str%num_ions

             do ns=1,str%num_symm

                do na=-nmax,nmax
                   do nb=-nmax,nmax
                      do nc=-nmax,nmax

                         v(:) = str%ion_positions(:,n) - matmul(str%symm_ops(1:3,1:3,ns),str%ion_positions(:,m)) - &
                              str%symm_ops(1,4,ns)*str%lattice_car(1:3,1) - &
                              str%symm_ops(2,4,ns)*str%lattice_car(1:3,2) - &
                              str%symm_ops(3,4,ns)*str%lattice_car(1:3,3) - &
                              na*str%lattice_car(:,1)-nb*str%lattice_car(:,2)-nc*str%lattice_car(:,3)

                         if(dot_product(v,v).lt.1e-4_gp) str%ion_equiv(ns,n)=m

                      end do
                   end do
                end do

             end do

          end do
       end do

       if (any(str%ion_equiv.lt.0)) then
          nmax=nmax+1
          if(nmax.gt.10) stop 'read_cell : failed to find str%ion_equiv'
          goto 88
       end if

    end if

    close(unit_cell)

    return

999 stop 'There is a problem reading the cell information. Stopping.'

  end subroutine read_cell

  subroutine write_cell(str)

    type(structure), intent(in) :: str

    write(cellfile,'(a)') trim(seedname)//"-out.cell"

    open(unit=unit_out_cell,file=cellfile,status="unknown",err=998)

    ! ** Write the lattice vectors

    write (unit_out_cell,'(a)') "%BLOCK LATTICE_ABC"
    write (unit_out_cell,'(3f20.15)') str%lattice_abc(1:3)
    write (unit_out_cell,'(3f20.15)') str%lattice_abc(4:6)
    write (unit_out_cell,'(a)') "%ENDBLOCK LATTICE_ABC"
    write (unit_out_cell,*)

    ! ** Write the cell contents

    write (unit_out_cell,'(a)') "%BLOCK POSITIONS_FRAC"
    do i=1,str%num_ions
       v(:) = matmul(str%lattice_rec,str%ion_positions(:,i)) 
       if((trim(str%ion_names(i)).ne.'Z')) then
          if(sum(abs(str%ion_spins)).gt.0.0_gp) then
             write (unit_out_cell,'(a4,3f20.15,a6,f10.5)') trim(str%ion_names(i)),v(:)-floor(v(:)),' SPIN=',str%ion_spins(i)
          else
             write (unit_out_cell,'(a4,3f20.15)') trim(str%ion_names(i)),v(:)-floor(v(:))
          end if
       end if
    end do
    write (unit_out_cell,'(a)') "%ENDBLOCK POSITIONS_FRAC"
    write (unit_out_cell,*)

    if(fix_cell) then
       write (unit_out_cell,'(a)') "FIX_ALL_CELL : true"
       write (unit_out_cell,*)
    end if

    ! ** Write the kpoint spacing, if required

    if(len_trim(kpspacing)>0) then
       write (unit_out_cell,'(a)') kpspacing
       write (unit_out_cell,*)
    end if

    ! ** Let CASTEP find the symmetry

    write (unit_out_cell,'(a)') 'SYMMETRY_GENERATE'
    write (unit_out_cell,'(a)') 'SNAP_TO_SYMMETRY'
    write (unit_out_cell,*)

    ! ** Write out the pseudopotential information

    if(num_pot>0) then
       write (unit_out_cell,'(a)') "%BLOCK SPECIES_POT"
       do i=1,num_pot
          write (unit_out_cell,'(a)') spec_pot(i)
       end do
       write (unit_out_cell,'(a)') "%ENDBLOCK SPECIES_POT"
       write (unit_out_cell,*)
    end if

    ! ** Write out the applied external pressure, convert to GPa

    write (unit_out_cell,'(a)') "%BLOCK EXTERNAL_PRESSURE"
    write (unit_out_cell,'(3f10.5)') str%external_pressure(1,1)*160.21766208_gp,str%external_pressure(1,2)*160.21766208_gp,&
         str%external_pressure(1,3)*160.21766208_gp
    write (unit_out_cell,'(2F10.5)') str%external_pressure(2,2)*160.21766208_gp,str%external_pressure(2,3)*160.21766208_gp
    write (unit_out_cell,'(1F10.5)') str%external_pressure(3,3)*160.21766208_gp
    write (unit_out_cell,'(a)') "%ENDBLOCK EXTERNAL_PRESSURE"
    write (unit_out_cell,*) 

    ! ** Write out the pseudopotential information

    if(num_hub>0) then
       write (unit_out_cell,'(a)') "%BLOCK HUBBARD_U"
       do i=1,num_hub
          write (unit_out_cell,'(a)') hubbard_u(i)
       end do
       write (unit_out_cell,'(a)') "%ENDBLOCK HUBBARD_U"
       write (unit_out_cell,*)
    end if

    flush (unit_out_cell)

    close(unit=unit_out_cell)

    return

998 stop 'There is a problem writing the cell information. Stopping.'

  end subroutine write_cell

  subroutine write_xyze(str)

    type(structure), intent(in) :: str

    integer :: ni
    character(len=10) :: ctemp

    if(quiet) return

    write (ctemp,'(i10)') str%num_ions
    write (unit_xyze,'(a)') trim(adjustl(ctemp))

    write (unit_xyze,'(a,9f20.6,a)') 'Lattice="',reshape(str%lattice_car,(/9/)),'" Properties=species:S:1:pos:R:3'
    do ni=1,str%num_ions
       write (unit_xyze,'(a,3f20.13)') trim(adjustl(str%ion_names(ni))),str%ion_positions(:,ni)
    end do

    flush(unit_xyze)

  end subroutine write_xyze

  subroutine update_cell(str)

    type(structure), intent(inout) :: str

    str%volume = cell_volume(str%lattice_car)
    
    str%lattice_abc=cell_abc(str%lattice_car)
    
    str%lattice_rec=cell_rec(str%lattice_car)

    ! ** Ensure volume is positive
    
    str%volume=abs(str%volume)

  end subroutine update_cell

  function cell_volume(l)

    real(kind=gp), dimension(3,3), intent(in) :: l

    real(kind=gp) :: cell_volume

    cell_volume = l(1,1)*(l(2,2)*l(3,3)-l(3,2)*l(2,3))+l(2,1)*(l(3,2)*l(1,3)-l(1,2)*l(3,3))+l(3,1)*(l(1,2)*l(2,3)-l(2,2)*l(1,3))
    
  end function cell_volume

  function cell_car(labc)

    real(kind=gp), dimension(6), intent(in) :: labc

    real(kind=gp) :: cell_car(3,3)

    cell_car(:,1) = (/labc(1),0.0_gp,0.0_gp/)
    cell_car(:,2) = (/labc(2)*cos(dgrd*labc(6)),labc(2)*sin(dgrd*labc(6)),0.0_gp/)
    cell_car(1,3) = labc(3)*cos(dgrd*labc(5))
    cell_car(2,3) = labc(3)*(cos(dgrd*labc(4))-cos(dgrd*labc(5))*cos(dgrd*labc(6)))/sin(dgrd*labc(6))
    cell_car(3,3) = sqrt(labc(3)**2-cell_car(1,3)**2-cell_car(2,3)**2)

  end function cell_car
  
  function cell_abc(l)

    real(kind=gp), dimension(3,3), intent(in) :: l

    real(kind=gp) :: cell_abc(6)

    cell_abc(1) = sqrt(l(1,1)**2+l(2,1)**2+l(3,1)**2)
    cell_abc(2) = sqrt(l(1,2)**2+l(2,2)**2+l(3,2)**2)
    cell_abc(3) = sqrt(l(1,3)**2+l(2,3)**2+l(3,3)**2)
    cell_abc(4) = acos(dot_product(l(:,2),l(:,3))/cell_abc(2)/cell_abc(3))/dgrd
    cell_abc(5) = acos(dot_product(l(:,1),l(:,3))/cell_abc(1)/cell_abc(3))/dgrd
    cell_abc(6) = acos(dot_product(l(:,1),l(:,2))/cell_abc(1)/cell_abc(2))/dgrd
    
  end function cell_abc

  function cell_rec(l)

    real(kind=gp), dimension(3,3), intent(in) :: l

    real(kind=gp) :: cell_rec(3,3)

    real(kind=gp) :: vol

    cell_rec(1,1)=l(2,2)*l(3,3)-l(3,2)*l(2,3)
    cell_rec(2,1)=l(2,3)*l(3,1)-l(3,3)*l(2,1)
    cell_rec(3,1)=l(2,1)*l(3,2)-l(3,1)*l(2,2)
    cell_rec(1,2)=l(3,2)*l(1,3)-l(1,2)*l(3,3)
    cell_rec(2,2)=l(3,3)*l(1,1)-l(1,3)*l(3,1)
    cell_rec(3,2)=l(3,1)*l(1,2)-l(1,1)*l(3,2)
    cell_rec(1,3)=l(1,2)*l(2,3)-l(2,2)*l(1,3)
    cell_rec(2,3)=l(1,3)*l(2,1)-l(2,3)*l(1,1)
    cell_rec(3,3)=l(1,1)*l(2,2)-l(2,1)*l(1,2)

    vol=l(1,1)*(l(2,2)*l(3,3)-l(3,2)*l(2,3))+l(2,1)*(l(3,2)*l(1,3)-l(1,2)*l(3,3))+l(3,1)*(l(1,2)*l(2,3)-l(2,2)*l(1,3))
    
    cell_rec(:,:)=cell_rec(:,:)/vol
    
  end function cell_rec
  
!   subroutine random_nearly_cubic_supercell(str,ncell,mat,cmin)
!
!     type(structure),         intent(in)  :: str
!     integer      ,           intent(in)  :: ncell
!     integer, dimension(3,3), intent(out) :: mat
!     real(kind=gp),           intent(out) :: cmin
!
!     real(kind=gp), allocatable :: cmins(:)
!     real(kind=gp) :: labc(6),cost
!     integer :: n(1),ncount,nt,nthread
!     integer, allocatable :: minmat(:,:,:)
!
!     !$omp parallel sections
!     nthread=omp_get_num_threads()
!     !$omp end parallel sections
!     allocate(minmat(3,3,nthread),cmins(nthread))
!
!     cmins=huge(1.0_gp)
!     !$omp parallel do private(mat,labc,cost) schedule(dynamic)
!     do ncount=1,scell_attempts
!        nt=omp_get_thread_num()+1
!        call random_supercell(ncell,mat)
!        labc=cell_abc(matmul(str%lattice_car,mat))
!        cost=sqrt(sum((labc(1:3)-sum(labc(1:3)/3))**2))/sum(labc(1:3)/3)+sqrt(sum(cos(dgrd*labc(4:6))**2))
!        if(cost.lt.cmins(nt)) then
!           cmins(nt)=cost
!           minmat(:,:,nt)=mat
!        end if
!
!     end do
!     !$omp end parallel do
!
!     n=minloc(cmins)
!     cmin=cmins(n(1))
!     mat=minmat(:,:,n(1))
!
!   end subroutine random_nearly_cubic_supercell
!
!   subroutine enumerate_nearly_cubic_supercell(str,ncell,mat,cmin)
!
!     type(structure),         intent(in)  :: str
!     integer      ,           intent(in)  :: ncell
!     integer, dimension(3,3), intent(out) :: mat
!     real(kind=gp),           intent(out) :: cmin
!
!     real(kind=gp), allocatable :: cmins(:)
!
!     real(kind=gp) :: labc(6),cost
!     real(kind=dp) :: sc(3,3)
!
!     integer, allocatable :: minmat(:,:,:)
!
!     integer :: n(1),nmax=300,a,b,c,d,e,f,info,nthread,nt
!
!     !$omp parallel sections
!     nthread=omp_get_num_threads()
!     !$omp end parallel sections
!     allocate(minmat(3,3,nthread),cmins(nthread))
!
!     cmins=huge(1.0_gp)
!     !$omp parallel do private(mat,sc,info,labc,nt,cost,b,c,d,e,f) schedule(dynamic)
!     do a=1,min(ncell,nmax)
!        nt=omp_get_thread_num()+1
!        do b=1,min(ncell/a,nmax)
!           do c=1,min(ncell/a/b,nmax)
!              if(a*b*c.ne.ncell) cycle
!              do d=0,b-1
!                 do e=0,c-1
!                    do f=0,c-1
!
!                       mat=0 ; mat(1,1)=a ; mat(2,2)=b ; mat(3,3)=c ; mat(1,2)=d ; mat(1,3)=e ; mat(2,3)=f
!
!                       sc=mat
!
!                       info=spg_niggli_reduce(sc,1e-4_dp)
!
!                       mat=nint(sc)
!
!                       labc=cell_abc(matmul(str%lattice_car,mat))
!
!                       cost=sqrt(sum((labc(1:3)-sum(labc(1:3)/3))**2))/sum(labc(1:3)/3)+sqrt(sum(cos(dgrd*labc(4:6))**2))
!
!                       if(cost.lt.cmins(nt)) then
!                          cmins(nt)=cost
!                          minmat(:,:,nt)=mat
!                       end if
!
!                    end do
!                 end do
!              end do
!           end do
!        end do
!     end do
!     !$omp end parallel do
!
!     n=minloc(cmins)
!     cmin=cmins(n(1))
!     mat=minmat(:,:,n(1))
!
!   end subroutine enumerate_nearly_cubic_supercell

  subroutine shake_nearly_cubic_supercell(str,ncell,mat,cmin)

    type(structure),         intent(in)  :: str
    integer      ,           intent(in)  :: ncell
    integer, dimension(3,3), intent(out) :: mat
    real(kind=gp),           intent(out) :: cmin

    real(kind=gp), allocatable :: cmins(:)

    real(kind=gp) :: labc(6),cost

    integer, allocatable :: minmat(:,:,:) 

    integer :: n(1),a,i,j,nthread,nt

    !$omp parallel sections
    nthread=omp_get_num_threads()
    !$omp end parallel sections
    allocate(minmat(3,3,nthread),cmins(nthread))

    minmat=0

    cmins=huge(1.0_gp)
    !$omp parallel do private(mat,labc,nt,cost) schedule(dynamic)
    do a=1,scell_attempts/nthread
       nt=omp_get_thread_num()+1

       mat=0
       do while(cell_volume(matmul(str%lattice_car,mat)).le.1e-4_gp)
          do i=1,3
             do j=1,3
                mat(i,j)=minmat(i,j,nt)+(random_integer(3)-2)
             end do
          end do
       end do

       labc=cell_abc(matmul(str%lattice_car,mat))

       cost=sqrt(sum((labc(1:3)-sum(labc(1:3)/3))**2))/sum(labc(1:3)/3)+sqrt(sum(cos(dgrd*labc(4:6))**2))+&
            real(abs(ncell-num_super_cell(mat)),gp)/ncell

       if(cost.lt.cmins(nt)) then
          cmins(nt)=cost
          minmat(:,:,nt)=mat
       end if

    end do
    !$omp end parallel do

    n=minloc(cmins)
    mat=minmat(:,:,n(1))
    labc=cell_abc(matmul(str%lattice_car,mat))
    cmin=sqrt(sum((labc(1:3)-sum(labc(1:3)/3))**2))/sum(labc(1:3)/3)+sqrt(sum(cos(dgrd*labc(4:6))**2))

  end subroutine shake_nearly_cubic_supercell
  
  subroutine random_supercell(ncell,mat)

    integer,                 intent(in)  :: ncell
    integer, dimension(3,3), intent(out) :: mat

    real(kind=dp) :: sc(3,3)

    integer       :: info

    info=0
    
    do while(info.ne.1)

       mat=0
       
       do while(mat(1,1)*mat(2,2)*mat(3,3).ne.ncell)
          
          mat(1,1)=random_integer(ncell)
          
          mat(2,2)=random_integer(ncell/mat(1,1)+1)
          
          mat(3,3)=ncell/mat(1,1)/mat(2,2)
          
       end do
       
       mat(1,2)=random_integer(mat(2,2))-1
       mat(1,3)=random_integer(mat(3,3))-1
       mat(2,3)=random_integer(mat(3,3))-1 
       
       sc=mat

       info=spg_niggli_reduce(sc,1e-4_dp)

       mat=nint(sc)

    end do

  end subroutine random_supercell

  integer function num_super_cell(mat)

    implicit none

    integer, dimension(3,3), intent(inout) :: mat

    integer :: det,i,j,k

    det=0
    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       det=det+mat(1,i)*(mat(2,j)*mat(3,k)-mat(2,k)*mat(3,j))
    end do

    num_super_cell = det

  end function num_super_cell

  function supercell(mat,str)

    integer, dimension(3,3), intent(in) :: mat
    type(structure),         intent(in) :: str

    type(structure) :: supercell

    real(kind=gp) :: vfrac(3),mtemp(3,3)

    integer :: i,j,k,nc,ni,nii,nmax

    supercell=str

    ! ** Make new cell vectors

    supercell%lattice_car = matmul(str%lattice_car,mat)
    
    ! ** Update lattice info

    call update_cell(supercell)

    supercell%num_cells=nint(supercell%volume/str%volume)

    supercell%num_ions=str%num_ions*supercell%num_cells

    ! ** Reallocate

    if(allocated(supercell%subcell_vec)) deallocate(supercell%subcell_vec)
    allocate(supercell%subcell_vec(3,supercell%num_cells))

    if(allocated(supercell%ion_positions)) deallocate(supercell%ion_positions)
    allocate(supercell%ion_positions(3,supercell%num_ions))

    if(allocated(supercell%external_forces)) deallocate(supercell%external_forces)
    allocate(supercell%external_forces(3,supercell%num_ions))

    if(allocated(supercell%ion_names)) deallocate(supercell%ion_names)
    allocate(supercell%ion_names(supercell%num_ions))

    if(allocated(supercell%ion_names_spec)) deallocate(supercell%ion_names_spec)
    allocate(supercell%ion_names_spec(supercell%num_ions))

    if(allocated(supercell%ion_species)) deallocate(supercell%ion_species)
    allocate(supercell%ion_species(supercell%num_ions))

    if(allocated(supercell%ion_rad)) deallocate(supercell%ion_rad)
    allocate(supercell%ion_rad(supercell%num_ions))

    if(allocated(supercell%ion_vdw)) deallocate(supercell%ion_vdw)
    allocate(supercell%ion_vdw(supercell%num_ions))

    if(allocated(supercell%ion_mass)) deallocate(supercell%ion_mass)
    allocate(supercell%ion_mass(supercell%num_ions))

    if(allocated(supercell%ion_spins)) deallocate(supercell%ion_spins)
    allocate(supercell%ion_spins(supercell%num_ions))

    if(allocated(supercell%ion_cons)) deallocate(supercell%ion_cons)
    allocate(supercell%ion_cons(supercell%num_ions))

    ! ** Identify sub cells

    nmax=1
    nc=0

    supercell%subcell_vec(:,1)=0.0_gp

    do while(nc.lt.supercell%num_cells)

       mtemp = matmul(supercell%lattice_rec,str%lattice_car)

       nc=1
       do i=-nmax,nmax 
          do j=-nmax,nmax
             do k=-nmax,nmax

                if((i==0).and.(j==0).and.(k==0)) cycle

                vfrac(:)=i*mtemp(:,1)+j*mtemp(:,2)+k*mtemp(:,3)   

                if(all(vfrac.gt.-1e-4_gp).and.all(vfrac.lt.1.0_gp-1e-4_gp)) then
                   nc=nc+1
                   if(nc.gt.supercell%num_cells) stop 'supercell: subvector generation failed - nc > num_cells'
                   supercell%subcell_vec(:,nc) = i*str%lattice_car(:,1)+j*str%lattice_car(:,2)+k*str%lattice_car(:,3)
                end if

             end do
          end do
       end do

       nmax=nmax*2

    end do

    ! ** Fill

    nii=0
    do nc=1,supercell%num_cells
       do ni=1,str%num_ions

          nii=nii+1

          supercell%ion_positions(:,nii) = str%ion_positions(:,ni)+supercell%subcell_vec(:,nc)
          supercell%external_forces(:,nii) = str%external_forces(:,ni)

          supercell%ion_names(nii) = str%ion_names(ni)
          supercell%ion_names_spec(nii) = str%ion_names_spec(ni)
          supercell%ion_species(nii) = str%ion_species(ni)
          supercell%ion_rad(nii) = str%ion_rad(ni)
          supercell%ion_vdw(nii) = str%ion_vdw(ni)
          supercell%ion_mass(nii) = str%ion_mass(ni)
          supercell%ion_spins(nii) = str%ion_spins(ni)
          supercell%ion_cons(nii) = str%ion_cons(ni)

       end do
    end do

    if(symmgen.and.(.not.any(supercell%ion_cons))) call findsymm_cell(supercell)

    ! ** Identify symmetry equivalent ions

    if(supercell%num_symm.gt.1) then

       if(allocated(supercell%ion_equiv)) deallocate(supercell%ion_equiv)

       allocate(supercell%ion_equiv(supercell%num_symm,supercell%num_ions))

       nmax=1
       supercell%ion_equiv=-1

       do while(any(supercell%ion_equiv.lt.0))

          supercell%ion_equiv=-1

          do n=1,supercell%num_ions
             do m=1,supercell%num_ions

                do ns=1,supercell%num_symm

                   do na=-nmax,nmax
                      do nb=-nmax,nmax
                         do nc=-nmax,nmax

                            v(:) = supercell%ion_positions(:,n) - &
                                 matmul(supercell%symm_ops(1:3,1:3,ns),supercell%ion_positions(:,m)) - &
                                 supercell%symm_ops(1,4,ns)*supercell%lattice_car(1:3,1) - &
                                 supercell%symm_ops(2,4,ns)*supercell%lattice_car(1:3,2) - &
                                 supercell%symm_ops(3,4,ns)*supercell%lattice_car(1:3,3) - &
                                 na*supercell%lattice_car(:,1)-nb*supercell%lattice_car(:,2)-nc*supercell%lattice_car(:,3)

                            if(dot_product(v,v).lt.1e-4_gp) supercell%ion_equiv(ns,n)=m

                         end do
                      end do
                   end do


                end do


             end do
          end do
          nmax=nmax+1
       end do

    end if

  end function supercell

  function car2habc(car)

    real(kind=gp), dimension(3,3), intent(in) :: car

    real(kind=gp), dimension(3) :: car2habc

    real(kind=gp), dimension(3) :: acb,acc,ccb
    
    real(kind=gp) :: vol

    vol=cell_volume(car)
  
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
  
  subroutine findsymm_cell(str)

    type(structure), intent(inout) :: str

    real(kind=gp) :: lattice_trans(3,3),ion_fractional(3,str%num_ions)

    type(SpglibDataset) :: dset

    ! ** spglib uses the transpose of str%lattice_car

    lattice_trans=transpose(str%lattice_car)

    ! ** fractional positions

    do i=1,str%num_ions
       ion_fractional(:,i)=matmul(str%lattice_rec,str%ion_positions(:,i)) 
    end do

    dset=spg_get_dataset(real(lattice_trans,dp),real(ion_fractional,dp),str%ion_species,str%num_ions,real(abs(symmtol),dp))

    str%num_symm=dset%n_operations

    ! ** Get the international symbol for the space group based on the provided tolerance

    if(dset%spacegroup_number/=0) then
       str%spacegroup=deunder(trim(dset%international_symbol))
    else
       str%spacegroup='P1'
    end if

    if(.not.quiet) then
       write (stdout,'(a,a)')  ' space group    : ',trim(str%spacegroup)
       write (stdout,'(a,i0)') ' num symmops    : ',str%num_symm
       write (stdout,'(a,i0)') ' num ions       : ',str%num_ions
    end if

    if(allocated(str%symm_ops)) deallocate(str%symm_ops)
    allocate(str%symm_ops(3,4,str%num_symm))

    do ns=1,str%num_symm
       str%symm_ops(1:3,1:3,ns)=matmul(str%lattice_car,matmul(transpose(dset%rotations(1:3,1:3,ns)),str%lattice_rec)) 
       str%symm_ops(1:3,4,ns)=dset%translations(1:3,ns) 
    end do

  end subroutine findsymm_cell

  function deunder(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=len(string))   :: deunder

    integer                      :: i

    deunder=string
    do while (scan(deunder,'_').gt.0)
       i=scan(deunder,'_')
       deunder=trim(deunder(1:i-1))//trim(deunder(i+1:))
    end do

  end function deunder

end module cell
