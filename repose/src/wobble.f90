!==================================================================================!
!                                      wobble                                      !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This is the main program of wobble - data derived potential lattice dynamics     !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program wobble

  use constants
  use omp_lib
  use rng
  use ld
  use cell
  use ddp

  implicit none

  real(kind=gp) :: total_time,start_time,end_time

  integer :: nomp
  
  type(structure) :: struct

  call cpu_time(start_time)

  ! ** Get random seed

  call init_pseudorandom()

  ! ** Get the command line arguments

  call get_arguments(seedname)

  ! ** Read the cell data

  call read_cell(struct)

  call update_cell(struct)

  ! ** Read the pair potential

  if(.not.dryrun) call read_ddp()

  ! ** Lattice dynamics
  
  call ld_calc(struct)

  call cpu_time(end_time)

  total_time=end_time-start_time

  if(.not.quiet) then

     write (stderr,*)
     write (stderr,'(a)') '------------------------------------'
     write (stderr,'(a,f13.4,a,f7.2,a)') 'total   : ',total_time,' sec',100.0_gp,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') 'ddp     : ',ball_time+feat_time+netw_time,' sec',&
          100*(ball_time+feat_time+netw_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -ball  : ',ball_time,' sec',100*(ball_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -feat  : ',feat_time,' sec',100*(feat_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -netw  : ',netw_time,' sec',100*(netw_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') 'ld      : ',ld_time,' sec',100*(ld_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -cell  : ',cell_time,' sec',100*(cell_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -dmat  : ',dm_time,' sec',100*(dm_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -therm : ',therm_time,' sec',100*(therm_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -map   : ',map_time,' sec',100*(map_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -dsyev : ',dsyev_time,' sec',100*(dsyev_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -kpath : ',kpath_time,' sec',100*(kpath_time)/total_time,' %'
     write (stderr,'(a,f13.4,a,f7.2,a)') ' -dos   : ',dos_time,' sec',100*(dos_time)/total_time,' %'

  end if

contains

  ! *******************************************************************************************

  subroutine banner()

    if(quiet) return

    write (stderr,'(a)') "                                              ,--,              "
    write (stderr,'(a)') "                          ,---,     ,---,   ,--.'|              "
    write (stderr,'(a)') "         .---.   ,---.  ,---.'|   ,---.'|   |  | :              "
    write (stderr,'(a)') "        /. ./|  '   ,'\ |   | :   |   | :   :  : '              "
    write (stderr,'(a)') "     .-'-. ' | /   /   |:   : :   :   : :   |  ' |      ,---.   "
    write (stderr,'(a)') "    /___/ \: |.   ; ,. ::     |,-.:     |,-.'  | |     /     \  "
    write (stderr,'(a)') " .-'.. '   ' .'   | |: :|   : '  ||   : '  ||  | :    /    /  | "
    write (stderr,'(a)') "/___/ \:     ''   | .; :|   |  / :|   |  / :'  : |__ .    ' / | "
    write (stderr,'(a)') ".   \  ' .\   |   :    |'   : |: |'   : |: ||  | '.'|'   ;   /| "
    write (stderr,'(a)') " \   \   ' \ | \   \  / |   | '/ :|   | '/ :;  :    ;'   |  / | "
    write (stderr,'(a)') "  \   \  |--'   `----'  |   :    ||   :    ||  ,   / |   :    | "
    write (stderr,'(a)') "   \   \ |              /    \  / /    \  /  ---`-'   \   \  /  "
    write (stderr,'(a)') "    '---'               `-'----'  `-'----'             `----'   "
    write (stderr,'(a)') "                                                                "
    write (stderr,'(a)') "                                                                "
    write (stderr,'(a)') "        Author: Chris J. Pickard, Cambridge, 2020-23 (c)        "
    write (stderr,'(a)') "         lattice dynamics using data derived potentials         "
    write (stderr,'(a)') "                                                                "
    write (stderr,'(a,*((i0,1x)))')  ' num threads    : ', omp_get_max_threads()
    write (stderr,'(a,*((i0,1x)))')  ' max threads    : ', omp_get_num_procs()
    write (stderr,*)

  end subroutine banner

  subroutine get_arguments(name)

    character(len=80), intent(out) :: name

    character(len=80), allocatable,dimension(:) :: cargs

    integer :: i,iargc,num_args

    num_args = iargc()

    if(num_args.eq.0) goto 99

    allocate(cargs(num_args))

    do i=1,num_args
       call getarg(i,cargs(i))
    end do

    name = trim(cargs(num_args))

    if(any(cargs.eq."-h")) goto 99

    nomp=omp_get_num_procs()
    if(any(cargs.eq."-ompnp")) then           ! * Number of omp threads (note - not mpi)
       do i=1,num_args-1
          if(cargs(i).eq."-ompnp") exit
       end do
       read(cargs(i+1),*) nomp
       nomp=min(nomp,omp_get_num_procs())
    end if
    call omp_set_num_threads(nomp) 

    if(any(cargs.eq."-q")) quiet=.true.       ! * Quiet, minimal output

    ncells=-1
    if(any(cargs.eq."-ncell")) then            ! * Number of cells in supercell
       do i=1,num_args-1
          if(cargs(i).eq."-ncell") exit
       end do
       read(cargs(i+1),*) ncells
       ncells=max(ncells,1)
    end if

    natoms=1024
    if(any(cargs.eq."-natom")) then            ! * Number of atoms in supercell
       do i=1,num_args-1
          if(cargs(i).eq."-natom") exit
       end do
       read(cargs(i+1),*) natoms
    end if

    if(any(cargs.eq."-scmax")) then            ! * Maximum number of attempts to build supercell
       do i=1,num_args-1
          if(cargs(i).eq."-scmax") exit
       end do
       read(cargs(i+1),*) scell_attempts
       scell_attempts=max(scell_attempts,1)
    end if

    ntemp=11
    if(any(cargs.eq."-ntemp")) then            ! * Number of temperatures to compute free energy
       do i=1,num_args-1
          if(cargs(i).eq."-ntemp") exit
       end do
       read(cargs(i+1),*) ntemp
       ntemp=max(ntemp,1)
    end if

    nsamp=10000
    if(any(cargs.eq."-nsamp")) then            ! * Number of samples for density of states
       do i=1,num_args-1
          if(cargs(i).eq."-nsamp") exit
       end do
       read(cargs(i+1),*) nsamp
       nsamp=max(nsamp,1)
    end if

    nkpts=101
    if(any(cargs.eq."-nkpts")) then            ! * Number of kpoints for dispersion plot
       do i=1,num_args-1
          if(cargs(i).eq."-nkpts") exit
       end do
       read(cargs(i+1),*) nkpts
       nkpts=max(nkpts,1)
    end if

    tmax=500.0_gp
    if(any(cargs.eq."-tmax")) then            ! * Maximum temperature
       do i=1,num_args-1
          if(cargs(i).eq."-tmax") exit
       end do
       read(cargs(i+1),*) tmax
    end if

    rmin=0.25_gp
    if(any(cargs.eq."-rmin")) then            ! * Minimum tolerated contact
       do i=1,num_args-1
          if(cargs(i).eq."-rmin") exit
       end do
       read(cargs(i+1),*) rmin
    end if
    
    acore=0.0_gp
    rcore=-1.0_gp

    unit='meV'
    if(any(cargs.eq."-unit")) then            ! * Energy units (meV,cm-1,THz)
       do i=1,num_args-1
          if(cargs(i).eq."-unit") exit
       end do
       read(cargs(i+1),*) unit
    end if

    if(any(cargs.eq."-c")) cluster=.true.     ! * No periodic boundary conditions
    
    if(any(cargs.eq."-therm")) thermo=.true.     ! * Compute thermodynamics
    
    if(any(cargs.eq."-disp")) dispersion=.true.     ! * Compute dispersion
    
    if(any(cargs.eq."-dos")) dos=.true.     ! * Compute density of states
    
    if(any(cargs.eq."-dryrun")) dryrun=.true.     ! * Just build supercell

    deallocate(cargs)

    call banner()

    return

99  write (stdout,*) 'Usage: wobble [-ompnp] [-q] [-c] [-thermo] [-disp] [-dos] [-ncell] [-natom] '
    write (stdout,*) '        [-scmax] [-ntemp] [-nsamp] [-nkpts] [-tmax] -rmin] [-unit] [-dryrun] [-h] <seedname>'
    write (stdout,*) 
    write (stdout,*) 'Example: `wobble C` will compute the lattice dynamics of the structure in C.cell'
    write (stdout,*)
    write (stdout,*) '  -ompnp  : Number of threads'
    write (stdout,*) '  -q      : Minimal output'
    write (stdout,*) '  -c      : Cluster'
    write (stdout,*) '  -therm  : Compute thermodynamics'
    write (stdout,*) '  -disp   : Compute dispersion'
    write (stdout,*) '  -dos    : Compute density of states'
    write (stdout,*) '  -ncell  : Number of cells in supercell'
    write (stdout,*) '  -natom  : Number of atoms in supercell'
    write (stdout,*) '  -scmax  : Maximum attempts to construct supercell'
    write (stdout,*) '  -ntemp  : Number of temperatures to compute free energy'
    write (stdout,*) '  -nsamp  : Number of samples for density of states '
    write (stdout,*) '  -nkpts  : Number of kpoints for dispersion plot'
    write (stdout,*) '  -tmax   : Maximum temperature'
    write (stdout,*) '  -rmin   : Minimum tolerated contact'
    write (stdout,*) '  -unit   : Energy units (meV,cm-1,THz)'
    write (stdout,*) '  -dryrun : Set up calculation only'
    write (stdout,*) '  -h      : Display this message'
    write (stdout,*)

    stop

  end subroutine get_arguments

  ! *******************************************************************************************

end program wobble
