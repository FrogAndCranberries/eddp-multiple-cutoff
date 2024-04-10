!==================================================================================!
!                                      repose                                      !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This is the main program of repose - a data derived potential structure optimiser!
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program repose_new

  use constants
  use omp_lib
  use rng
  use cell
  use ddp_new
  use opt

  implicit none

  real(kind=gp) :: etol,total_time,start_time,end_time

  real(kind=gp) :: dmin,dmax
  
  integer :: num_steps,maxsteps,nomp

  type(structure) :: struct

  call cpu_time(start_time)
  
  ! ** Get random seed

  call init_pseudorandom()

  ! ** Get the command line arguments

  call get_arguments(seedname)

  ! ** Read the cell data

  call read_cell(struct)

  call update_cell(struct)

  !if (optimise.and.(.not.fix_cell)) call compact_cell() !! SORT THIS OUT WITH SPGLIB

  ! ** Read the pair potential

  call read_ddp()

  ! ** Test gradients

  if(test) then
     call opt_test(struct) ; stop  
  end if
  
  ! ** Finish relaxation
  
  num_steps=maxsteps
  call opt_tpsd(struct,num_steps,etol,cell_fix=fix_cell)
  
  !if (optimise.and.(.not.fix_cell)) call compact_cell() !! SORT THIS OUT WITH SPGLIB

  if(sqrt(struct%variance)/struct%num_ions.lt.dmin) stop 'repose: deviation too small'
  if(sqrt(struct%variance)/struct%num_ions.gt.dmax) stop 'repose: deviation too large'
  
  if(.not.converged) then
     write(stdout,'(a,i6,a)') "Did not converge after ",total_steps," steps"
  else
     write(stdout,'(a,i6,a)') "Converged in ",total_steps," steps"
  end if
  
  write (stdout,*)
  write (stdout,'(a11,f25.10)') 'Volume:    ',  struct%volume
  write (stdout,'(a11,f25.10)') 'Pressure:  ',  (struct%stress(1,1)+struct%stress(2,2)+struct%stress(3,3))/3.0_gp*160.21766208_gp
  write (stdout,'(a11,f25.10)') 'Energy:    ',  struct%energy
  write (stdout,'(a11,f25.10)') 'Enthalpy:  ',  struct%enthalpy
  write (stdout,'(a11,f25.10)') 'Deviation:  ',  sqrt(struct%variance)
  flush (stdout)

  call cpu_time(end_time)

  total_time=end_time-start_time
  
  write (stderr,*)
  write (stderr,'(a)') '------------------------------------'
  write (stderr,'(a,f13.4,a,f7.2,a)') 'total  : ',total_time,' sec',100.0_gp,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') 'ddp    : ',ball_time+feat_time+netw_time,' sec',&
       100*(ball_time+feat_time+netw_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -ball : ',ball_time,' sec',100*(ball_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -feat : ',feat_time,' sec',100*(feat_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -netw : ',netw_time,' sec',100*(netw_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') 'opt    : ',prec_time,' sec',100*(prec_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -prec : ',prec_time,' sec',100*(prec_time)/total_time,' %'
  
  ! ** Write the final structure

  call write_cell(struct)

contains

  ! *******************************************************************************************

  subroutine banner()

    if(quiet) return

    write (stdout,'(a)') "                                                                   "
    write (stdout,'(a)') " ***  ****                 ****       ****       ****              "
    write (stdout,'(a)') "  **** **** *    ***      * ***  *   * ***  *   * **** *    ***    "
    write (stdout,'(a)') "   **   ****    * ***    *   ****   *   ****   **  ****    * ***   "
    write (stdout,'(a)') "   **          *   ***  **    **   **    **   ****        *   ***  "
    write (stdout,'(a)') "   **         **    *** **    **   **    **     ***      **    *** "
    write (stdout,'(a)') "   **         ********  **    **   **    **       ***    ********  "
    write (stdout,'(a)') "   **         *******   **    **   **    **         ***  *******   "
    write (stdout,'(a)') "   **         **        **    **   **    **    ****  **  **        "
    write (stdout,'(a)') "   ***        ****    * *******     ******    * **** *   ****    * "
    write (stdout,'(a)') "    ***        *******  ******       ****        ****     *******  "
    write (stdout,'(a)') "                *****   **                                 *****   "
    write (stdout,'(a)') "                        **                                         "
    write (stdout,'(a)') "                        **                                         "
    write (stdout,'(a)') "                         **                                        "                                                              
    write (stdout,'(a)') "                                                                   "
    write (stdout,'(a)') "         Author: Chris J. Pickard, Cambridge, 2020-23 (c)          "
    write (stdout,'(a)') "    relaxing positions and energies using data derived potentials  "
    write (stdout,'(a)') "                                                                   "
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
    if(any(cargs.eq."-ompnp")) then               ! * Number of omp threads (note - not mpi)
       do i=1,num_args-1
          if(cargs(i).eq."-ompnp") exit
       end do
       read(cargs(i+1),*) nomp
       nomp=min(nomp,omp_get_num_procs())
    end if
    call omp_set_num_threads(nomp) 
    
    if(any(cargs.eq."-test")) test=.true.     ! * Test derivatives
    if(any(cargs.eq."-t")) track=.true.       ! * Track the optimisation
    if(any(cargs.eq."-n")) optimise=.false.   ! * Do not relax the structure
    if(any(cargs.eq."-q")) quiet=.true.       ! * Quiet, minimal output
    maxsteps=100000
    if(any(cargs.eq."-m")) then               ! * Maximum number of optimisation steps
       do i=1,num_args-1
          if(cargs(i).eq."-m") exit
       end do
       read(cargs(i+1),*) maxsteps
    end if
    track_every=1
    if(any(cargs.eq."-te")) then              ! * Track every N steps
       do i=1,num_args-1
          if(cargs(i).eq."-te") exit
       end do
       read(cargs(i+1),*) track_every
    end if
    if(any(cargs.eq."-v")) then               ! * Volume override
       do i=1,num_args-1
          if(cargs(i).eq."-v") exit
       end do
       read(cargs(i+1),*) voveride
    end if
    if(any(cargs.eq."-p")) then               ! * Pressure override
       do i=1,num_args-1
          if(cargs(i).eq."-p") exit
       end do
       read(cargs(i+1),*) poveride
    end if
    poveride=poveride/160.21766208_gp ! Convert to eV/ang^3

    etol=1e-6_gp
    if(any(cargs.eq."-tol")) then             ! * Convergence tolerance
       do i=1,num_args-1
          if(cargs(i).eq."-tol") exit
       end do
       read(cargs(i+1),*) etol
    end if

    gamma=1.0_gp
    if(any(cargs.eq."-g")) then               ! * Gamma cell damping
       do i=1,num_args-1
          if(cargs(i).eq."-g") exit
       end do
       read(cargs(i+1),*) gamma
    end if

    rmin=0.25_gp
    if(any(cargs.eq."-rmin")) then            ! * Minimum tolerated contact
       do i=1,num_args-1
          if(cargs(i).eq."-rmin") exit
       end do
       read(cargs(i+1),*) rmin
    end if

    rcore=rmin*2
    if(any(cargs.eq."-r")) then               ! * Hard core radius
       do i=1,num_args-1
          if(cargs(i).eq."-r") exit
       end do
       read(cargs(i+1),*) rcore
    end if
    
    acore=1.0_gp
    if(any(cargs.eq."-a")) then               ! * Hard core strength
       do i=1,num_args-1
          if(cargs(i).eq."-a") exit
       end do
       read(cargs(i+1),*) acore
    end if

    dmin=-huge(1.0_gp)
    if(any(cargs.eq."-dmin")) then            ! * Minimum tolerated deviation
       do i=1,num_args-1
          if(cargs(i).eq."-dmin") exit
       end do
       read(cargs(i+1),*) dmin
    end if

    dmax=huge(1.0_gp)
    if(any(cargs.eq."-dmax")) then            ! * Maximum tolerated deviation
       do i=1,num_args-1
          if(cargs(i).eq."-dmax") exit
       end do
       read(cargs(i+1),*) dmax
    end if
    
    if(any(cargs.eq."-f")) fix_cell=.true.    ! * Fix the unit cell

    if(any(cargs.eq."-c")) cluster=.true.     ! * No periodic boundary conditions

    deallocate(cargs)

    call banner()

    return

99  write (stdout,*) 'Usage: repose [-ompnp] [-test] [-t] [-te] [-n] [-q] [-m] [-v] [-p] [-f] [-c] [-tol] [-r]'
    write (stdout,*) '             [-a] [-g] [-rmin] [-devmax] [-devmin] [-h] <seedname>'
    write (stdout,*) 
    write (stdout,*) 'Example: `repose -t C` will track a relaxation of the structure contained in C.cell'
    write (stdout,*)
    write (stdout,*) '  -ompnp : Number of threads'
    write (stdout,*) '  -test  : Test derivatives'
    write (stdout,*) '  -t     : Track relaxation'
    write (stdout,*) '  -te    : Track every N steps'
    write (stdout,*) '  -n     : No structural relaxation'
    write (stdout,*) '  -q     : Minimal output'
    write (stdout,*) '  -m     : Max relaxation steps'
    write (stdout,*) '  -v     : Volume overide'
    write (stdout,*) '  -p     : Pressure'
    write (stdout,*) '  -f     : Fix unit cell'
    write (stdout,*) '  -c     : Cluster'
    write (stdout,*) '  -tol   : Convergence thresold'
    write (stdout,*) '  -r     : Hard core radius'
    write (stdout,*) '  -a     : Hard core strength'
    write (stdout,*) '  -g     : Gamma cell damping'
    write (stdout,*) '  -rmin  : Minimum tolerated contact'
    write (stdout,*) '  -dmin  : Minimum tolerated deviation'
    write (stdout,*) '  -dmax  : Maximum tolerated deviation'
    write (stdout,*) '  -h     : Display this message'
    write (stdout,*)

    stop

  end subroutine get_arguments

  ! *******************************************************************************************

end program repose_new
