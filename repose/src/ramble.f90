!==================================================================================!
!                                      ramble                                      !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This is the main program of ramble - a data derived potential molecular dynamics !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program ramble

  use constants
  use omp_lib
  use rng
  use md
  use cell
  use ddp
  use opt

  implicit none

  real(kind=gp) :: etol,total_time,start_time,end_time
  real(kind=gp) :: cst
  
  integer :: num_steps,maxsteps,nomp,mat(3,3),i,j,k,ncells,natoms

  type(structure) :: struct,super
     
  logical :: optim=.false.,track_save

  call cpu_time(start_time)

  ! ** Get random seed

  call init_pseudorandom()

  ! ** Get the command line arguments

  call get_arguments(seedname)

  ! ** Read the cell data

  call read_cell(struct)

  call update_cell(struct)

  ! ** Read the pair potential

  call read_ddp()

  ! ** Construct supercell
  
  if(all(struct%supercell_matrix.eq.0)) then

     if(natoms.gt.0) then
        
        if(ncells.le.0) ncells=(natoms-1)/struct%num_ions+1
        
        call shake_nearly_cubic_supercell(struct,ncells,mat,cst)
        
        write (stderr,'(a,f8.4)') ' deviation from cubic: ',cst
        
     else

        if(ncells.gt.0) then
           
           call shake_nearly_cubic_supercell(struct,ncells,mat,cst)
        
           write (stderr,'(a,f8.4)') ' deviation from cubic: ',cst
           
        else

           mat=0
           mat(1,1)=1 ; mat(2,2)=1 ; mat(3,3)=1
        
        end if
        
     end if
  
  else

     mat=struct%supercell_matrix

  end if

  ncells=0
  do i=1,3
     j=mod(i,3)+1
     k=mod(j,3)+1
     ncells=ncells+mat(1,i)*(mat(2,j)*mat(3,k)-mat(2,k)*mat(3,j))
  end do

  super=supercell(mat,struct)
  call write_cell(super)

  write (stderr,'(a,i0,a,i0,a)') ' supercell matrix generating ',ncells,' cells containing ',super%num_ions,' atoms:'
  write (stderr,*)
  write (stderr,'(3i5)') mat
  write (stderr,*)
  write (stderr,'(1x,3f8.3)') super%lattice_abc(1:3)
  write (stderr,'(1x,3f8.2)') super%lattice_abc(4:6)
  write (stderr,*)

  if(optim) then
     track_save=track
     track=.false.
     num_steps=maxsteps
     call opt_tpsd(super,num_steps,tol=1e-6_gp,cell_fix=.true.)
     num_steps=maxsteps
     call opt_tpsd(super,num_steps,tol=1e-6_gp,cell_fix=fix_cell)
     track=track_save
     total_steps=0
  end if

  num_steps=maxsteps
  call md_cell(super,num_steps,fix_cell)

  if(optim) then
     track_save=track
     track=.false.
     num_steps=maxsteps
     call opt_tpsd(super,num_steps,tol=1e-6_gp,cell_fix=.true.)
     num_steps=maxsteps
     call opt_tpsd(super,num_steps,tol=1e-6_gp,cell_fix=fix_cell)
     track=track_save
  end if
  
  write (stdout,*)
  write (stdout,'(a11,f25.10)') 'Volume:    ',  super%volume
  write (stdout,'(a11,f25.10)') 'Pressure:  ',  (super%stress(1,1)+super%stress(2,2)+super%stress(3,3))/3.0_gp*160.21766208_gp
  write (stdout,'(a11,f25.10)') 'Energy:    ',  super%energy
  write (stdout,'(a11,f25.10)') 'Enthalpy:  ',  super%enthalpy
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
  write (stderr,'(a,f13.4,a,f7.2,a)') 'md     : ',md_time,' sec',100*(md_time)/total_time,' %'
  write (stderr,'(a,f13.4,a,f7.2,a)') ' -md : ',md_time,' sec',100*(md_time)/total_time,' %'

  ! ** Write the final structure

  call write_cell(super)

contains

  ! *******************************************************************************************

  subroutine banner()

    if(quiet) return

    write (stdout,'(a)') '                                                                   '
    write (stdout,'(a)') '                                     888      888                  '
    write (stdout,'(a)') '                                     888      888                  '
    write (stdout,'(a)') '                                     888      888                  '
    write (stdout,'(a)') '       888d888 8888b.  88888b.d88b.  88888b.  888  .d88b.          '
    write (stdout,'(a)') '       888P"      "88b 888 "888 "88b 888 "88b 888 d8P  Y8b         '
    write (stdout,'(a)') '       888    .d888888 888  888  888 888  888 888 88888888         '
    write (stdout,'(a)') '       888    888  888 888  888  888 888 d88P 888 Y8b.             '
    write (stdout,'(a)') '       888    "Y888888 888  888  888 88888P"  888  "Y8888          '                         
    write (stdout,'(a)') '                                                                   '
    write (stdout,'(a)') '        Author: Chris J. Pickard, Cambridge, 2020-23 (c)           '
    write (stdout,'(a)') '        molecular dynamics using data derived potentials           '
    write (stdout,'(a)') '                                                                   '
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

    if(any(cargs.eq."-t")) track=.true.       ! * Track the optimisation
    
    if(any(cargs.eq."-q")) quiet=.true.       ! * Quiet, minimal output
    maxsteps=huge(1)
    if(any(cargs.eq."-m")) then               ! * Maximum number of optimisation steps
       do i=1,num_args-1
          if(cargs(i).eq."-m") exit
       end do
       read(cargs(i+1),*) maxsteps
    end if
    num_run=1000
    if(any(cargs.eq."-n")) then               ! * Number of steps for running average
       do i=1,num_args-1
          if(cargs(i).eq."-n") exit
       end do
       read(cargs(i+1),*) num_run
    end if
    track_every=30
    if(any(cargs.eq."-te")) then              ! * Track every N steps
       do i=1,num_args-1
          if(cargs(i).eq."-te") exit
       end do
       read(cargs(i+1),*) track_every
    end if
    damping_rate=0.01_gp
    if(any(cargs.eq."-dr")) then              ! * Damping rate
       do i=1,num_args-1
          if(cargs(i).eq."-dr") exit
       end do
       read(cargs(i+1),*) damping_rate
    end if
    if(damping_rate.ge.1.0_gp) stop 'damping_rate should be less than one'
    
    damping_cell=0.9_gp
    if(any(cargs.eq."-dc")) then              ! * Damping cell
       do i=1,num_args-1
          if(cargs(i).eq."-dc") exit
       end do
       read(cargs(i+1),*) damping_cell
    end if
    if(damping_cell.ge.1.0_gp) stop 'damping_cell should be less than one'
    
    heat_rate=0.0_gp
    if(any(cargs.eq."-hr")) then              ! * Heat rate (K/ps)
       do i=1,num_args-1
          if(cargs(i).eq."-hr") exit
       end do
       read(cargs(i+1),*) heat_rate
    end if
    
    quench_rate=0.0_gp
    if(any(cargs.eq."-qr")) then              ! * Quench rate (K/ps)
       do i=1,num_args-1
          if(cargs(i).eq."-qr") exit
       end do
       read(cargs(i+1),*) quench_rate
    end if
    
    quench_time=-1.0_gp
    if(any(cargs.eq."-qt")) then              ! * Quench time (ps)
       do i=1,num_args-1
          if(cargs(i).eq."-qt") exit
       end do
       read(cargs(i+1),*) quench_time
    end if
    
    quench_end=0.0_gp
    if(any(cargs.eq."-qe")) then              ! * Quench end temperature (K)
       do i=1,num_args-1
          if(cargs(i).eq."-qe") exit
       end do
       read(cargs(i+1),*) quench_end
    end if
    
    stop_time=huge(1.0_gp)
    if(any(cargs.eq."-st")) then              ! * Stop time (ps)
       do i=1,num_args-1
          if(cargs(i).eq."-st") exit
       end do
       read(cargs(i+1),*) stop_time
    end if
    
    time_step=5e-2_gp
    if(any(cargs.eq."-ts")) then              ! * Time step
       do i=1,num_args-1
          if(cargs(i).eq."-ts") exit
       end do
       read(cargs(i+1),*) time_step
       time_step=time_step/ntu2ps/1000       ! convert from fs
    end if
    
    target_temp=300.0_gp
    if(any(cargs.eq."-tt")) then              ! * Target temperature
       do i=1,num_args-1
          if(cargs(i).eq."-tt") exit
       end do
       read(cargs(i+1),*) target_temp
    end if
    initial_temp=target_temp
    if(any(cargs.eq."-ti")) then              ! * Target temperature
       do i=1,num_args-1
          if(cargs(i).eq."-ti") exit
       end do
       read(cargs(i+1),*) initial_temp
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

    ncells=-1
    if(any(cargs.eq."-ncell")) then            ! * Number of cells in supercell
       do i=1,num_args-1
          if(cargs(i).eq."-ncell") exit
       end do
       read(cargs(i+1),*) ncells
       ncells=max(ncells,1)
    end if

    natoms=-1
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
    
    if(any(cargs.eq."-f")) fix_cell=.true.    ! * Fix the unit cell

    if(any(cargs.eq."-mc")) map_to_cell=.true.  ! * Map contents back to unit cell

    if(any(cargs.eq."-c")) cluster=.true.     ! * No periodic boundary conditions

    if(any(cargs.eq."-o")) optim=.true.    ! * Optimise before and after

    if(any(cargs.eq."-niggli")) niggli=.true.    ! * Niggli reduce during trajectory

    deallocate(cargs)

    call banner()

    return

99  write (stdout,*) 'Usage: ramble [-ompnp] [-t] [-te] [-dr] [-dc] [-hr] [-qr] [-qt] [-qe] [-st] [-ts] [-ti] [-tt] [-q] [-m]'
    write (stdout,*) '              [-n] [-v] [-p] [-f] [-mc] [-c] [-tol] [-r] [-a] [-rmin] [-h] [-ncell] [-natom] [-niggli]'
    write (stdout,*) '              [-scmax] [-o] <seedname>'
    write (stdout,*) 
    write (stdout,*) 'Example: `ramble -t C` will track the dynamics of the structure contained in C.cell'
    write (stdout,*)
    write (stdout,*) '  -ompnp  : Number of OMP threads'
    write (stdout,*) '  -t      : Track dynamics'
    write (stdout,*) '  -te     : Track every N steps'
    write (stdout,*) '  -dr     : Damping rate'
    write (stdout,*) '  -dc     : Damping cell'
    write (stdout,*) '  -ts     : Time step (fs)'
    write (stdout,*) '  -ti     : Initial temperature (K)'
    write (stdout,*) '  -tt     : Target temperature (K)'
    write (stdout,*) '  -hr     : Heat rate (K/ps)'
    write (stdout,*) '  -qr     : Quench rate (K/ps)'
    write (stdout,*) '  -qt     : Quench time (ps)'
    write (stdout,*) '  -qe     : Quench end temperature (K)'
    write (stdout,*) '  -st     : Stop time (ps)'
    write (stdout,*) '  -q      : Minimal output'
    write (stdout,*) '  -m      : Max number of steps'
    write (stdout,*) '  -n      : Number of steps for running average'
    write (stdout,*) '  -v      : Volume overide'
    write (stdout,*) '  -p      : Pressure'
    write (stdout,*) '  -f      : Fix unit cell'
    write (stdout,*) '  -mc     : Map to cell'
    write (stdout,*) '  -c      : Cluster (no PBC)'
    write (stdout,*) '  -tol    : Convergence thresold (eV)'
    write (stdout,*) '  -r      : Hard core radius (ang)'
    write (stdout,*) '  -a      : Hard core strength'
    write (stdout,*) '  -rmin   : Minimum tolerated contact (ang)'
    write (stdout,*) '  -ncell  : Number of cells in supercell'
    write (stdout,*) '  -natom  : Number of atoms in supercell'
    write (stdout,*) '  -scmax  : Maximum attempts to construct supercell'
    write (stdout,*) '  -o      : Optimise before and after'
    write (stdout,*) '  -niggli : Niggli reduce during trajectory'
    write (stdout,*) '  -h      : Display this message'
    write (stdout,*)

    stop

  end subroutine get_arguments

  ! *******************************************************************************************

end program ramble
