!==================================================================================!
!                                      md                                          !
!==================================================================================!
!                                                                                  !
!----------------------------------------------------------------------------------!
! This module performs the molecular dynamics                                      !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module md

  use constants
  use rng
  use cell
  use ddp
  use spglib_f08

  implicit none

  private

  public :: md_cell

  logical,       public :: map_to_cell=.false.
  real(kind=gp), public :: time_step
  real(kind=gp), public :: damping_rate
  real(kind=gp), public :: damping_cell
  real(kind=gp), public :: initial_temp
  real(kind=gp), public :: target_temp
  real(kind=gp), public :: heat_rate
  real(kind=gp), public :: quench_rate
  real(kind=gp), public :: quench_time  
  real(kind=gp), public :: quench_end
  real(kind=gp), public :: stop_time  
  integer,       public :: num_run
  integer,       public :: unit_track=21
  integer,       public :: unit_temp=22
  integer,       public :: unit_press=23
  integer,       public :: unit_quench=24
  integer,       public :: unit_stop=25
  
  real(kind=gp), public :: md_time=0.0_gp
  
  logical, public :: niggli=.false.
  
  !---------------------------------------------------!

  character(len=80)     :: trackfile,xyzefile,xyzerunfile
  
contains

  subroutine md_cell(str,steps,cell_fix)

    type(structure), intent(inout) :: str

    !-------------------------------------------------------------!

    integer,       intent(inout) :: steps
    logical,       intent(in)    :: cell_fix

    real(kind=gp) :: e,ke,p,sig(3,3),s(3,3),lat0(3,3),strain(3,3),temp0,p0(3),H0,trun,hrun,erun,vrun,prun,com(3),Ic(3,3)
    real(kind=gp) :: msd,sigrun(3,3),srun(3,3),lang(3),qi(3),pinst,tinst,tt_save

    real(kind=gp) :: car(3,3)
    real(kind=dp) :: car_dp(3,3)
    
    real(kind=gp) :: ion_positions_run(3,str%num_ions),ion_positions_run0(3,str%num_ions)

    real(kind=gp), save, allocatable :: g(:,:),v(:),x(:),xp(:),xp0(:),mass(:),x0(:)

    real(kind=gp) :: start_time,end_time,elapsed_time=0.0_gp
    
    integer :: ni,maxsteps,info

    logical :: tempfile
    logical :: pressfile
    logical :: quenchfile
    logical :: stopfile

    ! * Clean up steering files

    if(.not.quiet) then

    open(unit=unit_temp,file="TEMP",err=100)
100 close(unit_temp,status='delete')
    open(unit=unit_press,file="PRESS",err=101)
101 close(unit_press,status='delete')
    open(unit=unit_quench,file="QUENCH",err=102)
102 close(unit_quench,status='delete')
    open(unit=unit_stop,file="STOP",err=103)
103 close(unit_stop,status='delete')

    end if

    maxsteps=steps

    if(.not.allocated(g)) &
         allocate(g(3,str%num_ions),v(3*str%num_ions+3*3),x(3*str%num_ions+3*3),&
         xp(3*str%num_ions+3*3),xp0(3*str%num_ions+3*3),mass(3*str%num_ions+9))

    ! ** Open a file for the convergence data

    if(track) then

       write(trackfile,'(a)') trim(seedname)//".track"       
       open(unit=unit_track,file=trackfile,status="unknown",err=998)

       if(.not.quiet) then
          write(xyzefile,'(a)') trim(seedname)//".xyze"       
          open(unit=unit_xyze,file=xyzefile,status="unknown",err=998)
          call write_xyze(str)
          write(xyzerunfile,'(a)') trim(seedname)//"-run.xyze"       
          open(unit=unit_xyze_run,file=xyzerunfile,status="unknown",err=998)
       end if

    end if

    p=(str%external_pressure(1,1)+str%external_pressure(2,2)+str%external_pressure(3,3))/3.0_gp

    do ni=1,str%num_ions
       mass((ni-1)*3+1:(ni-1)*3+3)=str%ion_mass(ni)
    end do

    mass(3*str%num_ions+1:3*str%num_ions+9)=10*sum(str%ion_mass)
    
    ! * Initialise velocities

    call random_number(v(1:3*str%num_ions))
    v(1:3*str%num_ions)=(v(1:3*str%num_ions)-0.5)*2
    v(3*str%num_ions+1:3*str%num_ions+9)=0.0_gp ! zero the cell velocities

    ! * Set overall momentum to zero

    p0=0.0_gp
    do ni=1,str%num_ions
       p0(1:3)=p0(1:3)+v((ni-1)*3+1:(ni-1)*3+3)*mass((ni-1)*3+1:(ni-1)*3+3)
    end do

    p0=p0/str%num_ions

    do ni=1,str%num_ions
       v((ni-1)*3+1:(ni-1)*3+3)=v((ni-1)*3+1:(ni-1)*3+3)-p0(1:3)/mass((ni-1)*3+1:(ni-1)*3+3)
    end do

    ! * Set overall angular momentum to zero

    if(cluster) then

       com=0.0_gp
       do ni=1,str%num_ions
          com=com+str%ion_positions(:,ni)*str%ion_mass(ni)
       end do
       com=com/sum(str%ion_mass)

       lang=0.0_gp
       Ic=0.0_gp
       do ni=1,str%num_ions
          qi=str%ion_mass(ni)*v((ni-1)*3+1:(ni-1)*3+3)
          lang=lang+cross3(str%ion_positions(:,ni)-com,qi)
          Ic=Ic+str%ion_mass(ni)*ccmat(str%ion_positions(:,ni)-com)
       end do

       do ni=1,str%num_ions
          qi=cross3(matmul(inv(Ic),lang),str%ion_positions(:,ni)-com)
          v((ni-1)*3+1:(ni-1)*3+3)=v((ni-1)*3+1:(ni-1)*3+3)-qi          
       end do

    endif

    ! ** Rescale, assuming starting from relaxed state

    temp0=sum(mass(1:3*str%num_ions)*v(1:3*str%num_ions)**2)/size(v(1:3*str%num_ions))*ev2k

    v(1:3*str%num_ions)=v(1:3*str%num_ions)*sqrt(2*initial_temp/temp0)

    ! * Compute kinetic energy

    ke=sum(mass(1:3*str%num_ions)*v(1:3*str%num_ions)**2)/2

    ! ** Start up

    lat0=str%lattice_car

    strain=0.0_gp

    x=[reshape(str%ion_positions,(/3*str%num_ions/)),reshape(strain,(/3*3/))]

    x0=x

    s=0.0_gp
    call eval_ddp(str,e,g,sig)
    if(.not.cell_fix) then
       do ni=1,str%num_ions
          sig=sig+str%ion_mass(ni)*outer3(v((ni-1)*3+1:(ni-1)*3+3),v((ni-1)*3+1:(ni-1)*3+3))
       end do
       s = sig-str%external_pressure*str%volume
    end if
    
    H0=e+sum(mass(1:3*str%num_ions)*v(1:3*str%num_ions)**2)/2+p*str%volume

    ! ** Start up running averages

    trun=initial_temp!2*ke/size(v(1:3*str%num_ions))*ev2k
    hrun=(e+p*str%volume+ke)/str%num_ions
    erun=(e+ke)/str%num_ions
    vrun=str%volume/str%num_ions
    prun=(sig(1,1)+sig(2,2)+sig(3,3))/str%volume/3*evbyang3
    srun=s
    sigrun=sig
    ion_positions_run=str%ion_positions
    ion_positions_run0=str%ion_positions
    car=str%lattice_car

    xp=[reshape(matmul(ident+strain,g),(/3*str%num_ions/)),reshape(matmul(s,transpose(inv(ident+strain))),(/3*3/))]

    tt_save=target_temp

    steps = 0
    do while(steps.lt.maxsteps)

       if (.not.quiet) then
       inquire(file="STOP",exist=stopfile)
       if(stopfile) exit
       end if
       
       if(target_temp.lt.quench_end) exit
       
       if(elapsed_time.gt.stop_time) exit

       ! * Track convergence
       
       msd=sqrt(sum(((ion_positions_run-ion_positions_run0)/(time_step*ntu2ps))**2)/str%num_ions)

       if((mod(total_steps,track_every).eq.0).and.track) then

          ! ** Convert time to picoseconds (check - factor of 2 or sqrt2?), temperature to kelvin

          write (unit_track,'(9(1x,f20.10))') &
               elapsed_time,&
               str%volume/str%num_ions,&
               2*ke/size(v(1:3*str%num_ions))*ev2k,&
               target_temp,&
               (sig(1,1)+sig(2,2)+sig(3,3))/str%volume/3*evbyang3,&
               (e+p*str%volume)/str%num_ions,&
               p*evbyang3,&
               (e+ke+p*str%volume)/str%num_ions,&
               msd

          call flush(unit_track)

          call write_xyze(str)

          call write_xyze_run(str,car,ion_positions_run)

          if (.not.quiet) then

          inquire(file="TEMP",exist=tempfile)

          if(tempfile) then

             open(unit_temp,file="TEMP",form='formatted',err=998)

             read(unit_temp,*) target_temp

             close(unit_temp)

          end if

          inquire(file="PRESS",exist=pressfile)

          if(pressfile) then

             open(unit_press,file="PRESS",form='formatted',err=998)

             read(unit_press,*) p

             str%external_pressure=0.0_gp

             p=p/evbyang3

             str%external_pressure(1,1)=p
             str%external_pressure(2,2)=p
             str%external_pressure(3,3)=p

             close(unit_press)

          end if

          inquire(file="QUENCH",exist=quenchfile)

          if(quenchfile) then

             open(unit_quench,file="QUENCH",form='formatted',err=998)

             read(unit_quench,*) quench_rate
             
             quench_time=0.0_gp

             close(unit_quench)

          end if
          
          end if

       end if

       ! * Increment counters

       steps=steps+1
       total_steps=total_steps+1
       elapsed_time=total_steps*time_step*ntu2ps

       call cpu_time(start_time)

       xp0=xp

       ! * Integate

       x=x+v*time_step+xp0*time_step**2/mass/2

       ! * Update strain, unit cell and positions

       strain=reshape(x(3*str%num_ions+1:3*str%num_ions+3*3),(/3,3/))

       ! * Remove rotation

       strain=(strain+transpose(strain))/2
              
       if(fix_vol) then

         stop 'fix_vol implemention not completed'
          
       end if

       car_dp=transpose(matmul(ident+strain,lat0))

       if(niggli) info=spg_niggli_reduce(car_dp,1e-4_dp)

       str%lattice_car=transpose(car_dp)
              
       str%ion_positions=matmul(ident+strain,reshape(x(1:3*str%num_ions),(/3,str%num_ions/)))

       call update_cell(str)

       ! * Map to cell if requested
       
       if(map_to_cell) then
          do ni=1,str%num_ions
             str%ion_positions(:,ni) = matmul(str%lattice_rec,str%ion_positions(:,ni)) 
             str%ion_positions(:,ni) = str%ion_positions(:,ni) - real(floor(str%ion_positions(:,ni),8),gp)
             str%ion_positions(:,ni) = matmul(str%lattice_car,str%ion_positions(:,ni))
          end do
       end if

       ! * Compute gradient

       call eval_ddp(str,e,g,sig)
       
       if(.not.cell_fix) then 
          do ni=1,str%num_ions
             sig=sig+str%ion_mass(ni)*outer3(v((ni-1)*3+1:(ni-1)*3+3),v((ni-1)*3+1:(ni-1)*3+3))      
          end do
          s = sig-str%external_pressure*str%volume
       end if

       xp=[reshape(matmul(ident+strain,g),(/3*str%num_ions/)),reshape(matmul(s,transpose(inv(ident+strain))),(/3*3/))]

       ! ** Update velocities

       v=v+(xp+xp0)*time_step/mass/2

       ! * Compute kinetic energy

       ke=sum(mass(1:3*str%num_ions)*v(1:3*str%num_ions)**2)/2

       ! * Compute instantaneous pressure 

       pinst=(sig(1,1)+sig(2,2)+sig(3,3))/str%volume/3*evbyang3
       
       ! * Compute instantaneous temperature

       tinst=2*ke/size(v(1:3*str%num_ions))*ev2k

       ! * Heat/Quench

       if(elapsed_time.gt.quench_time) then
          if(quench_rate.gt.0.0_gp) then
             target_temp=target_temp-quench_rate*time_step*ntu2ps
          else
             if(quench_time.gt.0.0_gp) exit
          end if
          
       else 
         if(heat_rate.gt.0.0_gp) target_temp=min(tt_save,initial_temp+steps*heat_rate*time_step*ntu2ps)
       end if

       ! * Damp

       if(damping_rate.gt.0.0_gp) then
          v(1:3*str%num_ions)=v(1:3*str%num_ions)*(1-damping_rate*(tinst-target_temp)/tinst)
          v(3*str%num_ions+1:)=v(3*str%num_ions+1:)*(1-damping_cell*(pinst-p)/pinst)
       else
          v(1:3*str%num_ions)=v(1:3*str%num_ions)*(abs(H0-e-p*str%volume)/ke)**0.5_gp
          v(3*str%num_ions+1:)=v(3*str%num_ions+1:)*(1-damping_cell*(pinst-p)/pinst)
       end if

       ! * Set overall momentum to zero

       p0=0.0_gp
       do ni=1,str%num_ions
          p0(1:3)=p0(1:3)+v((ni-1)*3+1:(ni-1)*3+3)*mass((ni-1)*3+1:(ni-1)*3+3)
       end do

       p0=p0/str%num_ions

       do ni=1,str%num_ions
          v((ni-1)*3+1:(ni-1)*3+3)=v((ni-1)*3+1:(ni-1)*3+3)-p0(1:3)/mass((ni-1)*3+1:(ni-1)*3+3)
       end do

       ! * Set overall angular momentum to zero

       if(cluster) then

          com=0.0_gp
          do ni=1,str%num_ions
             com=com+str%ion_positions(:,ni)*str%ion_mass(ni)
          end do
          com=com/sum(str%ion_mass)

          lang=0.0_gp
          Ic=0.0_gp
          do ni=1,str%num_ions
             qi=str%ion_mass(ni)*v((ni-1)*3+1:(ni-1)*3+3)
             lang=lang+cross3(str%ion_positions(:,ni)-com,qi)
             Ic=Ic+str%ion_mass(ni)*ccmat(str%ion_positions(:,ni)-com)
          end do

          do ni=1,str%num_ions
             qi=cross3(matmul(inv(Ic),lang),str%ion_positions(:,ni)-com)
             v((ni-1)*3+1:(ni-1)*3+3)=v((ni-1)*3+1:(ni-1)*3+3)-qi          
          end do

       endif

       ! * Recompute kinetic energy

       ke=sum(mass(1:3*str%num_ions)*v(1:3*str%num_ions)**2)/2

       ! ** Update Running averages

       trun=trun+((2*ke/size(v(1:3*str%num_ions))*ev2k)-trun)/num_run
       hrun=hrun+((e+p*str%volume+ke)/str%num_ions-hrun)/num_run
       erun=erun+((e+ke)/str%num_ions-erun)/num_run
       vrun=vrun+(str%volume/str%num_ions-vrun)/num_run
       srun=srun+(s-srun)/num_run
       prun=prun+(pinst-prun)/num_run
       sigrun=sigrun+(sig-sigrun)/num_run

       car=car+(str%lattice_car-car)/num_run
       
       ion_positions_run0=ion_positions_run
       
       ion_positions_run=ion_positions_run+(str%ion_positions-ion_positions_run)/num_run

       ! ** Track running averages

       if((mod(total_steps,track_every).eq.0).and.track) &
             write (stderr,'(a,*(1x,f0.5))') trim(space_group(str,car,ion_positions_run)),prun,trun,vrun,erun,hrun,elapsed_time

       call cpu_time(end_time)

       md_time=md_time+end_time-start_time

    end do

    if(track) call write_xyze(str)

    str%volume=vrun*str%num_ions
    str%stress=sigrun/str%volume
    str%energy=erun*str%num_ions
    str%enthalpy=hrun*str%num_ions

    return

998 stop 'There is a problem opening the track/xyze/TEMP/PRESS file. Stopping.'

  end subroutine md_cell

  function space_group(str,car,pos)
     
     type(structure), intent(in) :: str
     real(kind=gp), intent(in) :: car(3,3)
     real(kind=gp), intent(in) :: pos(3,str%num_ions)
     
     character(len=20) :: space_group
     
     type(SpglibDataset) :: dset
     
     real(kind=gp) :: trans(3,3),rec(3,3),frac(3,str%num_ions)
     
     integer :: ni
     
     space_group='---'
     
     if(.not.symmgen) return
     
     ! ** spglib uses the transpose of lattice_car

     trans=transpose(car)

     ! ** Calculate the reciprocal lattice vectors

     rec(1,1)=car(2,2)*car(3,3)-car(3,2)*car(2,3)
     rec(2,1)=car(2,3)*car(3,1)-car(3,3)*car(2,1)
     rec(3,1)=car(2,1)*car(3,2)-car(3,1)*car(2,2)
     rec(1,2)=car(3,2)*car(1,3)-car(1,2)*car(3,3)
     rec(2,2)=car(3,3)*car(1,1)-car(1,3)*car(3,1)
     rec(3,2)=car(3,1)*car(1,2)-car(1,1)*car(3,2)
     rec(1,3)=car(1,2)*car(2,3)-car(2,2)*car(1,3)
     rec(2,3)=car(1,3)*car(2,1)-car(2,3)*car(1,1)
     rec(3,3)=car(1,1)*car(2,2)-car(2,1)*car(1,2)

     rec(:,:)=rec(:,:)/str%volume ! ** IS THIS THE RIGHT VOLUME

     ! ** fractional positions

     do ni=1,str%num_ions
        frac(:,ni)=matmul(rec,pos(:,ni)) 
     end do

     dset=spg_get_dataset(real(trans,dp),real(frac,dp),str%ion_species,str%num_ions,0.1_dp)

     ! ** Get the international symbol for the space group based on the provided tolerance

     if(dset%spacegroup_number/=0) space_group=deunder(trim(dset%international_symbol))
     
     
  end function space_group

  subroutine write_xyze_run(str,lcr,ipr)

    type(structure), intent(in) :: str
    real(kind=gp), dimension(3,3) :: lcr
    real(kind=gp), dimension(3,str%num_ions) :: ipr
    
    integer :: ni
    character(len=10) :: ctemp

    if(quiet) return
    
    write (ctemp,'(i10)') str%num_ions
    write (unit_xyze_run,'(a)') trim(adjustl(ctemp))

    write (unit_xyze_run,'(a,9f20.6,a)') 'Lattice="',reshape(lcr,(/9/)),'" Properties=species:S:1:pos:R:3'
    do ni=1,str%num_ions
       write (unit_xyze_run,'(a,3f20.13)') trim(adjustl(str%ion_names(ni))),ipr(:,ni)
    end do

    flush(unit_xyze_run)
    
  end subroutine write_xyze_run
  
  function inv(A) result(Ainv)

    real(kind=gp), dimension(:,:), intent(in) :: A
    real(kind=gp), dimension(size(A,1),size(A,2)) :: Ainv

    real(kind=gp), dimension(size(A,1)) :: work 
    integer,       dimension(size(A,1)) :: ipiv
    integer :: n, info

    external dgetrf,sgetrf
    external dgetri,sgetri

    Ainv=A
    n=size(A,1)

    if(gp.eq.sp) then
       call sgetrf(n,n,Ainv,n,ipiv,info)
    else
       call dgetrf(n,n,Ainv,n,ipiv,info)
    end if
    
    if (info.ne.0)  stop 'inv: matrix is numerically singular!'

    if(gp.eq.sp) then
       call sgetri(n,Ainv,n,ipiv,work,n,info)
    else
       call dgetri(n,Ainv,n,ipiv,work,n,info)
    end if
    
    if (info.ne.0) stop 'inv: matrix inversion failed!'

  end function inv

  function outer3(a,b) result(c)

    real(kind=gp), intent(in) :: a(3),b(3)

    real(kind=gp) :: c(3,3)

    c(1,1)=a(1)*b(1)
    c(2,1)=a(1)*b(2)
    c(3,1)=a(1)*b(3)

    c(1,2)=a(2)*b(1)
    c(2,2)=a(2)*b(2)
    c(3,2)=a(2)*b(3)

    c(1,3)=a(3)*b(1)
    c(2,3)=a(3)*b(2)
    c(3,3)=a(3)*b(3)
    
  end function outer3

  function cross3(a,b) result(c)
    
    real(kind=gp), intent(in) :: a(3),b(3)

    real(kind=gp) :: c(3)

    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    
  end function cross3
  
  function ccmat(a) result(c)
    
    real(kind=gp), intent(in) :: a(3)

    real(kind=gp) :: c(3,3)

    c(1,1)=a(2)*a(2)+a(3)*a(3)
    c(2,2)=a(1)*a(1)+a(3)*a(3)
    c(3,3)=a(1)*a(1)+a(2)*a(2)
    c(1,2)=-a(1)*a(2)
    c(2,1)=c(1,2)
    c(1,3)=-a(1)*a(3)
    c(3,1)=c(1,3)
    c(2,3)=-a(2)*a(3)
    c(3,2)=c(2,3)
    
  end function ccmat
  
end module md
