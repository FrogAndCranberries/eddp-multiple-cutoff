!==================================================================================!
!                                     nn_sine                                      !
!==================================================================================!
!                                                                                  !
! This file is part of the nn package.                                             !
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
! This program tests nn by fitting a sinusoidal function                           !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
program nn_sine

  use constants
  use nn

  implicit none

  type(network) :: net,newnet

  type(dataset) :: training,validation,testing

  integer, parameter :: ntrain=12,nvalid=7,ntest=101,unit_plot=21
  
  integer :: i
  
  real(kind=gp) xtrain(1,ntrain),ytrain(1,ntrain),xvalid(1,nvalid),yvalid(1,nvalid),xtest(1,ntest),ytest(1,ntest)
  
  ! -------------------------------------

  write (stderr,*) 'nn_sine test'
  
  ! ** Create network

  call nn_create(net,[1,5,1],actfn='tanh')
  
  write (stderr,*) 'Number of layers:  ',net%depth
  write (stderr,*) 'Number of weights: ',net%nwgts

  ! ** Generate training set

  do i=1,ntrain
     xtrain(1,i)=2.0_gp*(real(i-1,gp)/real(ntrain-1,gp)-0.5_gp)+5.0_gp
     ytrain(1,i)=sin(tpi*xtrain(1,i))+2.0_gp
  end do

  call nn_data(net,xtrain,ytrain,training)

  ! ** Normalise based on training data
  
  call nn_normalise(net,training,outonly=.false.)

   ! ** Generate validation set

  do i=1,nvalid
     xvalid(1,i)=2.0_gp*(real(i-1,gp)/real(nvalid-1,gp)-0.5_gp)+5.0_gp
     yvalid(1,i)=sin(tpi*xvalid(1,i))+2.0_gp
  end do

  call nn_data(net,xvalid,yvalid,validation)

  ! ** Normalise based on training data
  
  call nn_normalise(net,validation,outonly=.false.)
  
  ! ** Generate testing set

  do i=1,ntest
     xtest(1,i)=2.0_gp*(real(i-1,gp)/real(ntest-1,gp)-0.5_gp)+5.0_gp
     ytest(1,i)=sin(tpi*xtest(1,i))+2.0_gp
  end do

  call nn_data(net,xtest,ytest,testing)

  ! ** Normalise based on training data
  
  call nn_normalise(net,testing,outonly=.false.)

  if(.false.) then
     call nn_test_gradient(net,training)
     stop
  end if
  
  ! ** Train network

  call nn_train_galm(net,training,validation,thresh=epsilon(1.0_gp),track=.true.)
  
  ! ** Save trained network
  
  call nn_save(net,'trained.net')

  ! ** Load saved network
  
  call nn_load(newnet,'trained.net')
  
  ! ** Report costs
  
  write (stderr,*) 'Cost for training set:   ',nn_cost(newnet,training)
  
  write (stderr,*) 'Cost for validation set: ',nn_cost(newnet,validation)
  
  write (stderr,*) 'Cost for testing set:    ',nn_cost(newnet,testing)

  ! ** Plot results

  open(unit_plot,file='sine.plot',form='formatted',status='unknown')
  
  do i=1,ntrain
     write (unit_plot,*) xtrain(1,i),ytrain(1,i)
  end do
  write (unit_plot,*)
  do i=1,ntest
     write (unit_plot,*) xtest(1,i),ytest(1,i)
  end do
  write (unit_plot,*)
  do i=1,ntest
     write (unit_plot,*) xtest(1,i),nn_forward(newnet,xtest(1:1,i),scale=.true.)
  end do

  close(unit_plot)
  
end program nn_sine
