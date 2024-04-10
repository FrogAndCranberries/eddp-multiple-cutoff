!==================================================================================!
!                                     nn_poly                                      !
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
! This program tests nn by fitting a polynomial function                           !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
program nn_poly

  use constants
  use nn

  implicit none

  type(network) :: net,newnet

  type(dataset) :: testing,validation,training

  integer, parameter :: ntrain=11,nvalid=7,ntest=101,unit_plot=21
  
  integer :: i
  
  real(kind=gp) xtrain(1,ntrain),ytrain(1,ntrain),xvalid(1,nvalid),yvalid(1,nvalid),xtest(1,ntest),ytest(1,ntest)
  
  ! -------------------------------------

  write (stderr,*) 'nn_poly test'
  
  ! ** Create network

  call nn_create(net,[1,10,1],actfn='tanh')
    
  write (stderr,*) 'Number of layers:  ',net%depth
  write (stderr,*) 'Number of weights: ',net%nwgts

  ! ** Generate training set

  do i=1,ntrain
     xtrain(1,i)=15.0_gp*(real(i-1,gp)/real(ntrain-1,gp)-0.0_gp)
     ytrain(1,i)=polyfunc(xtrain(1,i))
  end do

  ! ** Generate training set
  
  do i=1,nvalid
     xvalid(1,i)=15.0_gp*(real(i-1,gp)/real(nvalid-1,gp)-0.0_gp)
     yvalid(1,i)=polyfunc(xvalid(1,i))
  end do

  ! ** Generate testing set

  do i=1,ntest
     xtest(1,i)=15.0_gp*(real(i-1,gp)/real(ntest-1,gp)-0.0_gp)
     ytest(1,i)=polyfunc(xtest(1,i))
  end do

  ! ** Create data sets
  
  call nn_data(net,xtrain,ytrain,training)
  call nn_data(net,xvalid,yvalid,validation)  
  call nn_data(net,xtest,ytest,testing)

  ! ** Normalise, according to training set

  call nn_normalise(net,training)
  call nn_normalise(net,validation)
  call nn_normalise(net,testing)
  
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
  
  write (stderr,*) 'Cost for training set: ',nn_cost(newnet,training)
  write (stderr,*) 'Cost for validation set: ',nn_cost(newnet,validation)  
  write (stderr,*) 'Cost for testing set:  ',nn_cost(newnet,testing)
  
  ! ** Plot results

  open(unit_plot,file='poly.plot',form='formatted',status='unknown')

  do i=1,ntrain
     write (unit_plot,*) training%in(1,i),training%out(:,i)
  end do
  write (unit_plot,*)
  do i=1,ntest
     write (unit_plot,*) testing%in(1,i),testing%out(:,i)
  end do
  write (unit_plot,*)
  do i=1,ntest
     write (unit_plot,*) testing%in(1,i),nn_forward(newnet,testing%in(:,i))
  end do

  close(unit_plot)
  
contains

  function polyfunc(x) result(y)
    
    real(kind=gp), intent(in) :: x

    real(kind=gp) :: y

    y=(x**6-36*x**5+450*x**4-2400*x**3+5400*x**2-4320*x+720)/720

    !y=log(abs(y))

    
  end function polyfunc
  
end program nn_poly
