!==================================================================================!
!                                  nn_fourier                                      !
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
! This program tests nn by learning fourier transforms                             !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2020-2023                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!
program nn_fourier

  use constants
  use nn

  implicit none

  type(network) :: net

  type(dataset) :: testing,training,validation

  integer, parameter :: ntrain=201,nvalid=101,ntest=1,unit_plot=21
  integer, parameter :: nin=10,nout=10

  integer :: i,n,m

  real(kind=gp) :: xtrain(nin,ntrain),ytrain(nout,ntrain),xtest(nin,ntest),ytest(nout,ntest),xvalid(nin,nvalid),yvalid(nout,nvalid)

  real(kind=gp) :: x

  ! -------------------------------------

  write (stderr,*) 'nn_fourier test'

  ! ** Create network

  call nn_create(net,[nin,10,nout],actfn='tanh')

  write (stderr,*) 'Number of inputs:  ',nin
  write (stderr,*) 'Number of outputs: ',nout
  write (stderr,*) 'Number of layers:  ',net%depth
  write (stderr,*) 'Number of weights: ',net%nwgts

  ! ** Generate training set

  ytrain=0.0_gp
  
  do i=1,ntrain

     call random_number(xtrain(:,i))

     do n=1,nout

        x=(real(n-1,gp)/real(nout-1,gp))

        do m=1,size(xtrain,1)

           ytrain(n,i)=ytrain(n,i)+xtrain(m,i)*cos((m-1)*pi*x)

        end do

     end do

  end do

  call nn_data(net,xtrain,ytrain,training)

  ! ** Generate validation set

  yvalid=0.0_gp
  
  do i=1,nvalid

     call random_number(xvalid(:,i))

     do n=1,nout

        x=(real(n-1,gp)/real(nout-1,gp)-0.0_gp)

        do m=1,size(xvalid,1)

           yvalid(n,i)=yvalid(n,i)+xvalid(m,i)*cos((m-1)*pi*x)

        end do

     end do

  end do

  call nn_data(net,xvalid,yvalid,validation)

  ! ** Generate testing set

  ytest=0.0_gp
  
  do i=1,ntest

     call random_number(xtest(:,i))

     do n=1,nout

        x=(real(n-1,gp)/real(nout-1,gp)-0.0_gp)

        do m=1,size(xtest,1)

           ytest(n,i)=ytest(n,i)+xtest(m,i)*cos((m-1)*pi*x)

        end do

     end do

  end do

  call nn_data(net,xtest,ytest,testing)

  ! ** Normalise based on training set

  call nn_normalise(net,training,outonly=.false.)
  call nn_normalise(net,validation,outonly=.false.)
  call nn_normalise(net,testing,outonly=.false.)
  
  if(.false.) then
     call nn_test_gradient(net,training)
     stop
  end if

  ! ** Train network

  call nn_train_galm(net,training,validation,thresh=epsilon(1.0_gp),track=.true.)
  
  ! ** Report costs

  write (stderr,*) 'Cost for training set:   ',nn_cost(net,training)
  
  write (stderr,*) 'Cost for validation set: ',nn_cost(net,validation)

  write (stderr,*) 'Cost for testing set:    ',nn_cost(net,testing)

  ! ** Plot results

  open(unit_plot,file='fourier.plot',form='formatted',status='unknown')
  
  do i=1,ntest

     write (unit_plot,'(f20.10)') testing%out(:,i)
     write (unit_plot,*)

     write (unit_plot,'(f20.10)') nn_forward(net,testing%in(:,i))
     write (unit_plot,*)

  end do
  write (unit_plot,*)


end program nn_fourier
