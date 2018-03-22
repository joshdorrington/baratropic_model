
!Joshua Dorrington 22/03/18 - University of Oxford
!This module contains all the parameters, physical and otherwise
module params
	
	use, intrinsic :: iso_fortran_env

    	implicit none
	!global constants
	integer, parameter :: dp = real64
	real(dp), parameter :: pi = 4._dp*ATAN(1._dp)
	
	!model timestepping params
	real(dp), parameter :: dt =0.0002_dp
	real(dp), parameter :: t_0=0._dp
	real(dp), parameter :: t_f=1000._dp
	real(dp), parameter :: sample_t =1._dp

	!derived parameters, don't touch directly
	integer, parameter :: step_num=int((t_f-t_0)/dt)
	integer, parameter :: sample_num=int((t_f-t_0)/sample_t)
	
	!model params
	integer, parameter :: dims=6
	integer, parameter :: coeff_num =11
	real(dp), parameter, dimension(dims) :: &
	init_con =(/0.927796_dp,0.160574_dp,-0.0969953_dp,-0.672758_dp,-0.118184_dp,0.214505_dp/)
	real(dp), parameter, dimension(dims) :: &
	xf=(/0.95_dp,0._dp,0._dp,-0.76095_dp,0._dp,0._dp/)

	!forcing magnitude
	real(dp), parameter :: C=0.1_dp
	!strength of coriolis terms
	real(dp), parameter :: beta=1.25_dp
	!geometry of box; ratio of channel length to width	
	real(dp), parameter :: b=0.5_dp
	!topographical height parameter
	real(dp), parameter :: g=0.2_dp
	
	!I/O parameters
	character (len=9) :: save_file="state.out"

end module params
