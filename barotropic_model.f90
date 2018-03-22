!Joshua Dorrington 22/03/18 - University of Oxford
!This is the main executable program.

!This model solves a 6d spectral truncation of a baratropic beta plane model as
! first presented in [Charney & DeVore 1979] 
!https://doi.org/10.1175/1520-0469(1979)036<1205:MFEITA>2.0.CO;2

!all model parameters are in params
!formulas for all spectral coefficients are in coeffs
!random number generation is in utils
!the model equations and integrator are in baratropic6d

program barotropic_model
	use, intrinsic :: iso_fortran_env
	use params
	use coeffs, only: get_coeffs, build_lin_op
	use utils, only: time_seed, white_noise
	use barotropic6d, only: run_model
	implicit none
	
	
	!Initialise data structures
	real(dp), dimension(sample_num,dims) :: state_vector
	real(dp), dimension(coeff_num) :: coeff
	real(dp), dimension(dims,dims) :: lin_op
	
	!Seed random numbers
	call time_seed()

	!Generate coefficients and the linear time evolution operator
	call get_coeffs(b,beta,g,coeff)
	call build_lin_op(lin_op,coeff)
	

	print*,"initialisation complete"
	call run_model(init_con,step_num,sample_num,state_vector,coeff,lin_op)
	
	OPEN( 10, FILE=save_file )
	WRITE(10, *) state_vector   
end program
