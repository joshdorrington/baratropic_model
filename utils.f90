
!Joshua Dorrington 22/03/18 - University of Oxford

!This module handles random number generation for stochastic runs
module utils
	use, intrinsic :: iso_fortran_env

	implicit none
	private
	public time_seed, white_noise
	   
	integer, parameter :: dp = real64
	contains
		!lifted wholesale from Sam Hatfield's Lorenz63-4D Var model
		subroutine time_seed()
            		integer :: i, n, clock
            		integer, allocatable :: seed(:)

            		call random_seed(size = n)
            		allocate(seed(n))

            		call system_clock(count=clock)

            		seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            		call random_seed(put = seed)

            		deallocate(seed)
		end subroutine

		
		subroutine white_noise(rand_arr,mean,std,n)
			real(dp) :: rand_arr(:)
			real(dp), intent(in) :: mean, std
			integer , intent(in):: n
			real(dp) :: rand1(n), rand2(n),u(n),v(n)

			!random univariates
			call random_number(rand1)
			call random_number(rand2)
			
			!Box-Muller
			u = sqrt(-2._dp * log(rand1))
           		v =   cos(8._dp*ATAN(1._dp) * rand2)
			rand_arr = mean + (std * u * v)
		end subroutine

end module utils
