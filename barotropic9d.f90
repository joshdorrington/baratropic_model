
!Joshua Dorrington 22/03/18 - University of Oxford
!this module contains the model equations and integrates
! them using a first order Euler Maruyama scheme

module barotropic6d
	use coeffs
	use params
	use utils, only: white_noise, red_noise
	implicit none
    	private
	public run_model


	contains
		subroutine run_model(init_con,stepnum,samplenum,output, coeff,lin_op,n_inner)

			integer(kind=int64), intent(in) :: stepnum
			integer, intent(in) :: samplenum, n_inner
			real(dp), intent(in) :: init_con(dims)
			real(dp), dimension(samplenum, dims) :: output
			real(dp), dimension(dims) :: k1, k2,x
			real(dp),intent(in), dimension(coeff_num) :: coeff
			real(dp), intent(in), dimension(dims,dims) :: lin_op
			integer :: i, j
			real(dp), dimension(n_inner,dims) :: stoch_arr

			!for holding coeffs:
			real(dp),dimension(6) :: a
			real(dp),dimension(3) :: beta
			real(dp), dimension(4) :: d
			real(dp), dimension(4) :: k
			real(dp), dimension(4) :: e


			!Unpacks the model coefficients
			a=coeff(1:6)
			beta=coeff(7:9)
			d=coeff(20:23)
			k=coeff(24:27)
			e=coeff(28:31)

			! First step
	       		output(1,:) = init_con
			x=init_con

			!all other steps, sampling at end of inner loop
			do i = 1, size(output,1)-1

				!generate random numbers for next inner loop
				if (noise_type=="w") then
					call white_noise(stoch_arr,0._dp,sqrt(dt),n_inner)
				else if (noise_type=="r") then
					call red_noise(stoch_arr,0._dp,sqrt(dt),n_inner)
				end if

				do j = 1, n_inner
		    			k1 = dt*dxdt(x)
		    			!k2 = dt*dxdt(x + k1)

		    			!x = x + 0.5_dp * (k1 + k2)
					x=x+k1+sigma*stoch_arr(j,:)
				end do
				output(i+1,:)=x
				!print*,(100._dp*(i-1))/size(output,1), " percent complete"
			end do

		contains

			!The o.d.e is here
			function dxdt(x)

				real(dp),intent(in) :: x(dims)
				real(dp) :: dxdt(dims)

				!applies linear matrix operator
				dxdt=matmul(lin_op,x)+C*xf
				!adds non linear terms
				dxdt=dxdt+(/e(2)*(x(3)*x(8)-x(2)*x(9)),&
				-x(3)*(a(1)*x(1)+a(4)*x(7))-d(1)*x(4)*x(6)-x(9)*(k(1)*x(1)+k(3)*x(7)),&
				 x(2)*(a(1)*x(1)+a(4)*x(7))+d(1)*x(4)*x(5)+x(8)*(k(1)*x(1)+k(3)*x(7)),&
				 e(1)*(x(2)*x(6)-x(3)*x(5))+e(3)*(x(5)*x(9)-x(6)*x(8)), &
				-x(6)*(a(2)*x(1)+a(5)*x(7))-x(4)*(d(2)*x(3)+d(3)*x(9)), &
				 x(5)*(a(2)*x(1)+a(5)*x(7))+x(4)*(d(2)*x(2)+d(3)*x(8)), &
				 e(4)*(x(2)*x(9)-x(3)*x(8)),&
				-x(9)*(a(3)*x(1)+a(6)*x(7))-d(4)*x(4)*x(6)-x(3)*(k(4)*x(7)-k(2)*x(1)),&
				 x(8)*(a(3)*x(1)+a(6)*x(7))+d(4)*x(4)*x(5)+x(2)*(k(4)*x(7)-k(2)*x(1))/)

			end function dxdt
		end subroutine run_model




end module barotropic6d
