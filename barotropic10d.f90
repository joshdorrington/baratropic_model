
!Joshua Dorrington 22/03/18 - University of Oxford
!this module contains the model equations and integrates
! them using a first order Euler Maruyama scheme

module barotropic10d
	use coeffs10d
	use params
	use utils, only: white_noise, red_noise
	implicit none
    	private
	public run_model


	contains
		!takes an initial condition, integrates it forwards by dt stepnum times,
		!sampling at samplenum evenly spaced increments.
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
			real(dp) :: a1,a2,a3,a4,B1,B2,B3,B4,d1,d2,d3,d4,p1,p2,p3, &
									g1,g2,gstar1,gstar2,gprime1,gprime2,e1,e2

			!Unpacks the model coefficients

			a1=coeff(1) !a_11
			a2=coeff(2) !a_12
			a3=coeff(3) !a_21
			a4=coeff(4) !a_22

			B1=coeff(5) !B_11
			B2=coeff(6) !B_12
			B3=coeff(7) !B_21
			B4=coeff(8) !B_22

			d1=coeff(9)  !delta_11
			d2=coeff(10) !delta_12
			d3=coeff(11) !delta_21
			d4=coeff(12) !delta_22

			p1=coeff(13) !rho_11
			p2=coeff(14) !rho_12
			p3=coeff(15) !rho_21

			g1=coeff(16) !gamma_11
			g2=coeff(17) !gamma_12

			gstar1=coeff(18) !gamma_star_11
			gstar2=coeff(19) !gamma_star_12

			gprime1=coeff(20) !gamma_prime_12
			gprime2=coeff(21) !gamma_prime_21

			e1=coeff(22) !epsilon_1
			e2=coeff(23) !epsilon_2


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
				else if (noise_type=="d") then
					stoch_arr=0._dp
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

				!adds non linear triad terms

				dxdt=dxdt+ &
				(/0._dp, & !x1
				 -a1*x(1)*x(3) -d1*x(4)*x(6) +toggle*(-p1*(x(5)*x(8)-x(6)*x(7))),  & !x2
				  a1*x(1)*x(2) +d1*x(4)*x(5) +toggle*(p1*(x(5)*x(7)+x(6)*x(8))),  & !x3
				  e1*(x(2)*x(6)-x(3)*x(5))   +toggle*(e2*(x(7)*x(10)-x(8)*x(9))), & !x4
				 -a2*x(1)*x(6) -d2*x(3)*x(4) +toggle*(p2*(x(2)*x(8)-x(3)*x(7))),  & !x5
				  a2*x(1)*x(5) +d2*x(2)*x(4) +toggle*(-p2*(x(2)*x(7)+x(3)*x(8))),  & !x6
				 toggle*(-a3*x(1)*x(8) -d3*x(4)*x(10)-p3*(x(2)*x(6)+x(3)*x(5))),  & !x7
				 toggle*(a3*x(1)*x(7) +d3*x(4)*x(9) +p3*(x(2)*x(5)-x(3)*x(6))),  & !x8
				 toggle*(-a4*x(1)*x(10)-d4*x(4)*x(8)), & !x9
				 toggle*(a4*x(1)*x(9) +d4*x(4)*x(7))/) !x10


			end function dxdt
		end subroutine run_model




end module barotropic10d
