!Joshua Dorrington 22/03/18 - University of Oxford
!This module generates the coefficients of the barotropic model,
!and also generates the linear time evolution operator

module coeffs10d
    use, intrinsic :: iso_fortran_env
    use params
    implicit none
    private
    public get_coeffs, build_lin_op


    contains

	    subroutine get_coeffs(b,beta,g,coeff_arr)

    		real(dp) :: coeff_arr(23)
    		real(dp), intent(in) :: b, beta, g
    		coeff_arr(1)=a_coeff(1,1,b) !a_11
    		coeff_arr(2)=a_coeff(1,2,b) !a_12
        coeff_arr(3)=a_coeff(2,1,b) !a_21
        coeff_arr(4)=a_coeff(2,2,b) !a_22

    		coeff_arr(5)=beta_coeff(1,1,b,beta) !B_11
    		coeff_arr(6)=beta_coeff(1,2,b,beta) !B_12
    		coeff_arr(7)=beta_coeff(2,1,b,beta) !B_21
    		coeff_arr(8)=beta_coeff(2,2,b,beta) !B_22


    		coeff_arr(9)=d_coeff(1,1,b) !delta_11
    		coeff_arr(10)=d_coeff(1,2,b)!delta_12
    		coeff_arr(11)=d_coeff(2,1,b)!delta_21
    		coeff_arr(12)=d_coeff(2,2,b)!delta_22

    		coeff_arr(13)=p_coeff(1,1,b) !p_11
        coeff_arr(14)=p_coeff(1,2,b) !p_12
        coeff_arr(15)=p_coeff(2,1,b) !p_21

        coeff_arr(16)=g_coeff(1,1,b,g) !g_11
        coeff_arr(17)=g_coeff(1,2,b,g) !g_12

        coeff_arr(18)=g_star_coeff(1,1,b,g) !gstar_11
        coeff_arr(19)=g_star_coeff(1,2,b,g) !gstar_12

        coeff_arr(20)=g_prime_coeff(1,2,b,g) !gprime_12
        coeff_arr(21)=g_prime_coeff(2,1,b,g) !gprime_21

        coeff_arr(22)=e_coeff(1,b) !epsilon_1
        coeff_arr(23)=e_coeff(2,b) !epsilon_2

	    end subroutine get_coeffs

            subroutine build_lin_op(mat,coeff)
	        real(dp), dimension(coeff_num),intent(in) :: coeff
		real(dp), dimension(dims,dims) :: mat
		integer :: i, j

		do j=1, dims
		  do i=1, dims
			mat(i,j)=0._dp
		  end do
		end do

		do i=1, dims
		  mat(i,i)=-C
		end do
    !sets all the linear terms as matrix elements

		mat(1,3)= coeff(18) !gstar_11
		mat(2,3)= coeff(5)  !B_11
		mat(3,1)=-coeff(16) !-g_11
		mat(3,2)=-coeff(5)  !-B_11
		mat(4,6)= coeff(19) !g_star_12
		mat(5,6)= coeff(6)  !B_12
    mat(5,8)= toggle*coeff(20) !gamma_prime_12
		mat(6,4)=-coeff(17) !-gamma_12
		mat(6,5)=-coeff(6)  !-B_12
    mat(6,7)=toggle*(-coeff(20)) !-gamma_prime_12
    mat(7,6)= coeff(21) !gamma_prime_21
    mat(7,8)= coeff(7)  !B21
    mat(8,5)=-coeff(21) !-gamma_prime_21
    mat(8,7)=-coeff(7)  !-B_21
    mat(9,10)= coeff(8) !B_22
    mat(10,9)=-coeff(8) !-B_22



	    end subroutine build_lin_op

	    function a_coeff(n,m,b) result(a)

	    	integer, intent(in) :: m,n
		    real(dp), intent(in) :: b
	    	real(dp) :: a

        a=(8._dp*sqrt(2._dp)*n*(m**2)*(n*(b**2)+(m**2)-1)) / (pi*(-1._dp + 4._dp*m**2)*((m**2)+(n*b)**2))
	    end function a_coeff

	    function beta_coeff(n,m,b,beta) result(output)

	    	integer, intent(in) :: m,n
		    real(dp), intent(in) :: beta, b
	    	real(dp) :: output

	    	output=(beta*n*(b**2))/((n*b)**2 + m**2)

	    end function beta_coeff

      function d_coeff(n,m,b) result(d)

 		    integer, intent(in) :: n,m
 		    real(dp), intent(in) :: b
 		    real(dp) :: d

        d=(64._dp*sqrt(2._dp)*n*((n*b)**2 - (-1._dp+(m**2))))/((15._dp*pi)*((m**2)+(n*b)**2))

 	    end function d_coeff

      function p_coeff(n,m,b) result(p)

        integer, intent(in) :: n,m
        real(dp), intent(in) :: b
        real(dp) :: p

        p=(9*(((n-2._dp)*b)**2 - (m-2._dp)**2))/(2*((n*b)**2+(m**2)))

      end function p_coeff

	    function g_coeff(n,m,b,g) result(gamma_)

	    	integer, intent(in) :: n,m
		    real(dp), intent(in) :: g, b
	    	real(dp) :: gamma_

        gamma_=(4._dp*sqrt(2._dp)*n*b*g*m**3)/((-1+4*m**2)*pi*((m**2)+((n*b)**2)))

	    end function g_coeff

      function g_star_coeff(n,m,b,g) result(gamma_)

	    	integer, intent(in) :: n,m
		    real(dp), intent(in) :: g, b
	    	real(dp) :: gamma_

        gamma_=(4._dp*sqrt(2._dp)*n*b*g*m)/(pi*(4._dp*(m**2)-1._dp))

	    end function g_star_coeff

      function g_prime_coeff(n,m,b,g) result(gamma_)

        integer, intent(in) :: n,m
        real(dp), intent(in) :: g, b
        real(dp) :: gamma_

        gamma_=(3._dp*b*g)/(4._dp*(((n*b)**2)*(m**2)))

      end function g_prime_coeff

      function e_coeff(m,b) result(e)

        integer, intent(in) :: m
        real(dp),intent(in) :: b
        real(dp) :: e

        e=16._dp*sqrt(2._dp)*m/(5._dp*pi)
      end function e_coeff



end module coeffs10d
