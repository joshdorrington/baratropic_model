!Joshua Dorrington 22/03/18 - University of Oxford
!This module generates the coefficients of the barotropic model,
!and also generates the linear time evolution operator

module coeffs
    use, intrinsic :: iso_fortran_env
    use params
    implicit none
    private
    public get_coeffs, build_lin_op


    contains

	    subroutine get_coeffs(b,beta,g,coeff_arr)
      		real(dp) :: coeff_arr(coeff_num)
      		real(dp), intent(in) :: b, beta, g
          ! the alpha advection coefficients
      		coeff_arr(1)=a_coeff(1,b)
      		coeff_arr(2)=a_coeff(2,b)
          coeff_arr(3)=a_coeff(3,b)
      		coeff_arr(4)=a_tilde_coeff(1,b)
          coeff_arr(5)=a_tilde_coeff(2,b)
      		coeff_arr(6)=a_tilde_coeff(3,b)
          !the beta coriolis coefficients
      		coeff_arr(7)=beta_coeff(1,b,beta)
      		coeff_arr(8)=beta_coeff(2,b,beta)
          coeff_arr(9)=beta_coeff(3,b,beta)
          !the gamma orographic coefficients
      		coeff_arr(10)=gamma_coeff(1,b,g)
      		coeff_arr(11)=gamma_coeff(2,b,g)
          coeff_arr(12)=gamma_coeff(3,b,g)
          coeff_arr(13)=gamma_tilde_coeff(1,b,g)
          coeff_arr(14)=gamma_tilde_coeff(2,b,g)
          coeff_arr(15)=gamma_tilde_coeff(3,b,g)
      		coeff_arr(16)=gamma_dagger_coeff(1,b,g)
      		coeff_arr(17)=gamma_dagger_coeff(3,b,g)
          coeff_arr(18)=gamma_star_coeff(1,b,g)
          coeff_arr(19)=gamma_star_coeff(3,b,g)
          !the delta advection coefficients
      		coeff_arr(20)=d_coeff(1,b)
      		coeff_arr(21)=d_coeff(2,b)
          coeff_arr(22)=d_tilde_coeff(2,b)
          coeff_arr(23)=d_tilde_coeff(3,b)
          !the kappa advection coefficients
          coeff_arr(24)=k_coeff(1,b)
          coeff_arr(25)=k_coeff(3,b)
          coeff_arr(26)=k_tilde_coeff(1,b)
          coeff_arr(27)=k_tilde_coeff(3,b)
          !the epsilon coeffincients
          !epsilon
      		coeff_arr(28)=16._dp*sqrt(2._dp)/(5._dp*pi)
          !epsilon_tilde
          coeff_arr(29)=64._dp*sqrt(2._dp)/(15._dp*pi)
          !epsilon_dagger
          coeff_arr(30)=80._dp*sqrt(2._dp)/(21._dp*pi)
          !epsilon*
          coeff_arr(31)=192._dp*sqrt(2._dp)/(35._dp*pi)
	    end subroutine get_coeffs

      subroutine build_lin_op(mat,coeff)
	        real(dp), dimension(coeff_num),intent(in) :: coeff
      		real(dp), dimension(dims,dims) :: mat
      		integer :: i, j

      		mat=0._dp

      		do i=1, dims
      		  mat(i,i)=-C
      		end do
      		mat(1,3)=coeff(13)!gamma_tilde_1
          mat(1,9)=-coeff(18)!-gamma_star_1

      		mat(2,3)=coeff(7)!beta_1

          mat(3,1)=-coeff(10)!-gamma_1
      		mat(3,2)=-coeff(7)!-beta_1
          mat(3,7)=coeff(16)!gamma_dagger_1

      		mat(4,6)=coeff(14)!gamma_tilde_2

      		mat(5,6)=coeff(8)!beta_2

      		mat(6,4)=-coeff(11)!-gamma_2
      		mat(6,5)=-coeff(8)!-beta_2

          mat(7,3)=-coeff(19)!-gamma_star_3
          mat(7,9)=coeff(15)!gamma_tilde_3

          mat(8,9)=coeff(9)!beta_3

          mat(9,7)=-coeff(9)!-beta_3
          mat(9,1)=coeff(17)!gamma_dagger_3
          mat(9,4)=-coeff(12)!-gamma_3
	    end subroutine build_lin_op
!THE ALPHA FUNCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    function a_coeff(m,b) result(a)
	    	integer, intent(in) :: m
		real(dp), intent(in) :: b
	    	real(dp) :: a
	    	a=8._dp*sqrt(2._dp)*(m**2)*(b**2+(m**2)-1)/(pi*(-1 +4*(m**2))*(b**2+m**2))
	    end function a_coeff

      function a_tilde_coeff(m,b) result(a)
        integer, intent(in) :: m
    real(dp), intent(in) :: b
        real(dp) :: a
        a=8._dp*sqrt(2._dp)*(m**2)*(b**2+(m**2)-9)/(pi*Abs((-9 +4*(m**2)))*(b**2+m**2))
      end function a_tilde_coeff
!THE BETA FUNC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    function beta_coeff(m,b,beta) result(output)
	    	integer, intent(in) :: m
		real(dp), intent(in) :: beta, b
	    	real(dp) :: output
	    	output=beta*(b**2)/(b**2 + m**2)
	    end function beta_coeff
!THE GAMMA FUNCs !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    function gamma_coeff(m,b,g) result(gamma_)
	    	integer, intent(in) :: m
		    real(dp), intent(in) :: g, b
	    	real(dp) :: gamma_
	    	gamma_=(g*4._dp*sqrt(2._dp)*b*m**3)/((4*(m**2)-1)*pi*(b**2+m**2))
	    end function gamma_coeff

	    function gamma_tilde_coeff(m,b,g) result(gamma_p)
		    integer, intent(in) :: m
		    real(dp), intent(in) :: g,b
		    real(dp) :: gamma_p
		    gamma_p=(g*4._dp*sqrt(2._dp)*b*m)/((4*(m**2)-1)*pi)
      end function gamma_tilde_coeff

     function gamma_dagger_coeff(m,b,g) result(gamma_)
       integer, intent(in) :: m
       real(dp), intent(in) :: g, b
       real(dp) :: gamma_
       gamma_=(g*4._dp*m*sqrt(2._dp)*b)/(15*pi*(b**2+m**2))
     end function gamma_dagger_coeff

     function gamma_star_coeff(m,b,g) result(gamma_)
       integer, intent(in) :: m
       real(dp), intent(in) :: g, b
       real(dp) :: gamma_
       gamma_=(g*4._dp*sqrt(2._dp)*b)/(15*pi*m)
     end function gamma_star_coeff

!THE DELTA FUNCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	   function d_coeff(m,b) result(d)
		   integer, intent(in) :: m
		   real(dp), intent(in) :: b
		   real(dp) :: d
		   d=(64._dp*sqrt(2._dp)*(b**2-(m**2) + 1))/(15._dp*pi*(b**2 + m**2))
	   end function d_coeff

     function d_tilde_coeff(m,b) result(d)
       integer, intent(in) :: m
       real(dp), intent(in) :: b
       real(dp) :: d
       d=(64._dp*sqrt(2._dp)*(b**2-(m**2) + 9))/(21._dp*pi*(b**2 + m**2))
     end function d_tilde_coeff

!THE KAPPA FUNCS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function k_coeff(m,b) result(k)
      integer, intent(in) :: m
      real(dp), intent(in) :: b
      real(dp) :: k
      k=(8._dp*sqrt(2._dp)*(b**2-(m**2) + 9))/(15._dp*pi*(b**2 + m**2))
    end function k_coeff

    function k_tilde_coeff(m,b) result(k)
      integer, intent(in) :: m
      real(dp), intent(in) :: b
      real(dp) :: k
      k=(216._dp*sqrt(2._dp)*(b**2-(m**2) + 1))/(35._dp*pi*(b**2 + m**2))
    end function k_tilde_coeff

end module coeffs
