module mod_rhs
    use iso_fortran_env, only: real64
    use mod_types
    implicit none
contains

  pure function rhs_function(x, y) result(r)
    real(kind=real64), intent(in) :: x, y
    real(kind=real64) :: r

    ! f = -Δu = 2*pi^2*sin(pi x) sin(pi y)
    r = 2.0_real64 * pi**2 * sin(pi * x) * sin(pi * y)
    !r = 1
  end function rhs_function
  
end module mod_rhs