module mod_rhs
    use iso_fortran_env, only: real64
    implicit none

contains

function src_const(xp, yp) result(val)
    real(real64), intent(in) :: xp, yp
    real(real64) :: val
    real(real64), parameter :: pi = 4.0_real64 * datan(1.0_real64)

    val = sqrt(xp * yp**2) * sin(xp * pi) * cos(-yp * pi)
end function src_const
end module mod_rhs
