module mod_types
    use iso_fortran_env, only: real64
    implicit none

  real(kind = real64), parameter :: zero = 0.0_real64
  real(kind = real64), parameter :: pi = 4.0_real64 * datan(1.0_real64)
end module mod_types