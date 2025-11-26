module mod_solver
    use iso_fortran_env, only: real64
    use mod_types
    implicit none
contains

  subroutine cg_solver(A, b, x, tol, maxit, iters, relres)
    ! Conjugate Gradient solver: solves A x = b for symmetric positive-definite A
    real(kind=real64), intent(in) :: A(:,:)
    real(kind=real64), intent(in) :: b(:)
    real(kind=real64), intent(inout) :: x(:)
    real(kind=real64), intent(in) :: tol
    integer, intent(in) :: maxit
    integer, intent(out) :: iters
    real(kind=real64), intent(out) :: relres

    integer :: n, k
    real(kind=real64), allocatable :: r(:), p(:), Ap(:)
    real(kind=real64) :: alpha, beta, rsold, rsnew

    n = size(A,1)
    allocate(r(n), p(n), Ap(n))

    r = b - matmul(A, x)
    p = r
    rsold = dot_product(r, r)

    write(*,'(A50)') '--------------------------------------------------'
    write(*,'(A6,3X,A14)') '  Iter', '  sqrt(r^{T}*r)'
    write(*,'(A6,3X,A14)') '  ----', '  -------------'
    write(*,'(A50)') '--------------------------------------------------'

    do k = 1, maxit
      Ap = matmul(A, p)
      alpha = rsold / dot_product(p, Ap)
      x = x + alpha * p
      r = r - alpha * Ap
      rsnew = dot_product(r, r)
      if (sqrt(rsnew) <= tol) then
        iters = k
        relres = sqrt(rsnew)
        write(*,'(I6,3X,1PE14.8)') k, sqrt(rsnew)
        write(*,'(A50)') '--------------------------------------------------'
        return
      end if
      beta = rsnew / rsold
      p = r + beta * p
      rsold = rsnew
      write(*,'(I6,3X,1PE14.8)') k, sqrt(rsnew)
    end do
    iters = maxit
    relres = sqrt(rsnew)
  end subroutine cg_solver

end module mod_solver