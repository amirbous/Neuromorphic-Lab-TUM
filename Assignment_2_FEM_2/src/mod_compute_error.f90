module mod_compute_error
  use mod_types
  implicit none
contains

  function compute_L2_error(x, y, conn, nelem, u) result(L2)
    ! Compute L2 error by element centroid sampling
    ! x, y: nodal coordinates (1..nnode)
    ! conn: connectivity (3 x nelem), 1-based node indices for triangles
    ! nelem: number of elements
    ! u: computed nodal solution

    real(kind=real64), intent(in) :: x(:), y(:)
    integer, intent(in) :: conn(:, :)
    integer, intent(in) :: nelem
    real(kind=real64), intent(in) :: u(:)
    real(kind=real64) :: L2

    integer :: t, n1, n2, n3
    real(kind=real64) :: x1, x2, x3, y1, y2, y3, area, uc, ue, err2

    L2 = 0.0_real64
    do t = 1, nelem
      n1 = conn(1,t)
      n2 = conn(2,t)
      n3 = conn(3,t)

      x1 = x(n1)
      x2 = x(n2)
      x3 = x(n3)

      y1 = y(n1)
      y2 = y(n2)
      y3 = y(n3)

      area = 0.5_real64 * ( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) )
      if (area <= 0.0_real64) then
        area = abs(area)
      end if

      uc = (u(n1) + u(n2) + u(n3)) / 3.0_real64
      ue = sin(pi * ((x1 + x2 + x3) / 3.0_real64)) * sin(pi * ((y1 + y2 + y3) / 3.0_real64 ))
      err2 = (uc - ue)**2
      L2 = L2 + err2 * area
    end do

    L2 = sqrt(L2)
  end function compute_L2_error

end module mod_compute_error
