module mod_mesh
    use iso_fortran_env, only: real64 
    use mod_types
    implicit none
contains

  subroutine build_structured_tri_mesh(nx, ny, x, y, nelem, conn)
    integer, intent(in) :: nx !< Number of nodes in x-direction
    integer, intent(in) :: ny !< Number of nodes in y-direction
    real(kind=real64), allocatable, intent(out) :: x(:) !< Allocated array of x-coordinates
    real(kind=real64), allocatable, intent(out) :: y(:) !< Allocated array of y-coordinates
    integer, intent(out) :: nelem !< Number of triangular elements
    integer, allocatable, intent(out) :: conn(:,:) !< Connectivity array (3 x nelem). Each column stores node indices of a triangle in counterclockwise-ish order.
    integer :: nnode, i, j, n, e, n1, n2, n3, n4

    if (nx < 2 .or. ny < 2) then
      write(*,*) "build_structured_tri_mesh: requires nx>=2 and ny>=2"
      stop 1
    end if

    nnode = nx * ny

    allocate(x(nnode), y(nnode))

    ! Fill node coordinates (row-major ordering: j=1..ny, i=1..nx)
    do j = 1, ny
      do i = 1, nx
        ! node index
        n = (j - 1) * nx + i
        x(n) = real(i - 1, kind=real64) / real(nx - 1, kind=real64)
        y(n) = real(j - 1, kind=real64) / real(ny - 1, kind=real64)
      end do
    end do

    ! Compute number of triangular elements and allocate connectivity
    nelem = 2 * (nx - 1) * (ny - 1)

    allocate(conn(3, nelem))

    ! Build connectivity. For each cell (i,j) with lower-left node n1
    e = 0
    do j = 1, ny - 1
      do i = 1, nx - 1
        n1 = (j - 1) * nx + i
        n2 = n1 + 1
        n3 = n1 + nx
        n4 = n3 + 1
        ! triangle 1: n1, n2, n4
        e = e + 1
        conn(:, e) = (/ n1, n2, n4 /)
        ! triangle 2: n1, n4, n3
        e = e + 1
        conn(:, e) = (/ n1, n4, n3 /)
      end do
    end do
  end subroutine build_structured_tri_mesh

end module mod_mesh
