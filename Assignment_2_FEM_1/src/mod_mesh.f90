module mod_mesh
    implicit none

contains

subroutine build_quad_mesh(nx, ny, x, y, nnode, nelem, conn)
    use iso_fortran_env, only: int32, real64
    implicit none

    integer(int32), intent(in) :: nx, ny
    real(real64), allocatable, intent(out) :: x(:), y(:)
    integer(int32), intent(out) :: nnode, nelem
    integer(int32), allocatable, intent(out) :: conn(:,:)

    integer(int32) :: i, j, n, e
    integer(int32) :: ll, lr, ur, ul

    if (nx < 2 .or. ny < 2) then
        error stop "build_quad_mesh: nx and ny must be >= 2"
    end if

    nnode = nx * ny
    allocate(x(nnode), y(nnode))

    ! node coordinates: i=1..nx, j=1..ny
    do j = 1, ny
        do i = 1, nx
            n = (j-1) * nx + i
            x(n) = real(i-1, real64) / real(nx-1, real64)
            y(n) = real(j-1, real64) / real(ny-1, real64)
        end do
    end do

    nelem = (nx - 1) * (ny - 1)
    allocate(conn(4, nelem))

    e = 0
    do j = 1, ny-1
        do i = 1, nx-1
            e = e + 1
            ll = (j-1) * nx + i        ! lower-left
            lr = ll + 1                ! lower-right
            ul = ll + nx               ! upper-left
            ur = ul + 1                ! upper-right

            conn(1,e) = ll
            conn(2,e) = lr
            conn(3,e) = ur
            conn(4,e) = ul
        end do
    end do
end subroutine build_quad_mesh

end module mod_mesh
