module mod_assemble
    use iso_fortran_env, only: int32, real64
    implicit none

    abstract interface
        function src_func(x, y) result(val)
            import real64
            real(real64), intent(in) :: x, y
            real(real64) :: val
        end function src_func
    end interface

contains

subroutine assemble_poisson_quad(x, y, conn, nelem, A, f, src, apply_dirichlet)
    implicit none
    ! inputs
    real(real64), intent(in) :: x(:)
    real(real64), intent(in) :: y(:)
    integer(int32), intent(in) :: conn(:, :)
    integer(int32), intent(in) :: nelem
    ! outputs
    real(real64), allocatable, intent(out) :: A(:,:)
    real(real64), allocatable, intent(out) :: f(:)
    ! optional user source
    procedure(src_func), optional :: src
    logical, optional, intent(in) :: apply_dirichlet

    ! locals
    integer(int32) :: nnode
    integer :: e
    integer :: i
    integer :: k
    integer :: l
    integer :: m
    integer, dimension(4) :: nodes
    real(real64), dimension(4,4) :: ke
    real(real64), dimension(4) :: fe
    real(real64) :: xi_gp, eta_gp, w
    real(real64), dimension(2) :: gp, gw
    real(real64), dimension(4) :: N, dN_dxi, dN_deta
    real(real64), dimension(2,2) :: J, invJ
    real(real64) :: detJ
    real(real64), dimension(4) :: xloc, yloc, dN_dx, dN_dy
    real(real64) :: xg, yg, s_val

    nnode = size(x)
    if (size(y) /= nnode) then
        error stop "assemble_poisson_quad: x and y must have same size"
    end if

    if (size(conn,1) /= 4) then
        error stop "assemble_poisson_quad: conn must have 4 rows (quad elements)"
    end if

    ! allocate and zero global structures (dense)
    allocate(A(nnode, nnode))
    allocate(f(nnode))
    A = 0.0_real64
    f = 0.0_real64

    ! 2x2 Gauss points on reference (-1,1)
    gp = (/ -1.0_real64/sqrt(3.0_real64), 1.0_real64/sqrt(3.0_real64) /)
    gw = (/ 1.0_real64, 1.0_real64 /)

    ! Loop over elements
    do e = 1, nelem
        nodes = conn(:, e) ! 1-based node indices in ordering ll, lr, ur, ul
        ! local coordinates
        do i = 1, 4
            xloc(i) = x(nodes(i))
        yloc(i) = y(nodes(i))
        end do

        ke = 0.0_real64
        fe = 0.0_real64

        ! quadrature
        do k = 1, 2
            xi_gp = gp(k)
            do l = 1, 2
            eta_gp = gp(l)
            w = gw(k) * gw(l)

            ! shape functions (Q1) on reference [-1,1] x [-1,1]
            N(1) = 0.25_real64 * (1.0_real64 - xi_gp) * (1.0_real64 - eta_gp)  ! N1: lower-left
            N(2) = 0.25_real64 * (1.0_real64 + xi_gp) * (1.0_real64 - eta_gp)  ! N2: lower-right
            N(3) = 0.25_real64 * (1.0_real64 + xi_gp) * (1.0_real64 + eta_gp)  ! N3: upper-right
            N(4) = 0.25_real64 * (1.0_real64 - xi_gp) * (1.0_real64 + eta_gp)  ! N4: upper-left

            ! derivatives wrt xi,eta
            dN_dxi(1)  = -0.25_real64 * (1.0_real64 - eta_gp)
            dN_deta(1) = -0.25_real64 * (1.0_real64 - xi_gp)

            dN_dxi(2)  =  0.25_real64 * (1.0_real64 - eta_gp)
            dN_deta(2) = -0.25_real64 * (1.0_real64 + xi_gp)

            dN_dxi(3)  =  0.25_real64 * (1.0_real64 + eta_gp)
            dN_deta(3) =  0.25_real64 * (1.0_real64 + xi_gp)

            dN_dxi(4)  = -0.25_real64 * (1.0_real64 + eta_gp)
            dN_deta(4) =  0.25_real64 * (1.0_real64 - xi_gp)

            ! Jacobian J = [ dx/dxi  dx/deta ; dy/dxi dy/deta ]
            J = 0.0_real64
            do i = 1, 4
                J(1,1) = J(1,1) + xloc(i) * dN_dxi(i)
                J(1,2) = J(1,2) + xloc(i) * dN_deta(i)
                J(2,1) = J(2,1) + yloc(i) * dN_dxi(i)
                J(2,2) = J(2,2) + yloc(i) * dN_deta(i)
            end do

            detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1)
            if (abs(detJ) <= 0.0_real64) then
                error stop "assemble_poisson_quad: singular element Jacobian (check element ordering)."
            end if

            invJ(1,1) =  J(2,2) / detJ
            invJ(1,2) = -J(1,2) / detJ
            invJ(2,1) = -J(2,1) / detJ
            invJ(2,2) =  J(1,1) / detJ

            ! compute physical gradients of shape functions
            do i = 1, 4
                dN_dx(i) = invJ(1,1) * dN_dxi(i) + invJ(1,2) * dN_deta(i)
                dN_dy(i) = invJ(2,1) * dN_dxi(i) + invJ(2,2) * dN_deta(i)
            end do

            ! evaluate source function at quadrature point (if provided)
            if (present(src)) then
                xg = 0.0_real64
                yg = 0.0_real64
                do i = 1, 4
                    xg = xg + xloc(i) * N(i)
                    yg = yg + yloc(i) * N(i)
                end do
                s_val = src(xg, yg)
            else
                s_val = 0.0_real64
            end if

            ! accumulate local stiffness and load
            do i = 1, 4
                do m = 1, 4
                    ke(i,m) = ke(i,m) + ( dN_dx(i) * dN_dx(m) + dN_dy(i) * dN_dy(m) ) * detJ * w
                end do
                fe(i) = fe(i) + N(i) * s_val * detJ * w
            end do

        end do
    end do

    ! assemble into global A and f
    do i = 1, 4
        do m = 1, 4
            A(nodes(i), nodes(m)) = A(nodes(i), nodes(m)) + ke(i,m)
        end do
        f(nodes(i)) = f(nodes(i)) + fe(i)
    end do

    end do  ! elements

    ! apply homogeneous Dirichlet BCs on boundary nodes if requested
    if (present(apply_dirichlet)) then
        if (apply_dirichlet) then
            do i = 1, nnode
                if ( x(i) <= 0.0_real64 + 1.0e-12_real64 .or. x(i) >= 1.0_real64 - 1.0e-12_real64 &
                .or. y(i) <= 0.0_real64 + 1.0e-12_real64 .or. y(i) >= 1.0_real64 - 1.0e-12_real64 ) then
                ! zero row and column
                    A(:, i) = 0.0_real64
                    A(i, :) = 0.0_real64
                    A(i, i) = 1.0_real64
                    f(i) = 0.0_real64
                end if
            end do
        end if
    end if

end subroutine assemble_poisson_quad

end module mod_assemble
