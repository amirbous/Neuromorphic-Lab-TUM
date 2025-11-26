module mod_assemble
    use iso_fortran_env, only: real64
    use mod_types
    use mod_rhs
    implicit none
contains

  subroutine assemble_poisson_tri(x, y, conn, nelem, A, b, is_bc)

    real(kind=real64), intent(in) :: x(:) !< x-coordinates of nodes
    real(kind=real64), intent(in) :: y(:) !< y-coordinates of nodes
    integer, intent(in) :: conn(:, :) !< Connectivity array
    integer, intent(in) :: nelem !< Number of triangular elements
    real(kind=real64), allocatable, intent(out) :: A(:,:) !< Global stiffness matrix
    real(kind=real64), allocatable, intent(out) :: b(:) !< Global load vector

    integer :: nnode
    ! is_bc is now an output argument (allocated inside this subroutine)
    logical, allocatable, intent(out) :: is_bc(:)

    integer :: t, i, j
    integer :: n1, n2, n3
    real(kind=real64) :: x1, x2, x3, y1, y2, y3
    real(kind=real64) :: area, twoArea
    real(kind=real64) :: bi(3), ci(3)
    real(kind=real64) :: ke(3,3), fe(3)
    real(kind=real64) :: cx, cy, fval

    ! Allocate and initialize global matrix and vector
    nnode = size(x)
    if (.not. allocated(A)) allocate(A(nnode, nnode))
    if (.not. allocated(b)) allocate(b(nnode))
    A = zero
    b = zero

    ! Loop over all elements
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

      ! triangle area
      area = 0.5_real64 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
      if (area <= zero) then
        ! Degenerate element; skip
        cycle
      end if

      twoArea = 2.0_real64 * area

      ! Gradients of linear shape functions (bi, ci)
      bi(1) = (y2 - y3) / twoArea
      bi(2) = (y3 - y1) / twoArea
      bi(3) = (y1 - y2) / twoArea
      ci(1) = (x3 - x2) / twoArea
      ci(2) = (x1 - x3) / twoArea
      ci(3) = (x2 - x1) / twoArea

      ! Element stiffness ke_ij = (b_i*b_j + c_i*c_j)*area
      ke = zero
      do i = 1,3
        do j = 1,3
          ke(i,j) = (bi(i) * bi(j) + ci(i) * ci(j)) * area
        end do
      end do

      ! Element load: use centroid one-point quadrature: fe_i = f(centroid)*area/3
      cx = (x1 + x2 + x3) / 3.0_real64
      cy = (y1 + y2 + y3) / 3.0_real64
      fval = rhs_function(cx, cy)
      fe = fval * area / 3.0_real64

      ! Assemble element stiffness into global matrix
      do i = 1,3
        do j = 1,3
          A(conn(i,t), conn(j,t)) = A(conn(i,t), conn(j,t)) + ke(i,j)
        end do
      end do

      ! Assemble load vector
      b(conn(1,t)) = b(conn(1,t)) + fe(1)
      b(conn(2,t)) = b(conn(2,t)) + fe(2)
      b(conn(3,t)) = b(conn(3,t)) + fe(3)
    end do

    ! ------------------------------------------------------------------
    ! Dirichlet BC: u = 0 on boundary (nodes with x=0, x=1, y=0, y=1)
    ! Compute boundary marker array and apply BCs by modifying A and b
    allocate(is_bc(nnode))
    is_bc = .false.
    do i = 1, nnode
      if ( x(i) <= 1.0e-12_real64 .or. x(i) >= 1.0_real64-1.0e-12_real64 &
        .or. y(i) <= 1.0e-12_real64 .or. y(i) >= 1.0_real64-1.0e-12_real64 ) then
        is_bc(i) = .true.
      end if
    end do

    ! Apply Dirichlet BCs by modifying the system (set rows/cols to zero,
    ! diagonal to 1 and RHS to zero) so that u(boundary) = 0 is enforced.
    do i = 1, nnode
      if (is_bc(i)) then
        A(i, :) = zero
        A(:, i) = zero
        A(i, i) = 1.0_real64
        b(i) = zero
      end if
    end do

  end subroutine assemble_poisson_tri
    
end module mod_assemble