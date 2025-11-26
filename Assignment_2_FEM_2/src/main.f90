!------------------------------------------------------------------------------
! High Performance Computing Center Stuttgart (HLRS)
!------------------------------------------------------------------------------
! MODULE: types
!
!> @author
!> Benjamin Schnabel
!
!> @brief:
!> Solves \f$ - \Delta u = f \f$ on unit square with \f$ u = 0 \f$ on the
!> boundary using linear triangular FEM and a CG solver.
!
! REVISION HISTORY:
! 24.10.2025 - Initial Version
! 24.11.2025 - Output files relocated to 'data/'. Code refactored into separate Fortran modules  and added subroutines for writing mesh and solution CSV files. (F.Ruhnke)
!------------------------------------------------------------------------------


program main

  use iso_fortran_env, only: real64
  use mod_mesh
  use mod_assemble
  use mod_rhs
  use mod_solver
  use mod_io
  use mod_types
  use mod_compute_error
  implicit none

  integer :: nx, ny, nnode, nelem 
  real(kind=real64), allocatable :: x(:), y(:)
  integer, allocatable :: conn(:,:)
  real(kind=real64), allocatable :: b(:)
  integer :: i, j, t
  integer :: ioerr
  logical, allocatable :: is_bc(:)
  integer :: nfree
  real(kind=real64), allocatable :: A(:,:)
  real(kind=real64), allocatable :: u(:)
  integer :: maxit, iters
  real(kind=real64) :: tol, relres
  real(kind=real64) :: L2err

  character(len=256) :: matrix_dir, mesh_dir, solution_dir
  character(len=256) :: nodes_filename, elements_filename
  character(len=256) :: solution_csv_filename
  character(len=256) :: A_matrix_filename


  ! -------------------------------------------------------------------------
  ! USER-CONFIGURABLE PARAMETERS
  !
  ! Edit the values below to change the experiment configuration. These are
  ! simple in-source defaults intended for students to modify before running.
  !
  !  - nx    : number of nodes in x-direction (integer)
  !  - ny    : number of nodes in y-direction (integer)
  !  - tol   : solver tolerance for the CG solver (real, e.g. 1.0e-10)
  !  - maxit : maximum number of CG iterations (integer)
  ! -------------------------------------------------------------------------
  nx = 30                 
  ny = 30                 
  tol = 1.0e-10_real64    
  maxit = nx * ny              ! default: nx * ny (scaled with problem size)
  ! -------------------------------------------------------------------------

  ! Create output directories
  matrix_dir = 'data/matrix'
  mesh_dir = 'data/mesh'
  solution_dir = 'data/solution'
  call execute_command_line('mkdir -p ' // trim(matrix_dir))      ! Create matrix output directory (if not existing)
  call execute_command_line('mkdir -p ' // trim(mesh_dir))        ! Create mesh output directory (if not existing)
  call execute_command_line('mkdir -p ' // trim(solution_dir))    ! Create solution output directory (if not existing)
  A_matrix_filename = trim(matrix_dir) // '/A.mtx'
  nodes_filename = trim(mesh_dir) // '/nodes.csv'
  elements_filename = trim(mesh_dir) // '/elements.csv'
  solution_csv_filename = trim(solution_dir) // '/solution.csv'


  write(*,"(A50)") "##################################################"
  write(*,"(A50)") "# Building mesh                                  #"
  write(*,"(A50)") "##################################################"

  call build_structured_tri_mesh(nx, ny, x, y, nelem, conn)
  nnode = nx * ny

  write(*,"(A6,I4)") "  nx: ", nx
  write(*,"(A6,I4)") "  ny: ", ny
  write(*,"(A9,I4)") "  Nodes: ", nnode
  write(*,"(A12,I4)") "  Elements: ", nelem


  write(*,"(A50)") "##################################################"
  write(*,"(A50)") "# Assemble system                                #"
  write(*,"(A50)") "##################################################"

  ! Assemble system (assemble_poisson_tri now also computes/apply Dirichlet BCs
  ! and returns the boolean marker array `is_bc`).
  call assemble_poisson_tri(x, y, conn, nelem, A, b, is_bc)

  ! Compute number of free DOFs for info
  nfree = 0
  do i = 1, nnode
    if (.not. is_bc(i)) then
      nfree = nfree + 1
    end if
  end do

  ! Initial guess 0
  allocate(u(nnode))
  u = zero

  write(*,"(A13,I4)") "  Free DOFs: ", nfree
  write(*,"(A18,I4)") "  Dirichlet DOFs: ", nnode - nfree


  write(*,"(A50)") "##################################################"
  write(*,"(A50)") "# Solve system with CG method                    #"
  write(*,"(A50)") "##################################################"

  call cg_solver(A, b, u, tol, maxit, iters, relres)

  write(*,"(A13,E10.3)") "  Tolerance: ", tol
  write(*,"(A18,I4)") "  Max iterations: ", nnode
  write(*,"(A31,I4)") "  CG finished after iteration: ", iters
  write(*,"(A10,E10.3)") "  relres: ", relres


  write(*,"(A50)") "##################################################"
  write(*,"(A50)") "# Exporting results                              #"
  write(*,"(A50)") "##################################################"

  call write_mesh_csv(nodes_filename, elements_filename, x, y, nnode, nelem, conn)
  call write_A_matrix_market(A, A_matrix_filename, ioerr=ioerr)
  call write_solution_csv(nodes_filename, elements_filename, solution_csv_filename, x, y, conn, nnode, nelem, u)


  write(*,"(A50)") "##################################################"
  write(*,"(A50)") "# Compare analytical and numerical result        #"
  write(*,"(A50)") "##################################################"

  ! compute L2 error using element centroid sampling
  L2err = compute_L2_error(x, y, conn, nelem, u)
  write(*,"(A12,E10.3)") "  L2 error: ", L2err


end program main
