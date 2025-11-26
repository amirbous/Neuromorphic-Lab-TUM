!------------------------------------------------------------------------------
! High Performance Computing Center Stuttgart (HLRS)
!------------------------------------------------------------------------------
!
!> @author
!> Benjamin Schnabel
!
!> @brief
!> Simple FEM application.
!
! REVISION HISTORY:
! 01.05.2024 - Initial Version
! 18.11.2025 - Move outputs out of src/ folder to data/ folder and including subroutine write_A_matrix_market. (F. Ruhnke)
!------------------------------------------------------------------------------

program main
    use mod_mesh
    use mod_assemble
    use mod_rhs
    use mod_solver
    use mod_io

    use iso_fortran_env, only: int32, real64

    implicit none

    integer(int32) :: nx, ny, nnode, nelem
    real(real64), allocatable :: x(:), y(:)
    real(real64), allocatable :: A(:,:), f(:)
    integer(int32), allocatable :: conn(:,:)
    real(kind=real64), allocatable :: u(:)
    real(real64) :: tol, relres
    integer :: maxit, iters
    integer :: ioerr

    character(len=256) :: matrix_dir, mesh_dir, solution_dir
    character(len=256) :: nodes_filename, elements_filename
    character(len=256) :: vtk_mesh_filename, vtk_solution_filename
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
    vtk_mesh_filename = trim(mesh_dir) // '/mesh.vtk'
    vtk_solution_filename = trim(solution_dir) // '/solution.vtk'
    solution_csv_filename = trim(solution_dir) // '/solution.csv'



    write(*,"(A50)") "##################################################"
    write(*,"(A50)") "# Building mesh                                  #"
    write(*,"(A50)") "##################################################"

    call build_quad_mesh(nx, ny, x, y, nnode, nelem, conn)
    nnode = nx * ny

    write(*,"(A6,I4)") "  nx: ", nx
    write(*,"(A6,I4)") "  ny: ", ny
    write(*,"(A9,I4)") "  Nodes: ", nnode
    write(*,"(A12,I4)") "  Elements: ", nelem



    write(*,"(A50)") "##################################################"
    write(*,"(A50)") "# Assemble system                                #"
    write(*,"(A50)") "##################################################"

    call assemble_poisson_quad(x, y, conn, nelem, A, f, src=src_const, apply_dirichlet=.true.)




    write(*,"(A50)") "##################################################"
    write(*,"(A50)") "# Solve system with CG method                    #"
    write(*,"(A50)") "##################################################"

    ! Allocate solution vector and set initial guess to zero before calling solver
    allocate(u(nx * ny))
    u = 0.0_real64
    call cg_solver(A, f, u, tol, maxit, iters, relres)

    write(*,"(A13,E10.3)") "  Tolerance: ", tol
    write(*,"(A18,I4)") "  Max iterations: ", nnode
    write(*,"(A31,I4)") "  CG finished after iteration: ", iters
    write(*,"(A10,E10.3)") "  relres: ", relres



    write(*,"(A50)") "##################################################"
    write(*,"(A50)") "# Exporting results                              #"
    write(*,"(A50)") "##################################################"

    call write_mesh_csv(nodes_filename, elements_filename, x, y, nnode, nelem, conn)
    call write_mesh_vtk(vtk_mesh_filename, x, y, nnode, nelem, conn)
    call write_A_matrix_market(A, A_matrix_filename, ioerr=ioerr)
    call write_solution_csv(nodes_filename, elements_filename, solution_csv_filename, x, y, conn, nnode, nelem, u)
    call write_solution_vtk(vtk_solution_filename, x, y, conn, nnode, nelem, u)

end program main
