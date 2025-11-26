module mod_io
    use iso_fortran_env, only: real64, int32
    use mod_types
    implicit none
contains


  subroutine write_mesh_csv(nodes_filename, elems_filename, x, y, nnode, nelem, conn)
    use iso_fortran_env, only: int32, real64
    implicit none

    character(len=*), intent(in) :: nodes_filename, elems_filename
    real(real64), intent(in) :: x(:), y(:)
    integer(int32), intent(in) :: nnode, nelem
    integer(int32), intent(in) :: conn(:,:)
    integer :: i, e
    integer :: ios

    open(unit=10, file=nodes_filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR opening nodes file: ", trim(nodes_filename)
      return
    end if
    write(10,'(A)') 'idx,x,y'
    do i = 1, nnode
      write(10,'(I0,"," ,F0.12, "," ,F0.12)') i, x(i), y(i)
    end do
    close(10)

    open(unit=11, file=elems_filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR opening elems file: ", trim(elems_filename)
      return
    end if
    write(11,'(A)') 'idx,n1,n2,n3'
    do e = 1, nelem
      write(11,'(I0,"," ,I0,"," ,I0,"," ,I0)') e, conn(1,e), conn(2,e), conn(3,e)
    end do
    close(11)
  end subroutine write_mesh_csv


  subroutine write_solution_csv(nodes_file, elems_file, sol_file, x, y, conn, nnode, nelem, u)
    use iso_fortran_env, only: real64
    implicit none

    character(len=*), intent(in) :: nodes_file, elems_file, sol_file
    real(kind=real64), intent(in) :: x(:), y(:)
    integer, intent(in) :: conn(:, :)     !
    integer, intent(in) :: nnode, nelem
    real(kind=real64), intent(in) :: u(:)
    integer :: i, e, ios
    integer :: unit_nodes, unit_elems, unit_sol

    ! optionally rewrite nodes CSV for consistency (safe)
    unit_nodes = 10
    open(unit=unit_nodes, file=nodes_file, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: opening nodes CSV ', trim(nodes_file)
      return
    end if
    write(unit_nodes,'(A)') 'idx,x,y'
    do i = 1, nnode
      write(unit_nodes,'(I0,"," ,F0.12, "," ,F0.12)') i, x(i), y(i)
    end do
    close(unit_nodes)

    ! write elements CSV for TRIANGLES (3 node indices per element)
    unit_elems = 11
    open(unit=unit_elems, file=elems_file, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: opening elems CSV ', trim(elems_file)
      return
    end if
    write(unit_elems,'(A)') 'idx,n1,n2,n3'
    do e = 1, nelem
      write(unit_elems,'(I0,"," ,I0,"," ,I0,"," ,I0)') e, conn(1,e), conn(2,e), conn(3,e)
    end do
    close(unit_elems)

    ! write solution CSV
    unit_sol = 12
    open(unit=unit_sol, file=sol_file, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: opening solution CSV ', trim(sol_file)
      return
    end if
    write(unit_sol,'(A)') 'idx,u'
    do i = 1, nnode
      write(unit_sol,'(I0,"," ,F0.12)') i, u(i)
    end do
    close(unit_sol)

  end subroutine write_solution_csv


  !------------------------------------------------------------------------------
  !> @brief  Write matrix A to MatrixMarket (coordinate) file
  !> @param[in]  A(:,:)      real(kind=real64)   matrix to write
  !> @param[in]  tol         real(kind=real64)   threshold below which entries are treated as zero
  !> @param[in]  filename    character(len=*)    output filename
  !> @param[out] ioerr       integer             iostat (0 on success)
  !------------------------------------------------------------------------------
  subroutine write_A_matrix_market(A, filename, tol, ioerr)
    use iso_fortran_env, only: real64
    implicit none
    real(kind=real64), intent(in) :: A(:, :)
    character(len=*), intent(in) :: filename
    real(kind=real64), intent(in), optional :: tol
    integer, intent(out) :: ioerr

    integer :: nrow, ncol
    integer :: i, j, nnz
    integer :: unit
    real(kind=real64) :: thr

    if (present(tol)) then
      thr = tol
    else
      thr = 1.0e-14_real64
    end if

    nrow = size(A,1)
    ncol = size(A,2)

    ! count nonzeros above threshold
    nnz = 0
    do i = 1, nrow
      do j = 1, ncol
        if (abs(A(i,j)) > thr) nnz = nnz + 1
      end do
    end do

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ioerr)
    if (ioerr /= 0) return

    write(unit,'(A)') '%%MatrixMarket matrix coordinate real general'
    write(unit,'(A)') '% Matrix exported by mod_io::export_A_matrix_market (1-based indices)'
    write(unit,'(I0,1x,I0,1x,I0)') nrow, ncol, nnz

    do i = 1, nrow
      do j = 1, ncol
        if (abs(A(i,j)) > thr) then
          write(unit,'(I0,1x,I0,1x,1PE25.16)') i, j, A(i,j)
        end if
      end do
    end do

    close(unit)
    ioerr = 0
  end subroutine write_A_matrix_market

end module mod_io