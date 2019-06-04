program total_sum
    !
    ! Program that initializes and sums the elements of an array. It paralelizes 
    ! the calculation using MPI.
    !
    ! Author : Edgardo Doerner, edoerner@fis.puc.cl
    !

    use mpi
    use iso_fortran_env, only: int32, int64, real64


implicit none

    integer(kind=int64), parameter :: nsize = 100000   ! size of the array
    integer(kind=int64), allocatable :: a(:)            ! array whose elements will be added
    integer(kind=int64) :: local_size

    real(kind=real64) :: start_time, end_time, rend_time

    real(kind=real64) :: local_res, total_res = 0
    integer(kind=int64) :: i, istart, iend, rest
    integer(kind=int32) :: mpi_rank, mpi_size, ierr

    call mpi_init(ierr)

    ! Identify MPI process and MPI size.
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, ierr)

    start_time = mpi_wtime()

    ! Print header with program information.
    if (mpi_rank .eq. 0) then
        write(6,'(A)') 'Function that calculates the addition of the elements of an array '
        write(6,'(A, I3)') 'Number of MPI processes : ', mpi_size
        write(6,'(A)') '*****************************************************************'
    endif

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    
    ! Get lower and upper array indices for each process.
    local_size = nsize/mpi_size
    
    istart = 1 + mpi_rank*local_size
    iend = (mpi_rank+1)*local_size
    
    ! Adjust if the array size and the number of mpi processes does not evenly 
    ! divide.
    rest = mod(nsize, mpi_size)
    if(mpi_rank .eq. mpi_size - 1) then
        iend = iend + rest
    endif

    write(6,'(A, I10)') 'Array : a(i) = i**2, i = ', istart, ', ', iend

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    ! Allocate and initialize array.
    a = [(i**2, i=istart,iend)]

    ! Perform addition using Fortran intrinsic sum function.
    local_res = 1.0_real64*sum(a)

    ! Print local result to console.
    write(6,'(A, F20.1)') 'Local sum of array elements = ', local_res

    ! Now reduce the calculation, output results.
    call mpi_reduce(local_res, total_res, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        0, MPI_COMM_WORLD, ierr)

    if (mpi_rank .eq. 0) then
        write(6,'(A)') '*******************************************************'
        write(6,'(A, F20.1)') 'Total sum of array elements = ', total_res
    endif

    ! Get elapsed time
    end_time = mpi_wtime()
    call mpi_reduce(end_time, rend_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
        0, MPI_COMM_WORLD, ierr)
    
    if(mpi_rank .eq. 0) then
        write(*,'(A,F15.5)') 'Elapsed time (s) : ', rend_time - start_time
    endif

    call mpi_finalize(ierr)

end program total_sum
