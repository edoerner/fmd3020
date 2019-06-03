program serial_sum
    !
    ! Program that initializes and sums the elements of an array.
    !
    ! Author : Edgardo Doerner, edoerner@fis.puc.cl
    !

    use iso_fortran_env, only: int64, real64


implicit none

    integer(kind=int64), parameter :: nsize = 100000 ! size of the array
    integer(kind=int64), dimension(nsize) :: a = 0  ! array whose elements will be added

    real(kind=real64) :: start_time, end_time
    real(kind=real64) :: res = 0.0
    integer(kind=int64) :: i

    call cpu_time(start_time)

    ! Print header with program information.
    write(6,'(A)') 'Function that calculates the addition of the elements of an array '
    write(6,'(A, I10)') 'Array : a(i) = i**2, i = 1, ', nsize
    write(6,'(A)') '*****************************************************************'

    ! Initialize array.
    a = [(i**2, i=1,nsize)]

    ! Perform addition using Fortran intrinsic sum function.
    res = 1.0_real64*sum(a)

    ! Print result to console.
    write(6,'(A, F20.1)') 'Sum of array elements = ', res

    ! Get elapsed time
    call cpu_time(end_time)
    write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time

end program serial_sum

