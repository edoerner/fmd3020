program dice
    !
    ! Program that generates uniform integer numbers between 1 and 6, simulating
    ! a dice throwing.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    implicit none

    integer, parameter :: nsamples = 1000000
    integer, dimension(nsamples) :: n          ! dice throwings

    integer :: i
    real :: rnno

    ! Use intrinsic procedures random_seed and random_get.
    ! random_get generates real numbers within [0,1)
    call random_seed

    do i = 1,nsamples
        call random_number(rnno)
        n(i) = 1 + floor(6.0*rnno)
    enddo

    open(unit=1, file='dice.txt')
    do i = 1,nsamples
        write(unit=1, fmt='(I1)') n(i)
    enddo

    close(1)

end program dice