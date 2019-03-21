program dices
    !
    ! Program that calculates the sum of two dices
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    implicit none

    integer, parameter :: nsamples = 10000
    integer, dimension(nsamples) :: sum

    integer :: i
    real :: rnno
    integer :: dice1, dice2

    ! Initialize PRNG
    call random_seed

    do i = 1,nsamples
        ! Sample first dice
        call random_number(rnno)
        dice1 = 1 + floor(rnno*6.0)

        ! Sample second dice
        call random_number(rnno)
        dice2 = 1 + floor(rnno*6.0)

        sum(i) = dice1 + dice2
    enddo

    ! Output results to file
    open(unit=1, file='sum.txt')
    do i = 1,nsamples
        write(unit=1, fmt='(I2)') sum(i)
    enddo

    close(1)

end program dices