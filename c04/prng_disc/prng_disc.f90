program prng_disc
    !
    ! Program that creates random numbers in [0,1), based on a Monte Carlo model of a 
    ! spinning disc.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    implicit none

    integer, parameter :: nsamples = 10000
    real, dimension(nsamples) :: samples

    integer :: i
    real :: rnno

    ! Initialize PRNG
    call prng_init

    do i = 1,nsamples
        call prng_get(rnno)
        samples(i) = rnno
    enddo

    ! Output result to file
    open(unit=1, file='prng_disc.txt')
    do i = 1,nsamples
        write(unit=1, fmt='(F5.3)') samples(i)
    enddo

    close(1)

contains

    subroutine prng_init

        call random_seed

    end subroutine prng_init

    subroutine prng_get(rnno)
        real,intent(out)  :: rnno
        real :: pi, phi

        pi = 4.0*atan(1.0)
        call random_number(rnno)

        phi = 2.0*pi*rnno
        rnno = phi/(2.0*pi)
        
    end subroutine prng_get

end program prng_disc