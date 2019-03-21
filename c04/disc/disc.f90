program disc
    !
    ! Program that creates random angles, based on a Monte Carlo model of a 
    ! spinning disc.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !

    implicit none

    integer, parameter :: nsamples = 10000
    real, dimension(nsamples) :: samples

    real :: pi
    real :: rnno
    integer :: i

    pi = 4.0*atan(1.0)
    write(6,*) 'Fortran internal value of pi : ', pi

    ! Initialize PRNG
    call random_seed

    do i = 1,nsamples
        call random_number(rnno)
        samples(i) = rnno*2.0*pi
        
    enddo

    ! Output result to file
    open(unit=1, file='disc.txt')
    do i = 1,nsamples
        write(unit=1, fmt='(F5.2)') samples(i)
    enddo

    close(1)

end program disc