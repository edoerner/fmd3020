program parabola
    !
    ! Program that samples random numbers from the distribution f(x) = x^2 from 
    ! random numbers in [0,1).
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
implicit none

    integer, parameter :: nsamples = 50000

    real, dimension (nsamples) :: rnnos
    real, dimension (nsamples) :: samples
    real :: mean

    ! Initialize PRNG
    call random_seed

    ! We will use whole-array arithmetic. Get an array of random numbers
    call random_number(rnnos)

    ! Get the samples from the ffmc of the cuadratic function distribution
    samples = ffmc(rnnos)

    ! Calculate the mean value of the samples
    mean = sum(samples)/max(1,size(samples))
    write(6,'(A,F5.1)') 'Average of the samples : ', mean

    ! Save results to file
    open(unit=1, file='parabola.txt')
    write(unit=1, fmt='(F10.5)') samples
    close(1)

contains

    pure elemental real function ffmc(rnno)
        real, intent(in) :: rnno
        
        ffmc = 3.0*(rnno)**(1.0/3.0)
        
    end function ffmc

end program parabola
