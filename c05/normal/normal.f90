program normal
    !
    ! Program that samples random numbers from the distribution f(x) = e(-x^2/2) from 
    ! random numbers in [0,1) using the numerical inversion method.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
implicit none

    integer, parameter :: nsamples = 100000
    real, dimension (nsamples) :: samples

    integer, parameter :: nintervals = 30
    real, parameter :: a=0.0, b=3.0         ! boundaries of the domain

    real, dimension(nintervals+1) :: xi = 0.0
    real, dimension(nintervals) :: p_i = 0.0  ! average pdf of interval

    integer :: i, j
    real :: dx      ! size of interval
    real :: rnno1, rnno2

    ! First determine the xi's
    xi(1) = a

    xi_loop : do i = 1,nintervals
        ! First guess of the step
        dx = (b-a)/nintervals        

        xi_sampling: do
            
            xi(i+1) = xi(i) + dx
            p_i(i) = 1.0/(nintervals*(xi(i+1) - xi(i)))

            ! Here the if sentence depends on the monotonicity of the function.
            ! The normal distribution is always decreasing.
            if(p_i(i) <= norm_dist(xi(i)) .and. p_i(i) >= norm_dist(xi(i+1))) then
                exit    ! i.e. accept the results
            else
                ! We must check if we have underestimated or overestimated the average
                ! pdf of the interval to decide if we increase or decrease the step
                if(p_i(i) < norm_dist((xi(i+1)))) then
                    dx = 0.5*dx ! underestimation, decrease step size
                else 
                    dx = 1.5*dx ! overestimation, increase step size
                endif
                
            endif
            
        end do xi_sampling
        
    enddo xi_loop

    ! Set the last xi to the right boundary of the domain.
    xi(nintervals+1) = b

    ! Now start the sampling procedure from the normal distribution.

    sample_loop: do i = 1,nsamples
        ! sample the interval
        call random_number(rnno1)
        j = int(nintervals*rnno1) + 1

        ! With the interval sample the random variable
        call random_number(rnno2)
        samples(i) = xi(j) + rnno2*(xi(j+1) - xi(j))
    enddo sample_loop

    ! Save results to file
    open(unit=1, file='normal.txt')
    write(unit=1, fmt='(F10.5)') samples
    close(1)

contains

    real function norm_dist(x)
        real, intent(in) :: x
        real :: pi

        pi = 4.0*atan(1.0)

        norm_dist = sqrt(2.0/pi)*exp(-0.5*x**2)
        
    end function norm_dist

end program normal