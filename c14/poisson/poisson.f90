program poisson_dist
    !
    ! Program that calculates the density and cumulative functions for a 
    ! poisson distribution.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int32, real64

    implicit none

    integer(kind=int32), parameter :: nexp = 20 ! number of experiments
    real(kind=real64) :: mean = 10 ! mean value for the Poisson distribution

    real(kind=real64), dimension(0:nexp) :: pdf = 0.0  ! pdf
    real(kind=real64), dimension(0:nexp) :: cdf = 0.0  ! cdf

    integer(kind=int32) :: i
    
    write(*,'(A,F10.2,A,I5)') 'Poisson distribution for m = ', mean, &
        ' and N = ', nexp

    ! Calculate the functions values for zero success.
    pdf(0) = poisson(0, mean)
    cdf(0) = pdf(0)

    ! Now apply a recursive formulation to calculate the distribution and its
    ! cumulative function.
    do i = 1,nexp
        ! First calculate the pdf.
        pdf(i) = poisson(i, mean)

        ! Now determine the cdf.
        cdf(i) = sum(pdf(0:i))
        
    enddo

    ! Save results to a file.
    open(unit=1, file='poisson.txt')
    do i = 0,nexp
        write(unit=1, fmt='(I5, F10.5, F10.5)') i, pdf(i), cdf(i)        
    enddo
    close(1)

    ! Get the position of the maximum value of the pdf. 1 is substracted 
    ! because pdf index starts at zero.
    write(*,'(A,I5)') 'nmax = ', maxloc(pdf)-1

    contains

    real(kind=real64) function poisson(n, m)
        ! This function calculates the pdf of the poisson distribution, given 
        ! a number n of trials and the mean value m.
        integer(kind=int32), intent(in) :: n
        real(kind=real64), intent(in) :: m

        integer(kind=int32) :: itr
        
        ! Initialize recursive formulation.
        poisson = exp(-m)

        ! Calculate the pdf for the desired number of experiments.
        do itr = 1,n
            poisson = (m/itr)*poisson            
        enddo
        
    end function poisson

end program poisson_dist