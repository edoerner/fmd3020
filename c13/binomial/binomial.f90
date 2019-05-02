program binomial
    !
    ! Program that calculates the density and cumulative functions for a 
    ! binomial distribution, given the probability of success and a 
    ! number of experiments.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int32, real64

    implicit none

    integer(kind=int32), parameter :: nexp = 100 ! number of experiments
    real(kind=real64), parameter :: p = 0.7      ! probabily of success

    real(kind=real64), dimension(0:nexp) :: pdf = 0.0  ! pdf
    real(kind=real64), dimension(0:nexp) :: cdf = 0.0  ! cdf

    integer(kind=int32) :: i

    write(*,'(A,F10.5,A,I5)') 'Binomial distribution for p = ', p, &
        ' and N = ', nexp

    ! Calculate the functions values for zero success.
    pdf(0) = (1.0 - p)**nexp
    cdf(0) = pdf(0)

    ! Now apply a recursive formulation to calculate the distribution and its
    ! cumulative function.
    do i = 1,nexp
        ! First calculate the pdf.
        pdf(i) = (p/(1.0 - p))*((nexp - i + 1.0)/i)*pdf(i-1)

        ! Now determine the cdf.
        cdf(i) = sum(pdf(0:i))
        
    enddo

    ! Save results to a file.
    open(unit=1, file='binomial.txt')
    do i = 0,nexp
        write(unit=1, fmt='(I5, F10.5, F10.5)') i, pdf(i), cdf(i)        
    enddo
    close(1)

    ! Get the position of the maximum value of the pdf. 1 is substracted 
    ! because pdf index starts at zero.
    write(*,'(A,I5)') 'nmax = ', maxloc(pdf)-1

end program binomial