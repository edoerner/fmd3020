program normal_dist
    !
    ! Program that calculates the density and cumulative functions for a 
    ! normal distribution.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int32, real64

    implicit none

    integer(kind=int32), parameter :: nexp = 30 ! number of experiments
    real(kind=real64) :: mean = 10  ! mean of the normal distribution
    real(kind=real64) :: var = 10   ! variance of the normal distribution

    real(kind=real64), dimension(0:nexp) :: pdf = 0.0  ! pdf
    real(kind=real64), dimension(0:nexp) :: cdf = 0.0  ! cdf

    integer(kind=int32) :: i
    
    write(*,'(A,F10.2,A,F10.2)') 'Normal distribution for mean = ', mean, &
        ' and variance = ', var

    ! Calculate the density function for the normal distribution.
    do i = 0,nexp
        ! First calculate the pdf.
        pdf(i) = normal(1.0_8*i, mean, var)

        ! Now determine the cdf.
        cdf(i) = sum(pdf(0:i))
        
    enddo

    ! Save results to a file.
    open(unit=1, file='normal.txt')
    do i = 0,nexp
        write(unit=1, fmt='(I5, F10.5, F10.5)') i, pdf(i), cdf(i)        
    enddo
    close(1)

    ! Get the position of the maximum value of the pdf. 1 is substracted 
    ! because pdf index starts at zero.
    write(*,'(A,I5)') 'nmax = ', maxloc(pdf)-1

    contains

    real(kind=real64) function normal(x, mean, var)
        ! This function calculates the pdf of the normal distribution.
        real(kind=real64), intent(in) :: x
        real(kind=real64), intent(in) :: mean
        real(kind=real64), intent(in) :: var
        real(kind=real64), parameter :: pi = 4.0*atan(1.0) 

        normal = 1.0/sqrt(2*pi*var)*exp(-(x-mean)**2/(2.0*var))        
        
        
    end function normal

end program normal_dist