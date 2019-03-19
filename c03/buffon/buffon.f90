program buffon
    !
    ! Program that calculates pi using the Monte Carlo method through 
    ! Buffon's needle problem
    !

    implicit none

    real, parameter :: l = 0.5          ! needle's length (cm)
    real, parameter :: d = 1.0          ! distance between vertical lines (cm)
    integer, parameter :: nl = 100      ! number of throwings

    real :: x               ! x coordinate of needle's center of mass (cm)
    real :: theta           ! orientation of the needle (rad) within (0, pi)
    real :: pi              ! internal value of pi
    real :: pi_buffon       ! estimated value of pi
    real :: xl, xr          ! position of needle's tips (cm)
    real :: rnno            ! random number within [0,1)

    integer :: ni = 0       ! number of intersections
    integer :: i            ! loop index

    pi = 4.0*atan(1.0)
    write(6,*) 'Fortran internal value of pi : ', pi

    ! Initialize PRNG
    call random_seed

    do i = 1,nl
        ! needle's cm coordinate
        call random_number(rnno)
        x = d*rnno

        ! needle's orientation
        call random_number(rnno)
        theta = pi*rnno-pi/2.0

        ! needle's tips coordinates
        xl = x - (l/2.0)*cos(theta)
        xr = x + (l/2.0)*cos(theta)

        ! check if needle crossed a vertical line
        if ((xl .le. 0.0) .or. (xr .ge. d)) then
            ni = ni + 1
        end if
    enddo

    ! Finally, calculate value of pi based of the probability that a needle crosses 
    ! a vertical line.
    pi_buffon = (real(nl)/real(ni))*(2*l/d)

    write(6,*) 'Approximated value of pi : ', pi_buffon
    write(6,*) 'Number of throwings : ', nl
    write(6,*) 'Number of intersections : ', ni

end program buffon