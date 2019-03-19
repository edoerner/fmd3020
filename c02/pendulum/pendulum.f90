program pendulum
!
! This program solves the movement equation for a simple pendulum using the 
! explicit Euler method.
!
    implicit none

    integer, parameter :: nm = 1000 ! number of time steps
    real, parameter :: dt = 0.01    ! time step (s)
    real, parameter :: g = 9.806    ! standard gravitational acceleration (m/s^2)
    real, parameter :: l = 1.0      ! legth of the pendulum (m)

    ! Initial values for angle theta and angular velocity omega
    real, parameter :: thetain = 0.79   ! aproximately pi/4
    real, parameter :: omegain = 0.0

    integer :: n    ! indices in time

    real, dimension(nm) :: theta, omega

    theta = 0.0
    omega = 0.0

    ! Apply initial values
    theta(1) = thetain
    omega(1) = omegain

    print *, 'Solution of the movement equation for a simple pendulum using ' 
    print *, 'an Euler explicit approach.'
    print *, 'Problem parameters: '
    print *, 'length of the pendulum (m) : ', l
    print *, 'time step (s) : ', dt
    print *, 'number of steps : ', nm
    print *, 'initial angle (rad) : ', thetain
    print *, 'initial angular velocity (rad/s) : ', omegain

    ! Solution of the movement equation using the Euler explicit method.
    time_loop: do n = 1,nm-1
        omega(n+1) = omega(n) - g/l*sin(theta(n))*dt
        theta(n+1) = theta(n) + omega(n)*dt
    enddo time_loop

    ! Save solution in a text file
    open(unit=1, file='pendulum.txt')
    output_loop: do n = 1,nm
        write(unit=1,fmt=*) (n-1)*dt, theta(n), omega(n)
    enddo output_loop

    close(unit=1)

end program pendulum