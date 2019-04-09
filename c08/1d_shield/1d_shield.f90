program shield_1d
    !
    ! Program that simulates the particle transport through a 1D shield using Monte Carlo methods.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_mfp
    use mod_scatt

implicit none

    ! Transport parameters
    real, parameter :: sigma_t = 1.0    ! Total interaction cross section (cm-1)
    real, parameter :: sigma_a = 0.9    ! Absorption cross section (cm-1)
    real, parameter :: sigma_s = 0.1    ! Scattering cross section (cm-1)

    ! Geometry parameters
    real, parameter :: xthick = 5.0     ! thickness of the 1D shield (cm)

    ! Particle parameters
    real, parameter :: xin = 0.0
    real, parameter :: uin = 1.0
    real, parameter :: wtin = 1.0

    real :: x, u, wt

    ! Simulation parameters
    integer, parameter :: nhist = 10000000

    ! Scoring variables
    real, dimension(3) :: score = 0.0   ! score(1) : reflection
                                        ! score(2) : absorption
                                        ! score(3) : transmission

    integer :: ihist 
    real :: pstep   ! distance to next interaction  
    real :: rnno 
    real :: start_time, end_time

    call cpu_time(start_time)

    ! Initialize RNG
    call rng_init(20180815)

    ! Print some info to console
    write(*,'(A, I15)') 'Number of histories : ', nhist
    write(*,'(A, F15.5)') 'Total interaction cross-section (cm-1): ', sigma_t
    write(*,'(A, F15.5)') 'Absoption interaction cross-section (cm-1): ', sigma_a
    write(*,'(A, F15.5)') 'Scattering interaction cross-section (cm-1): ', sigma_s
    write(*,'(A, F15.5)') '1D shield thickness (cm): ', xthick

    ihist_loop: do ihist = 1,nhist
        
        ! Initialize history
        x = xin
        u = uin
        wt = wtin

        ! Enter transport process
        transport_loop: do
            
            ! Distance to the next interaction.
            pstep = mfp(sigma_t)

            ! Update position of particle.
            x = x + pstep*u

            ! Check position of particle with geometry.
            if(x .lt. 0.0) then
                ! The particle was reflected.
                score(1) = score(1) + wt
                exit
            else if(x .gt. xthick) then
                ! The particle was transmitted.
                score(3) = score(3) + wt
                exit
            endif

            ! Determine interaction type.
            rnno = rng_set()
            if(rnno .le. sigma_a/sigma_t) then
                ! The particle was absorbed.
                score(2) = score(2) + wt
                exit
            else
                ! The particle was scattered, get new direction.
                u = scatt(u)
            endif
            
        enddo transport_loop
    enddo ihist_loop

    ! Print results to console
    write(*,'(A20,F15.5)') 'Reflection : ', score(1)/real(nhist)
    write(*,'(A20,F15.5)') 'Absorption : ', score(2)/real(nhist)
    write(*,'(A20,F15.5)') 'Transmission : ', score(3)/real(nhist)

    ! Get elapsed time
    call cpu_time(end_time)
    write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time

end program shield_1d