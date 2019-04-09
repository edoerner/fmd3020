program scatt_angle
    !
    ! Program that samples the scattering angle of the particle, based on 
    ! the scattering interaction cross section.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_scatt

implicit none

    ! Transport parameters
    integer, parameter :: nhist = 100000  ! number of particle histories simulated

    ! Geometry parameters
    real, parameter :: sigma_s = 0.1   ! scattering interaction cross section (cm-1)

    ! Particle parameters
    real, parameter :: uin = 1.0

    integer :: i
    real, dimension(nhist) :: u ! direction (cosine in x) of the particle

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Get distance to next interaction
    do i = 1,nhist
        u(i) = scatt(uin)        
    enddo

    ! Save results to file
    open(unit=1, file='scatt_sampling.txt')
    write(unit=1, fmt='(F10.5)') u
    close(1)

end program scatt_angle