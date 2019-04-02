program int_selector
    !
    ! Program that samples which interaction the particle will follow, based on 
    ! the total and interaction cross sections.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
implicit none

    ! Transport parameters
    integer, parameter :: nhist = 1000  ! number of particle histories simulated

    ! Geometry parameters
    real, parameter :: sigma_t = 1.0            ! total interaction cross section (cm-1)
    real, parameter :: sigma_a = 0.9*sigma_t    ! absorption cross section (cm-1)
    real, parameter :: sigma_s = 0.1*sigma_t    ! scattering cross section (cm-1)

    integer :: i
    real :: rnno
    integer, dimension(2) :: ninteraction  = 0     ! counter holding the number of interactions for each process

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Get distance to next interaction
    do i = 1,nhist
         rnno = rng_set()
         
         ! Select interaction
         if(rnno <= sigma_a/sigma_t) then
            ! absorption
            ninteraction(1) = ninteraction(1) + 1
         else
            ! scattering
            ninteraction(2) = ninteraction(2) + 1
         endif

    enddo

    ! Print results
    write(*,*) 'Fraction of particles absorbed = ', real(ninteraction(1))/sum(ninteraction)
    write(*,*) 'Fraction of particles scattered = ', real(ninteraction(2))/sum(ninteraction)

end program int_selector