program mfp_sampling
    !
    ! Program that samples the free flight, i.e. distance to the next interaction.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_mfp

    implicit none

    ! Transport parameters
    integer, parameter :: nhist = 1000  ! number of particle histories simulated

    ! Geometry parameters
    real :: sigma_t = 1.0   ! total interaction cross section (cm-1)

    integer :: i
    real, dimension(nhist) :: tstep       ! distance to next interaction, before check with region boundaries
    
    ! Initialize the PRNG
    call rng_init(20180815)

    ! Get distance to next interaction
    do i = 1,nhist
        tstep(i) = mfp(sigma_t)        
    enddo
    
    ! Save results to file
    open(unit=1, file='mfp_sampling.txt')
    write(unit=1, fmt='(F10.5)') tstep
    close(1)

end program mfp_sampling