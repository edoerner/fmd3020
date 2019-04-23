program mfp_sampling
    !
    ! Program that samples the free flight, i.e. distance to the next interaction. and estimates its 
    ! mean value and variance for different number of histories.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_mfp

    implicit none

    ! Transport parameters
    integer, parameter :: nbatch = 10   ! number of statistical batches
    integer, parameter :: ncase = 9    ! number of cases that will be studied
    integer, parameter, dimension(ncase) :: nperbatch = [1, 2, 5, 10, 25, 50, 100, 200, 500]

    ! Geometry parameters
    real :: sigma = 2.0   ! total interaction cross section (cm-1)

    integer :: icase, ibatch, isample
    real :: path, score, score2           ! auxiliary variables used for scoring
    real, allocatable :: step(:)    ! distance to next interaction (cm)
    real, allocatable :: step2(:)   ! distance to next interaction (squared)

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Allocate the arrays that will hold the samples.
    allocate(step(size(nperbatch)))
    allocate(step2(size(nperbatch)))

    step = 0.0
    step2 = 0.0

    ! Start the sampling process. We calculate the mean and variance of the sample for 
    ! each nperbatch value given by the user.
    ncase_loop: do icase = 1,ncase        
        
        ! Proceed with the sampling process. Start the batch loop
        batch_loop: do ibatch = 1,nbatch   
        
            ! Initialize scoring variable
            score = 0.0
            score2 = 0.0
            
            sampling_loop: do isample = 1,nperbatch(icase)
                ! Sampling process. Accumulate the sample and its squared value.
                path = mfp(sigma)
                score = score + path
                score2 = score2 + path**2
                
            enddo sampling_loop

            ! Accumulate results and proceed to next batch.
            step(icase) = step(icase) + score
            step2(icase) = step2(icase) + score2
            
        enddo batch_loop

        ! Statistical analysis.
        step(icase) = step(icase)/(nbatch*nperbatch(icase))
        step2(icase) = (step2(icase) - (nbatch*nperbatch(icase))*step(icase)**2)/(((nbatch*nperbatch(icase))-1))

        write(unit=6, fmt='(I2, I5, F10.5, F10.5)') icase, nperbatch(icase), step(icase), step2(icase)
        
    enddo ncase_loop
    
    ! Save results to file
    open(unit=1, file='mfp_sampling.txt')
    do icase = 1,ncase
        write(unit=1, fmt='(I5, F10.5, F10.5)') nperbatch(icase), step(icase), step2(icase)        
    enddo
    close(1)

end program mfp_sampling