program dice_sampling
    !
    ! Program that generates uniform integer numbers between 1 and 6, simulating
    ! a dice throwing. This program is used to test the central limit theorem.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng

    implicit none

    integer, parameter :: nbatch = 1500     ! number of batches
    integer, parameter :: nperbatch = 10    ! throwings per batch

    real, dimension(nbatch) :: sum_score    ! sum of the throwings in each batch
    real, dimension(nbatch) :: sum2_score   ! sum of the throwings in each batch (squared)

    integer :: ibatch, isample
    real :: score, score2
    real :: throw

    ! Initialize the PRNG
    call rng_init(20180815)

    ! Proceed with the sampling process. Start the batch loop
    batch_loop: do ibatch = 1,nbatch   
        
        ! Initialize scoring variable
        score = 0.0
        score2 = 0.0
    
        sampling_loop: do isample = 1,nperbatch
            ! Sampling process. Accumulate the sample and its squared value.
            throw = dice()
            score = score + throw
            score2 = score2 + throw**2
            
        enddo sampling_loop

        ! Accumulate results and proceed to next batch.
        sum_score(ibatch) = score
        sum2_score(ibatch) = score2

        ! Print mean and variance of current batch.
        write(unit=6, fmt='(I5, F10.5, F10.5)') ibatch, score/nperbatch, &
            (score2 - nperbatch*(score/nperbatch)**2)/(nperbatch-1)
    
    enddo batch_loop

    ! Save results to a file
    open(unit=1, file='dice_sampling.txt')
    do ibatch = 1,nbatch
        write(unit=1, fmt='(F10.5)') sum_score(ibatch)       
    enddo
    close(1)

    ! Now obtain the mean across all batches and the associated uncertainty.
    score = sum(sum_score/nperbatch)/nbatch
    score2 = sum((sum_score/nperbatch - score)**2)/(nbatch-1)
    write(unit=6, fmt='(A, F10.5, A, F10.5)') 'sample mean = ', score, &
            ', sample variance = ', score2
    
    
    contains

    integer function dice()
        ! Function that generates uniform integer numbers between 1 and 6, simulating
        ! a dice throwing.    

        real :: rnno

        ! Generate an uniform random number in [0,1)
        rnno = rng_set()

        dice = 1 + floor(6.0*rnno)
        
    end function dice
end program dice_sampling