program rng_sampling
    !
    ! Program used to sample the RANDU PRNG.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int64, real64
    
    use mod_rng, only: rng_init => lcg_init, rng_get => lcg_get, rng_get_int => lcg_get_int

    implicit none

    integer(kind=int64), parameter :: seed = 1          ! parameters of the rng
    integer(kind=int64), parameter :: a_in = 65539 
    integer(kind=int64), parameter :: b_in = 0
    integer(kind=int64), parameter :: m_in = 16777216   ! 2**24

    ! We want to store integer and real versions of the random numbers. Set sequence 
    ! size to maximum theoretical value.
    integer(kind=int64), parameter :: nsamples = m_in
    integer(kind=int64), dimension(nsamples) :: int_samples = 0
    real(kind=real64), dimension(nsamples) :: samples = 0.0

    integer(kind=int64) :: isample, nseeds, npoints 
    real(kind=real64) :: sample_mean

    ! Initialize the PRNG.
    call rng_init(a_in, b_in, m_in, seed)

    write(6,'(A, I10)') 'Number of samples : ', nsamples

    ! Get samples.
    do isample = 1,nsamples
        
        samples(isample) = rng_get()
        
    enddo

    ! Calculate mean of the sample.
    sample_mean = sum(samples)/size(samples)

    write(6,'(A, F10.3)') 'Random sequence mean : ', sample_mean

    ! Now arrange the generated random numbers in 3-tuples and output to file.
    npoints = floor(nsamples/3.0)

    open(unit=1, file='rng_sequence_plot.txt')

    do isample = 1,npoints,3
        write(unit=1, fmt='(F10.5, F10.5, F10.5)') samples(isample), samples(isample+1), samples(isample+2)        
    enddo

    close(unit=1)

    ! Re-Initialize the PRNG and get integer random numbers.
    call rng_init(a_in, b_in, m_in, seed)

    ! Get samples.
    do isample = 1,nsamples
        
        int_samples(isample) = rng_get_int()
        
    enddo

    ! Obtain period of the sequence. For this, get the number of times that the 
    ! seed is found on the random sequence.
    nseeds = count(int_samples .eq. seed)

    write(6,'(A, I10)') 'Number of times the seed is repeated : ', nseeds
    write(6,'(A, I10, A, F7.3, A)') 'Period of the random sequence : ', m_in/nseeds, ' (', 100.0/nseeds, ' % M )'


end program rng_sampling