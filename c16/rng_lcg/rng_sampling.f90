program rng_sampling
    !
    ! Program that generates uniform integer numbers between 1 and 6, simulating
    ! a dice throwing. This program is used to test the central limit theorem.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int32
    
    use mod_rng

    implicit none

    integer(kind=int32), parameter :: nsamples = 16     ! number of batches.
    integer(kind=int32), dimension(4) :: lcg_params        ! parameters of the rng.
    integer(kind=int32) :: isample, sample

    ! Define the parameters that will define the PRNG.
    lcg_params(1) = 1   ! seed
    lcg_params(2) = 5   ! a
    lcg_params(3) = 1   ! b
    lcg_params(4) = 16  ! M

    ! Initialize the PRNG
    call rng_init(lcg_params)

    ! Sampling loop. Print sample to screen.
    do isample = 1,nsamples
        sample = rng_set()
        write(*,'(A,I2,A,I5)') 'x', isample, ' = ', sample
    enddo

end program rng_sampling