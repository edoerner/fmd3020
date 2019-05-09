program rng_sampling
    !
    ! Program used to sample the RANDU PRNG.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int64, real64
    
    use mod_rng, only: rng_init => randu_init, rng_get => randu_eta_get

    implicit none

    integer(kind=int64), parameter :: nsamples = 10000  ! number of pairs that will be sampled.
    integer(kind=int64), parameter :: seed = 1          ! parameters of the rng.
    integer(kind=int64) :: isample 

    real(kind=real64) :: xi, yi, zi

    ! Initialize the PRNG
    call rng_init(seed)

    ! Sampling loop. Print sample to screen.
    open(unit=1, file='randu_2D.txt')

    do isample = 1,nsamples
        xi = rng_get()
        yi = rng_get()
        write(unit=1,fmt='(F10.5, F10.5)') xi, yi
    enddo

    close(unit=1)

    open(unit=2, file='randu_3D.txt')

    do isample = 1,nsamples
        xi = rng_get()
        yi = rng_get()
        zi = rng_get()
        write(unit=2,fmt='(F10.5, F10.5, F10.5)') xi, yi, zi
    enddo

    close(unit=2)

end program rng_sampling