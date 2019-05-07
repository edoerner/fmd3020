module mod_rng
    !
    ! Module used to define pseudo-random numer generators.
    !

    use iso_fortran_env, only: int32

implicit none

    integer(kind=int32) :: xold     ! last random number in the sequence.
    integer(kind=int32) :: a, b     ! integer parameters for LCG generators.
    integer(kind=int32) :: M

contains
    !
    ! Module used to define pseudo-random numer generators.
    !
    subroutine rng_init(params)
        ! Initialize the rng. It saves the seed and parameters to 
        ! internal variables in the module.
        integer(kind=int32), intent(in)  :: params(:)
        
        ! The first element in params is always the seed of the PRNG
        xold = params(1)

        ! The rest of the arguments are the integer parameters of the PRNG.
        a = params(2)
        b = params(3)
        M = params(4)

        ! Print some information to console.
        write(*,'(A)') 'Linear congruential generator with parameters: '
        write(*,'(A,I5,A,I5,A,I5)') 'a = ', a, ', b = ', b, ', M = ', M
        write(*,'(A,I5)') 'x 0 = ', xold

    end subroutine rng_init

    integer(kind=int32) function rng_set()
    
        integer(kind=int32) :: xnew

        xnew = mod(a*xold + b, M)
        xold = xnew
        rng_set = xold
        
    end function rng_set

end module mod_rng