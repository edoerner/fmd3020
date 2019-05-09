module mod_rng
    !
    ! Module used to define pseudo-random numer generators.
    !

    use iso_fortran_env, only: int64, real64

implicit none

    integer(kind=int64) :: xold     ! last random number in the sequence.
    integer(kind=int64) :: a, b     ! integer parameters for LCG generators.
    integer(kind=int64) :: M
    
contains
    !
    ! Module used to define pseudo-random numer generators.
    !
    subroutine randu_init(seed)
        ! Initialize the rng. It saves the seed and parameters to 
        ! internal variables in the module.
        integer(kind=int64), intent(in)  :: seed
                
        ! The first element in params is always the seed of the PRNG
        xold = seed

        ! The rest of the arguments are the integer parameters of the PRNG.
        a = 65539
        b = 0
        M = 2147483648

        ! Print some information to console.
        write(*,'(A)') 'RANDU PRNG with parameters: '
        write(*,'(A,I10,A,I10,A,I15)') 'a = ', a, ', b = ', b, ', M = ', M
        write(*,'(A,I10)') 'x 0 = ', xold

    end subroutine randu_init

    integer(kind=int64) function randu_get()
    
        integer(kind=int64) :: xnew

        xnew = mod(a*xold + b, M)
        xold = xnew
        randu_get = xold
        
    end function randu_get

    real(kind=real64) function randu_eta_get()

        integer(kind=int64) :: xnew
        
        xnew = randu_get()
        randu_eta_get = (1.0_real64*xnew)/(1.0_real64*M)

    end function randu_eta_get

end module mod_rng