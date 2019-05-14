module mod_rng
    !
    ! Module used to define pseudo-random numer generators.
    !

    use iso_fortran_env, only: int64, real64

implicit none

    integer(kind=int64) :: xold     ! last random number in the sequence.

    private :: a, b, M              ! limit access of LCG parameters within module.
    integer(kind=int64) :: a, b     ! integer parameters for LCG generators.
    integer(kind=int64) :: M
    
contains
    !
    ! Module used to define pseudo-random numer generators.
    !
    subroutine lcg_init(a_in, b_in, m_in, seed)
        ! Initialize the rng. It saves the seed and parameters to 
        ! internal variables in the module.
        integer(kind=int64), intent(in)  :: a_in, b_in, m_in, seed
                
        ! The first element in params is always the seed of the PRNG
        xold = seed

        ! The rest of the arguments are the integer parameters of the PRNG.
        a = a_in
        b = b_in
        M = m_in

        ! Print some information to console.
        write(*,'(A)') 'LCG PRNG with parameters: '
        write(*,'(A,I10,A,I10,A,I15)') 'a = ', a, ', b = ', b, ', M = ', M
        write(*,'(A,I10)') 'x 0 = ', xold

    end subroutine lcg_init

    integer(kind=int64) function lcg_get_int()
    
        integer(kind=int64) :: xnew

        xnew = mod(a*xold + b, M)
        xold = xnew
        lcg_get_int = xold
        
    end function lcg_get_int

    real(kind=real64) function lcg_get()

        integer(kind=int64) :: xnew
        
        xnew = lcg_get_int()
        lcg_get = (1.0_real64*xnew)/(1.0_real64*M)

    end function lcg_get

end module mod_rng