module mod_rng
    !
    ! Definition of the pseudo-random number generator. Fortran intrinsic RNG.
    !
implicit none

contains

    subroutine rng_init(seed)
        ! Initialize the rng. Obtained from https://stackoverflow.com/a/51898794/10919904
        integer,intent(in)  :: seed
        
        integer, parameter :: wp = selected_real_kind(15,307)
        integer, parameter :: n_discard = 100

        integer :: state_size, i
        integer, allocatable, dimension(:) :: state
        real(wp) :: ran, oldran

        call random_seed( size=state_size )
        write(*,*) 'Initializing intrinsic PRNG -- state size is: ', state_size
        write(*,*) 'Number of discards : ', n_discard

        allocate(state(state_size))

        ! Simple method of initializing seed from single scalar
        state = seed
        call random_seed( put=state )

        ! 'Prime' the generator by pulling the first few numbers
        ! In reality, these would be discarded but I will print them for demonstration
        ran = 0.5_wp
        do i=1,n_discard
            oldran = ran
            call random_number(ran)
            write(*,'(a,i3,2es26.18)') 'iter, ran, diff: ', i, ran, ran-oldran
        enddo

    end subroutine rng_init

    real function rng_set()
    
        call random_number(rng_set)
        
    end function rng_set

end module mod_rng