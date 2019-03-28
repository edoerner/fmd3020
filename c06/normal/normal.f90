program normal
    !
    ! Program that samples random numbers from the distribution f(x) = e(-x^2/2) from 
    ! random numbers in [0,1) using the rejection method.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
implicit none

    integer, parameter :: nsamples = 100000
    real, dimension (nsamples) :: samples_box, samples_exp

    real, parameter :: a = 0.0, b = 10.0

    integer :: i, ntries_box = 0, ntries_exp    ! loop counters
    real :: xsample                             ! sampling candidate
    real :: rnno1, rnno2

    real :: pi, pmax, emax

    pi = 4.0*atan(1.0)
    pmax = sqrt(2.0/pi)
    emax = sqrt(2.0*exp(1.0)/pi)

    ! Initialize PRNG
    call random_seed

    ! First sample using a box x in [a,b], y in [0,pmax]
    sample_box: do i = 1,nsamples
        reject_box: do
            ntries_box = ntries_box + 1
            ! Get sampling candidate
            call random_number(rnno1)
            xsample = a + rnno1*(b-a)

            ! Decide if candidate is accepted or not
            call random_number(rnno2)
            if(pmax*rnno2 <= norm_dist(xsample)) then
                ! accept sample                
                exit
            endif           
        enddo reject_box
        ! Store sample
        samples_box(i) = xsample
    enddo sample_box

    ! Now sample using an exponencial as envelope
    sample_exp: do i = 1,nsamples
        reject_exp: do
            ntries_exp = ntries_exp + 1
            ! Get sampling candidate
            call random_number(rnno1)
            xsample = -log(1.0-rnno1)

            ! Decide if candidate is accepted or not
            call random_number(rnno2)
            if(rnno2*emax*exp(-xsample) <= norm_dist(xsample)) then
                ! accept sample                
                exit
            endif           
        enddo reject_exp
        ! Store sample
        samples_exp(i) = xsample
    enddo sample_exp

    ! Save results to file
    open(unit=1, file='normal_box.txt')
    write(unit=1, fmt='(F10.5)') samples_box
    close(1)

    ! Save results to file
    open(unit=2, file='normal_exp.txt')
    write(unit=2, fmt='(F10.5)') samples_exp
    close(2)

    ! Print some info to console
    print *, 'Number of samples : ', nsamples
    print *, 'Number of tries (box): ', ntries_box
    print *, 'Ratio : ', real(nsamples)/ntries_box
    print *, 'Number of tries (exp): ', ntries_exp
    print *, 'Ratio : ', real(nsamples)/ntries_exp

contains

    real function norm_dist(x)
        real, intent(in) :: x
        
        norm_dist = sqrt(2.0/pi)*exp(-0.5*x**2)
        
    end function norm_dist

end program normal