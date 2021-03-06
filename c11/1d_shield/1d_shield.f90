program shield_1d
    !
    ! Program that simulates the particle transport through a 1D shield using Monte Carlo methods.
    ! Multi-region version
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use mod_rng
    use mod_mfp
    use mod_scatt

implicit none

    ! Geometry parameters
    integer, parameter :: nreg = 3                              ! number of region in the 1D shield
    real, dimension(nreg) :: xthick = (/0.1, 0.1, 0.1/)         ! thickness of the 1D shield (cm)
    real, dimension(nreg+1) :: xbounds = 0.0                    ! boundaries of the shield regions (cm)

    ! Transport parameters
    real, dimension(nreg) :: sigma_t = (/0.1, 10.0, 100.0/)     ! Total interaction cross section (cm-1)
    real, dimension(nreg) :: sigma_a = (/0.01, 0.1, 10.0/)      ! Absorption cross section (cm-1)
    real, dimension(nreg) :: sigma_s = (/0.09, 9.9, 90.0/)      ! Scattering cross section (cm-1)

    ! Particle parameters
    real, parameter :: xin = 0.0        ! initial position (cm)
    real, parameter :: uin = 1.0        ! initial direction
    real, parameter :: wtin = 1.0       ! statistical weight
    integer, parameter :: irin = 1         ! initial region

    integer :: ir
    real :: x, u, wt

    ! Simulation parameters
    integer, parameter :: nperbatch = 500    
    integer, parameter :: nbatch = 10               ! number of statistical batches
    integer, parameter :: nhist = nbatch*nperbatch  ! number of simulated histories

    ! Scoring variables
    real, dimension(0:nreg+1, nbatch) :: score = 0.0    ! score(0) : reflection
                                                        ! score(1:nreg) : absorption
                                                        ! score(nreg+1) : transmission

    real, dimension(0:nreg+1) :: mean_score = 0.0       ! mean value of scoring array
    real, dimension(0:nreg+1) :: unc_score = 0.0        ! uncertainty values of scoring array

    integer :: i, ihist, ibatch     ! loop counters 
    logical :: pdisc                ! flag to discard a particle
    integer :: irnew                ! index of new region 
    real :: pstep                   ! distance to next interaction 
    real :: dist                    ! distance to boundary along particle direction
    real :: rnno 
    real :: start_time, end_time

    call cpu_time(start_time)

    ! Initialize RNG
    call rng_init(20180815)

    ! Print some info to console
    write(*,'(A, I15)') 'Number of histories : ', nhist
    write(*,'(A, I15)') 'Number of batches : ', nbatch
    write(*,'(A, I15)') 'Number of histories per batch batches : ', nperbatch

    do i = 1,nreg
        write(*,'(A, I5, A)') 'For region ', i, ':'
        write(*,'(A, F15.5)') 'Total interaction cross-section (cm-1): ', sigma_t(i)
        write(*,'(A, F15.5)') 'Absoption interaction cross-section (cm-1): ', sigma_a(i)
        write(*,'(A, F15.5)') 'Scattering interaction cross-section (cm-1): ', sigma_s(i)
        write(*,'(A, F15.5)') '1D shield thickness (cm): ', xthick(i)
    enddo
    
    ! Initialize simulation geometry.
    write(*,'(A)') 'Region boundaries (cm):'
    write(*,'(F15.5)') xbounds(1)
    do i = 1,nreg
        xbounds(i+1) = xbounds(i) + xthick(i)
        write(*,'(F15.5)') xbounds(i+1)
    enddo
    
    ibatch_loop: do ibatch = 1,nbatch
        ihist_loop: do ihist = 1,nperbatch
            
            ! Initialize history
            x = xin
            u = uin
            wt = wtin
            ir = irin
            
            ! Set flag used for particle discard.
            pdisc = .false.

            ! Enter transport process
            particle_loop: do
            
                ptrans_loop: do

                    ! Distance to the next interaction.
                    if(ir == 0 .or. ir == nreg+1) then
                        ! Vacuum step
                        pstep = 1.0E8
                    else
                        pstep = mfp(sigma_t(ir))
                    endif

                    ! Save particle current region.
                    irnew = ir

                    ! Check expected particle step with geometry.
                    if(u .lt. 0.0) then
                        if(irnew == 0) then
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                        else
                            ! The particle goes to the front face of the shield.
                            dist = (xbounds(ir) - x)/u

                            ! Now check if the particle leaves current region.
                            if(dist < pstep) then
                                pstep = dist                      
                                irnew = irnew - 1
                            endif
                        endif
                        
                    else if(u .gt. 0.0) then
                        if(irnew == nreg+1) then
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                        else
                            ! The particle goes to the back face of the shield.
                            dist = (xbounds(ir+1) - x)/u

                            ! Now check if the particle leaves current region.
                            if(dist < pstep) then
                                pstep = dist
                                irnew = irnew + 1
                            endif
                        endif
                        
                    endif

                    ! Check if particle has been discarded.
                    if(pdisc .eqv. .true.) then
                        exit
                    endif

                    ! Update position of particle.
                    x = x + pstep*u

                    ! Reached this point, if the particle has not changed region, it must interact.
                    if(ir == irnew) then
                        exit
                    else
                        ! Update particle region index and continue transport process.
                        ir = irnew
                    endif

                enddo ptrans_loop    
            
                if(pdisc .eqv. .true.) then
                    ! The particle has been discarded. Score it and end tracking.
                    score(ir, ibatch) = score(ir, ibatch) + wt
                    exit
                endif

                ! Determine interaction type.
                rnno = rng_set()
                if(rnno .le. sigma_a(ir)/sigma_t(ir)) then
                    ! The particle was absorbed.
                    score(ir, ibatch) = score(ir, ibatch) + wt
                    exit
                else
                    ! The particle was scattered, get new direction to re-enter transport loop.
                    u = scatt(u)
                endif
            
            enddo particle_loop
        enddo ihist_loop

    enddo ibatch_loop

    ! Save absorption probability to file.
    open(unit=1, file='absorption.txt')
    do i=1,nbatch
        write(unit=1, fmt='(I3, F10.5)') i, sum(score(1:nreg,i))/nperbatch
    enddo
    close(1)

    ! Statistical procesing 
    mean_score = sum(score, dim = 2)/nbatch

    unc_score = sum(score**2, dim = 2) - nbatch*mean_score**2
    unc_score = unc_score/(nbatch*(nbatch - 1))
    unc_score = sqrt(unc_score)

    mean_score = mean_score/nperbatch
    unc_score = unc_score/nperbatch

    ! Print results to console.
    write(*,'(A,F15.5,A,F15.5,A)') 'Reflection : ', mean_score(0), ' +/-', &
        100*unc_score(0)/mean_score(0), '%'

    ! To calculate the uncertainty of the absorption we need to combine the uncertainty of the deposition 
    ! in each region.
    write(*,'(A,F15.5,A,F15.5,A)') 'Absorption : ', sum(mean_score(1:nreg)), ' +/-', &
        100*sqrt(sum(unc_score(1:nreg)**2))/sum(mean_score(1:nreg)), '%'
    
    write(*,'(A,F15.5,A,F15.5,A)') 'Transmission : ', mean_score(nreg+1), ' +/-', &
        100*unc_score(nreg+1)/mean_score(nreg+1), '%'

    ! Get elapsed time
    call cpu_time(end_time)
    write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time

end program shield_1d