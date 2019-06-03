program shield_1d
    !
    ! Program that simulates the particle transport through a 1D shield using Monte Carlo methods.
    ! Multi-region version. This version incorporates geometrical splitting.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !
    use iso_fortran_env, only: int32, real64

    use mod_rng
    use mod_mfp
    use mod_scatt

implicit none

    ! Geometry parameters
    integer(kind=int32), parameter :: nreg = 5              ! number of region in the 1D shield
    real(kind=real64), dimension(nreg) :: xthick = 2.0      ! thickness of each region (cm)
    real(kind=real64), dimension(nreg+1) :: xbounds = 0.0   ! boundaries of the shield regions (cm)

    ! Transport parameters
    real(kind=real64), dimension(nreg) :: sigma_t = 1.0     ! Total interaction cross section (cm-1)
    real(kind=real64), dimension(nreg) :: sigma_a = 0.8     ! Absorption cross section (cm-1)
    real(kind=real64), dimension(nreg) :: sigma_s = 0.2     ! Scattering cross section (cm-1)

    ! Particle parameters
    real(kind=real64), parameter :: xin = 0.0       ! initial position (cm)
    real(kind=real64), parameter :: uin = 1.0       ! initial direction
    real(kind=real64), parameter :: wtin = 1.0      ! statistical weight
    integer(kind=int32), parameter :: irin = 1      ! initial region

    ! Now we introduce the concept of stack, therefore we expand each particle 
    ! parameter to a size of the stack.
    integer(kind=int32), parameter :: nstack = 100  ! maximum number of particles that can hold the stack
    integer(kind=int32) :: ip                       ! current index of the stack
    integer(kind=int32), dimension(nstack) :: ir
    real(kind=real64), dimension(nstack) :: x, u, wt

    ! Geometrical splitting VRT parameters
    real(kind=real64) :: gs_r = 2.63                ! ratio of region importances.
    integer(kind=int32) :: gs_n                     ! integer part of gs_r
    real(kind=real64), dimension(nreg) :: gs_i      ! importance of each region
    
    ! Simulation parameters
    integer(kind=int32), parameter :: nperbatch = 1E6    
    integer(kind=int32), parameter :: nbatch = 10               ! number of statistical batches
    integer(kind=int32), parameter :: nhist = nbatch*nperbatch  ! number of simulated histories

    ! Scoring variables
    real(kind=real64), dimension(0:nreg+1) :: score, score2 = 0.0    ! score(0) : reflection
                                                        ! score(1:nreg) : absorption
                                                        ! score(nreg+1) : transmission

    real(kind=real64), dimension(0:nreg+1) :: mean_score = 0.0   ! mean value of scoring array
    real(kind=real64), dimension(0:nreg+1) :: unc_score = 0.0    ! uncertainty values of scoring array

    integer(kind=int32) :: i, ihist, ibatch     ! loop counters 
    logical :: pdisc                ! flag to discard a particle
    integer(kind=int32) :: irnew                ! index of new region 
    real(kind=real64) :: pstep                   ! distance to next interaction 
    real(kind=real64) :: dist                    ! distance to boundary along particle direction
    real(kind=real64) :: rnno 
    real(kind=real64) :: max_var, fom
    real(kind=real64) :: start_time, end_time

    call cpu_time(start_time)

    ! Initialize RNG
    call rng_init(20180815)

    ! Print some info to console
    write(*,'(A, I15)') 'Number of histories : ', nhist
    write(*,'(A, I15)') 'Number of batches : ', nbatch
    write(*,'(A, I15)') 'Number of histories per batch : ', nperbatch

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

    ! Initialize data for geometrical splitting VRT.
    gs_i(1) = gs_r
    do i = 2,nreg
        gs_i(i) = gs_r*gs_i(i-1)            
    enddo

    write(*,'(A)') '(GS VRT ) Importances for each region:'
    do i = 1,nreg        
        write(*,'(A, I3, A, F15.5)') 'nreg = ', i, ', I = ', gs_i(i)
    enddo
    
    ibatch_loop: do ibatch = 1,nbatch
        ihist_loop: do ihist = 1,nperbatch
            
            ! Initialize stack
            ip = 1
            x(ip) = xin
            u(ip) = uin
            wt(ip) = wtin
            ir(ip) = irin
            
            ! Set flag used for particle discard.
            pdisc = .false.

            ! Enter transport process
            particle_loop: do
            
                ptrans_loop: do

                    ! Distance to the next interaction.
                    if(ir(ip) == 0 .or. ir(ip) == nreg+1) then
                        ! Vacuum step
                        pstep = 1.0E8
                    else
                        pstep = mfp(sigma_t(ir(ip)))
                    endif

                    ! Save particle current region.
                    irnew = ir(ip)

                    ! Check expected particle step with geometry.
                    if(u(ip) .lt. 0.0) then
                        if(irnew == 0) then
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                        else
                            ! The particle goes to the front face of the shield.
                            dist = (xbounds(ir(ip)) - x(ip))/u(ip)

                            ! Now check if the particle leaves current region.
                            if(dist < pstep) then
                                pstep = dist                      
                                irnew = irnew - 1
                            endif
                        endif
                        
                    else if(u(ip) .gt. 0.0) then
                        if(irnew == nreg+1) then
                            ! The particle is leaving the geometry, discard it.
                            pdisc = .true.
                        else
                            ! The particle goes to the back face of the shield.
                            dist = (xbounds(ir(ip)+1) - x(ip))/u(ip)

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
                    x(ip) = x(ip) + pstep*u(ip)

                    ! Reached this point, if the particle has not changed region, it must interact.
                    if(ir(ip) == irnew) then
                        exit
                    else
                        ! Implementation of geometrical splitting VRT.
                        gs_r = gs_i(irnew)/gs_i(ir(ip))

                        if(gs_r > 1.0) then
                            ! Particle is going to the ROI. Perform particle splitting.
                            gs_n = floor(gs_r)
                            rnno = rng_set()

                            if(rnno .le. 1.0 - (gs_r - gs_n)) then
                                ! Divide the particle in n particles

                            else
                                ! Divide the particle in n+1 particles

                            endif
                        else
                            ! Particle is going backwards the ROI. Perform russian roulette.
                            rnno = rng_set()
                            if(rnno .gt. gs_r) then
                                ! particle does not survive, finish particle
                                ip = ip - 1
                                cycle
                            else
                                ! particle survives, adjusts weight accordingly and continue transport.
                                wt(ip) = (1.0/gs_r)*wt(ip)                        
                            endif
                        endif
                        
                        ! Update particle region index and continue transport process.
                        ir(ip) = irnew
                    endif

                enddo ptrans_loop    
            
                if(pdisc .eqv. .true.) then
                    ! The particle has been discarded. Score it and decrease stack.
                    score(ir(ip)) = score(ir(ip)) + wt(ip)
                    ip = ip - 1

                    if(ip .eq. 0) then
                        ! Finish particle history.
                        exit
                    else
                        ! Start transport of the next particle.
                        cycle
                    endif
                endif

                ! Determine interaction type.
                rnno = rng_set()
                if(rnno .le. sigma_a(ir(ip))/sigma_t(ir(ip))) then
                    ! The particle was absorbed. Score it and decrease stack.
                    score(ir(ip)) = score(ir(ip)) + wt(ip)
                    ip = ip - 1

                    if(ip .eq. 0) then
                        ! Finish particle history.
                        exit
                    else
                        ! Start transport of the next particle.
                        cycle
                    endif                   
                else
                    ! The particle was scattered, get new direction to re-enter transport loop.
                    u(ip) = scatt(u(ip))
                endif
            
            enddo particle_loop
            
        enddo ihist_loop

        ! Accumulate results and reset score arrays.
        mean_score = mean_score + score
        unc_score = unc_score + score**2

        score = 0.0
        score2 = 0.0

    enddo ibatch_loop

    ! Statistical procesing
    mean_score = mean_score/nbatch
    unc_score = (unc_score - nbatch*mean_score**2)/(nbatch-1)
    unc_score = unc_score/nbatch
    unc_score = sqrt(unc_score)

    ! Calculate relative uncertainty, to be used for the FOM.
    unc_score = unc_score/mean_score

    ! Print results to console.
    write(*,'(A,F10.5,A,F10.5,A)') 'Reflection : ', mean_score(0)/nperbatch, ' +/-', 100.0*unc_score(0), '%'

    ! To calculate the uncertainty of the absorption we need to combine the uncertainty of the deposition 
    ! in each region.
    write(*,'(A,F10.5,A,F10.5,A)') 'Absorption : ', sum(mean_score(1:nreg))/nperbatch, ' +/-', & 
        100.0*sum(unc_score(1:nreg)), '%'
    
    write(*,'(A,F10.5,A,F10.5,A)') 'Transmission : ', mean_score(nreg+1)/nperbatch, ' +/-', &
        100.0*unc_score(nreg+1), '%'

    ! Get elapsed time
    call cpu_time(end_time)
    write(*,'(A,F15.5)') 'Elapsed time (s) : ', end_time - start_time

    ! Now calculate figure of merit.
    max_var = maxval(unc_score)
    fom = 1.0/(max_var**2*(end_time - start_time))
    write(*,'(A,F15.5)') 'Figure of merit (FOM) : ', fom

end program shield_1d