program shield_1d
    !
    ! Program that simulates the particle transport through a 1D shield using Monte Carlo methods.
    ! Multi-region version. Parallel version using MPI.
    !
    ! Author: Edgardo Doerner, edoerner@fis.puc.cl
    !

    use mpi

    use iso_fortran_env, only: int32, real64

    use mod_rng
    use mod_mfp
    use mod_scatt

implicit none

    ! Geometry parameters
    integer(kind=int32), parameter :: nreg = 1              ! number of region in the 1D shield
    real(kind=real64), dimension(nreg) :: xthick = (/10.0/)  ! thickness of the 1D shield (cm)
    real(kind=real64), dimension(nreg+1) :: xbounds = 0.0    ! boundaries of the shield regions (cm)

    ! Transport parameters
    real(kind=real64), dimension(nreg) :: sigma_t = (/1.0/)  ! Total interaction cross section (cm-1)
    real(kind=real64), dimension(nreg) :: sigma_a = (/0.8/)  ! Absorption cross section (cm-1)
    real(kind=real64), dimension(nreg) :: sigma_s = (/0.2/)  ! Scattering cross section (cm-1)

    ! Particle parameters
    real(kind=real64), parameter :: xin = 0.0    ! initial position (cm)
    real(kind=real64), parameter :: uin = 1.0    ! initial direction
    real(kind=real64), parameter :: wtin = 1.0   ! statistical weight
    integer(kind=int32), parameter :: irin = 1  ! initial region

    integer(kind=int32) :: ir
    real(kind=real64) :: x, u, wt

    ! Simulation parameters
    integer(kind=int32) :: nperbatch = 1E7    
    integer(kind=int32) :: nbatch = 10          ! number of statistical batches

    ! Scoring variables
    real(kind=real64), dimension(0:nreg+1) :: score, score2 = 0.0    ! score(0) : reflection
                                                        ! score(1:nreg) : absorption
                                                        ! score(nreg+1) : transmission

    real(kind=real64), dimension(0:nreg+1) :: mean_score = 0.0   ! mean value of scoring array
    real(kind=real64), dimension(0:nreg+1) :: unc_score = 0.0    ! uncertainty values of scoring array

    integer(kind=int32) :: i, ihist, ibatch     ! loop counters 
    logical :: pdisc                            ! flag to discard a particle
    integer(kind=int32) :: irnew                ! index of new region 
    real(kind=real64) :: pstep                  ! distance to next interaction 
    real(kind=real64) :: dist                   ! distance to boundary along particle direction
    real(kind=real64) :: rnno 
    real(kind=real64) :: max_var, fom
    real(kind=real64) :: start_time, end_time, r_end_time

    ! Definitions needed for MPI execution.
    integer(kind=int32) :: ierr
    integer(kind=int32) :: mpi_rank, mpi_size

    real(kind=real64), dimension(0:nreg+1) :: r_mean_score = 0.0   ! mean value of scoring array
    real(kind=real64), dimension(0:nreg+1) :: r_unc_score = 0.0    ! uncertainty values of scoring array

    call mpi_init(ierr)

    start_time = mpi_wtime()

    ! Identify MPI process and MPI size.
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, ierr)

    ! Initialize RNG
    call rng_init(20180815 + mpi_rank)

    ! Initialize simulation geometry.
    do i = 1,nreg
        xbounds(i+1) = xbounds(i) + xthick(i)
    enddo

    ! Print some info to console
    if (mpi_rank .eq. 0) then
        write(6,'(A)') '*****************************************************************'
        write(6,'(A)') 'Monte Carlo simulation of 1D shield with MPI support '
        write(6,'(A)') '*****************************************************************'
        write(6,'(A, I3)') 'Number of MPI processes : ', mpi_size
        write(*,'(A, I15)') 'Number of histories : ', nperbatch*nbatch
        write(*,'(A, I15)') 'Number of batches : ', nbatch
        write(*,'(A, I15)') 'Number of histories per batch : ', nperbatch
        write(6,'(A)') '*****************************************************************'

        do i = 1,nreg
            write(*,'(A, I5, A)') 'For region ', i, ':'
            write(*,'(A, F15.5)') 'Total interaction cross-section (cm-1): ', sigma_t(i)
            write(*,'(A, F15.5)') 'Absoption interaction cross-section (cm-1): ', sigma_a(i)
            write(*,'(A, F15.5)') 'Scattering interaction cross-section (cm-1): ', sigma_s(i)
            write(*,'(A, F15.5)') '1D shield thickness (cm): ', xthick(i)
        enddo
        write(6,'(A)') '*****************************************************************'
    endif
    
    ! Adjust nperbatch according to number of MPI processes.
    nperbatch = nperbatch/mpi_size

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
                    score(ir) = score(ir) + wt
                    exit
                endif

                ! Determine interaction type.
                rnno = rng_set()
                if(rnno .le. sigma_a(ir)/sigma_t(ir)) then
                    ! The particle was absorbed.
                    score(ir) = score(ir) + wt
                    exit
                else
                    ! The particle was scattered, get new direction to re-enter transport loop.
                    u = scatt(u)
                endif
            
            enddo particle_loop
        enddo ihist_loop

        ! Accumulate results and reset score arrays.
        mean_score = mean_score + score
        unc_score = unc_score + score**2

        score = 0.0
        score2 = 0.0

    enddo ibatch_loop

    ! Combine results across MPI processes.
    call mpi_reduce(mean_score, r_mean_score, size(mean_score), MPI_DOUBLE_PRECISION, MPI_SUM, &
        0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(unc_score, r_unc_score, size(unc_score), MPI_DOUBLE_PRECISION, MPI_SUM, &
        0, MPI_COMM_WORLD, ierr)

    ! Statistical procesing. Note that now we hace a total of nbatch*mpi_size 
    ! statistical batches across MPI processes.
    if(mpi_rank .eq. 0) then
        r_mean_score = r_mean_score/(nbatch*mpi_size)
        r_unc_score = (r_unc_score - nbatch*mpi_size*r_mean_score**2)/(nbatch*mpi_size-1)
        r_unc_score = r_unc_score/(nbatch*mpi_size)
        r_unc_score = sqrt(r_unc_score)

        ! Calculate relative uncertainty, to be used for the FOM.
        r_unc_score = r_unc_score/r_mean_score

        ! Print results to console.
        write(*,'(A,F10.5,A,F10.5,A)') 'Reflection : ', r_mean_score(0)/nperbatch, ' +/-', 100.0*r_unc_score(0), '%'

        ! To calculate the uncertainty of the absorption we need to combine the uncertainty of the deposition 
        ! in each region.
        write(*,'(A,F10.5,A,F10.5,A)') 'Absorption : ', sum(r_mean_score(1:nreg))/nperbatch, ' +/-', & 
            100.0*sum(r_unc_score(1:nreg)), '%'
    
        write(*,'(A,F10.5,A,F10.5,A)') 'Transmission : ', r_mean_score(nreg+1)/nperbatch, ' +/-', &
            100.0*r_unc_score(nreg+1), '%'
    endif

    ! Get elapsed time
    end_time = mpi_wtime()
    call mpi_reduce(end_time, r_end_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
        0, MPI_COMM_WORLD, ierr)
    
    if(mpi_rank .eq. 0) then
        write(*,'(A,F15.5)') 'Elapsed time (s) : ', r_end_time - start_time

        ! Now calculate figure of merit.
        max_var = maxval(r_unc_score)
        fom = 1.0/(max_var**2*(end_time - start_time))
        write(*,'(A,F15.5)') 'Figure of merit (FOM) : ', fom
    endif    

    call mpi_finalize(ierr)

end program shield_1d