module mod_scatt
    !
    ! This module handles the calculation of the new direction of a particle that suffers 
    ! a scattering process.
    !
    use iso_fortran_env, only: int32, real64

    use mod_rng
implicit none

contains

    real(kind=real64) function scatt(uin)
        !
        ! This function samples the scattering angle and determines the new direction 
        ! of a particle that suffers a scattering process. Isotropic scattering.
        !
        real(kind=real64), intent(in) :: uin ! direction of the particle before interaction
        real(kind=real64) :: u0, phi0        ! scattering angles
        real(kind=real64) :: rnno, pi

        pi = 4.0_real64*atan(1.0_real64)

        ! First sample polar angle (as directional cosine).
        rnno = rng_set()
        u0 = 2.0*rnno-1.0

        ! Then sample azimuthal angle.
        rnno = rng_set()
        phi0 = 2.0*pi*rnno

        ! Finally calculate the new direction of the particle.
        scatt = uin*u0 + sqrt(1.0-uin**2)*sqrt(1.0-u0**2)*cos(phi0)

    end function scatt

end module mod_scatt