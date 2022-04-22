!***********************************************************************
! Module description
!***********************************************************************
!
! JWF : this subroutine computes the exchange energy of SCAN
!  and some of its regularised variants.
!  Wrapper routines have been created to abuse the spin-polarized
!  routine for computation of the closed-shell results, and then
!  use those same routines for the spin-polarized case.
!
!  Three Integer flags are required to select SCAN variant:
!   IALPHA  :  iso-orbital indicator
!               0: alpha        (scan)
!               1: alpha'       (rscan)
!               2: \bar{alpha}  (r++scan, r2scan, r4scan)
!
!   IINTERP :  interpolation function
!               0: scan
!               1: rscan, r++scan, r2scan, r4scan
!
!   IDELFX  :  gradient expansion correction
!               0: scan (scan, rscan, r++scan)
!               1: 2nd order (r2scan)
!               2: 4th order (r4scan)
!
!
!  Hence, functionals are accessed by passing:
!   SCAN:       (0, 0, 0)
!   rSCAN:      (1, 1, 0)
!   r++SCAN:    (2, 1, 0)
!   r2SCAN:     (2, 1, 1)
!   r4SCAN:     (2, 1, 2)
!
! Interface follows the original SCAN interface
! Author: Jefferson E Bates
! eMail : jeb@temple.edu
! Date  : 23.06.2016
!
! Updated:
! Author: James W. Furness
! eMail : jfurness@tulane.edu (james.w.furness.1@gmail.com)
! Date  : 24/06/2020
!
! This work is made available under the CC0 1.0 Universal (CC0 1.0)
! Public Domain Dedication.
! https://creativecommons.org/publicdomain/zero/1.0/
!
! The person who associated a work with this deed has dedicated the work
! to the public domain by waiving all of his or her rights to the work
! worldwide under copyright law, including all related and neighboring
! rights, to the extent allowed by law.
!
! You can copy, modify, distribute and perform the work, even for
! commercial purposes, all without asking permission. See Other Information
! below.
!
! Other Information:
!
! In no way are the patent or trademark rights of any person affected by CC0,
! nor are the rights that other persons may have in the work or in how the work
! is used, such as publicity or privacy rights.
!
! Unless expressly stated otherwise, the person who associated a work with this
! deed makes no warranties about the work, and disclaims liability for all uses
! of the work, to the fullest extent permitted by applicable law.
!
! When using or citing the work, you should not imply endorsement by the author
! or the affirmer.
!
! While we have made every effort to ensure the code's correctness, it is provided
! as is and no warranties or guarantees are given.
!***********************************************************************

subroutine xscan_r0(d0, g00, t0, f, IALPHA, IINTERP, IDELFX)
! SCAN exchange functional, 0th derivative
    IMPLICIT NONE
!    include 'fortrankinds.h'
    integer, PARAMETER :: kdp = KIND(1.0d0)

    real(kdp), intent(in) :: d0, g00, t0
    real(kdp), intent(out) :: f
    integer, intent(in) :: IALPHA, IINTERP, IDELFX

    real(kdp) :: vxd0 = 0.0d0
    real(kdp) :: vxg0 = 0.0d0
    real(kdp) :: vxt0 = 0.0d0

    ! Input densities are already twice the actual spin 0 density

    call eps_SCANx(d0, sqrt(g00), t0, f, vxd0, vxg0, vxt0, &
   &               IALPHA, IINTERP, IDELFX)

end subroutine xscan_r0

subroutine xscan_r1(d0, g00, t0, vxd0, vxg00, vxt0, f, &
   &                IALPHA, IINTERP, IDELFX)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)
! SCAN exchange functional, 1st derivatives

    real(kdp), intent(in) :: d0, g00, t0
    real(kdp), intent(out) :: f, vxd0, vxg00, vxt0
    integer, intent(in) :: IALPHA, IINTERP, IDELFX

    real(kdp) :: vxg0 = 0.0d0

!    call eps_SCANx(2.0d0*d0, 2.0d0*sqrt(g00), 2.0d0*t0, f, vxd0, vxg0, vxt0, &
!   &               IALPHA, IINTERP, IDELFX)
    call eps_SCANx(d0, g00, t0, f, vxd0, vxg0, vxt0, &
   &               IALPHA, IINTERP, IDELFX)

    vxd0 = vxd0
    !vxg00 = vxg0/(2.0d0*sqrt(g00))
    vxg00 = vxg0/g00 ! Edited where input g00 is 2*sqrt(g00)
    vxt0 = vxt0

end subroutine xscan_r1

subroutine xscan_u0(d0, d1, g00, g11, t0, t1, f, &
   &                IALPHA, IINTERP, IDELFX)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)
! SCAN exchange functional, 1st derivatives

    real(kdp), intent(in) :: d0, d1, g00, g11, t0, t1
    real(kdp), intent(out) :: f
    integer, intent(in) :: IALPHA, IINTERP, IDELFX

    real(kdp) :: vxd0, vxg0, vxt0
    real(kdp) :: den, grd, tau
    real(kdp) :: ex0, ex1

    ! Dummy potential
    vxd0 = 0.0d0
    vxg0 = 0.0d0
    vxt0 = 0.0d0

    den = 2.0d0*d0
    grd = 2.0d0*sqrt(g00)
    tau = 2.0d0*t0

    call eps_SCANx(den, grd, tau, ex0, vxd0, vxg0, vxt0, &
   &               IALPHA, IINTERP, IDELFX)

    den = 2.0d0*d1
    grd = 2.0d0*sqrt(g11)
    tau = 2.0d0*t1

    call eps_SCANx(den, grd, tau, ex1, vxd0, vxg0, vxt0, &
   &               IALPHA, IINTERP, IDELFX)

    f = (ex0 + ex1)/2.0d0

end subroutine xscan_u0

subroutine xscan_u1(d0, d1, g00, g11, t0, t1, vxd0, vxd1, vxg00, vxg11, vxt0, vxt1, f, &
   &                IALPHA, IINTERP, IDELFX)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)
! SCAN exchange functional, 1st derivatives

    real(kdp), intent(in) :: d0, d1, g00, g11, t0, t1
    real(kdp), intent(out) :: f, vxd0, vxd1, vxg00, vxg11, vxt0, vxt1
    integer, intent(in) :: IALPHA, IINTERP, IDELFX

    real(kdp) :: vxg0, vxg1

    real(kdp) :: den, grd, tau
    real(kdp) :: ex0, ex1

    den = 2.0d0*d0
    grd = 2.0d0*sqrt(g00)
    tau = 2.0d0*t0

    call eps_SCANx(den, grd, tau, ex0, vxd0, vxg0, vxt0, &
   &               IALPHA, IINTERP, IDELFX)

    den = 2.0d0*d1
    grd = 2.0d0*sqrt(g11)
    tau = 2.0d0*t1

    call eps_SCANx(den, grd, tau, ex1, vxd1, vxg1, vxt1, &
   &               IALPHA, IINTERP, IDELFX)

    f = (ex0 + ex1)/2.0d0

    vxg00 = vxg0/(2.0d0*sqrt(g00))
    vxg11 = vxg1/(2.0d0*sqrt(g11))
end subroutine xscan_u1

subroutine eps_SCANx(den, grd, tau, eps_x, dedd, dedg, dedt, &
   &                 IALPHA, IINTERP, IDELFX)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)

    real(kdp), intent(in) :: den, grd, tau
    real(kdp), intent(out) :: eps_x, dedd, dedg, dedt
    integer, intent(in) :: IALPHA, IINTERP, IDELFX

    real(kdp) :: cp, den83, den53
    real(kdp) :: exlda, dexldadd, fx, dfxdp, dfxda
    real(kdp) :: p, dpdd, dpdg
    real(kdp) :: alpha, dadd, dadg, dadt, reg, dreg
    real(kdp) :: tueg, dtuegdd, tueg_con
    real(kdp) :: tauw, dtauwdd, dtauwdg

    real(kdp), parameter :: AX = -0.7385587663820224058842300326808360d0
    real(kdp), parameter :: PI =  3.141592653589793238462643383279503d0 !3.1415927d0
    real(kdp), parameter :: PI2 = PI*PI
    real(kdp), parameter :: ETA = 1.0d-3
    real(kdp), parameter :: TAU_R = 1.0d-4
    real(kdp), parameter :: A_REG = 1.0d-3

!       Reduced density gradient [FD]
    den83 = den**(8.0d0/3.0d0)
    cp = 4.0d0*(3.0d0*PI2)**(2.0d0/3.0d0)
    p = grd**2/(cp*den83)
    dpdd = -8.0d0/3.0d0*p/den
    dpdg = 2.0d0*p/grd

!       Regularised Alpha [FD]
    tueg_con = 3.0d0/10.0d0*(3.0d0*PI2)**(2.0d0/3.0d0)
    den53 = den**(5.0d0/3.0d0)
    if (IALPHA .eq. 1) then ! regularised tau_ueg for rSCAN \tilde{alpha}
        tueg = tueg_con*den53 + TAU_R
    else                    ! Unregularised tau_ueg
        tueg = tueg_con*den53
    end if
    !dtuegdd = 5.0d0/3.0d0*tueg/den
    dtuegdd = 5.0d0/3.0d0*tueg_con*den53/den

    tauw = grd**2/8.0d0/den
    dtauwdd = -grd**2/(8.0d0*den**2)
    dtauwdg = 2.0d0*grd/(8.0d0*den)

    if (IALPHA .eq. 0) then ! SCAN unregularised alpha
        alpha = (tau - tauw)/tueg
        dadd = -(tau - tauw)*dtuegdd/tueg**2 - dtauwdd/tueg
        dadg = -dtauwdg/tueg
        dadt = 1.0d0/tueg

    else if (IALPHA .eq. 1) then ! rSCAN alpha'
        alpha = (tau - tauw)/tueg
        dadd = ((tauw - tau)*dtuegdd - tueg*dtauwdd)/tueg**2
        dadg = -grd/(4.0d0*den*tueg)
        dadt = 1.0d0/tueg

        if (A_REG .gt. 0.0d0) then
            reg = alpha**3/(alpha**2 + A_REG)
            dreg = (alpha**4 + 3.0d0*alpha**2*A_REG)/(alpha**2 + A_REG)**2
            dadd = dadd*dreg
            dadg = dadg*dreg
            dadt = dadt*dreg
            alpha = reg
        endif

    else if (IALPHA .eq. 2) then ! \bar{alpha} for r2scan and r4scan
        alpha = (tau - tauw)/(tueg + ETA*tauw)
        dadd = -dtauwdd/(tueg + ETA*tauw) &
       &    - (tau - tauw)*(dtuegdd + ETA*dtauwdd)/(tueg + ETA*tauw)**2
        dadg = -ETA*(tau - tauw)*dtauwdg/(tueg + ETA*tauw)**2 &
       &    - dtauwdg/(tueg + ETA*tauw)
        dadt = 1.0d0/(tueg + ETA*tauw)

    else
        !call quit('ERROR: Unknown IALPHA in SCAN')
        print *,'ERROR: Unknown IALPHA in SCAN'
        call stopit
    end if

!       UEG exchange density [FD]
    !exlda = AX*den**(4.0d0/3.0d0)
    exlda = AX*den**(1.0d0/3.0d0) !YY Ex per particle
    dexldadd = 4.0d0*AX*den**(1.0d0/3.0d0)/3.0d0

    call exchange_enhancement(p, alpha, fx, dfxdp, dfxda, ETA, IINTERP, IDELFX)

!       Exchange Energy[FD]
    eps_x = exlda*fx
    !YY multiplying density back to exlda
    exlda=exlda*den ! This is needed for derivatives
    dedd = dexldadd*fx + exlda*(dfxdp*dpdd + dfxda*dadd)
    dedg = exlda*(dfxdp*dpdg + dfxda*dadg)
    dedt = exlda*dfxda*dadt

end subroutine

subroutine exchange_enhancement(p, alpha, Fx, dfxdp, dfxda, ETA, IINTERP, IDELFX)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)

    real(kdp), intent(in) :: p, alpha
    real(kdp), intent(out) :: Fx, dfxdp, dfxda
    real(kdp), intent(in) :: ETA
    integer, intent(in) :: IINTERP, IDELFX


    real(kdp) :: ief, diefda, oma
    real(kdp) :: h0x
    real(kdp) :: h1x, dh1xdp, dh1xda
    real(kdp) :: gx, dgxdp
    real(kdp) :: del_f2, C2, ALPHA_GE, damp, ddampdp
    real(kdp) :: del_fx, ddel_fxdp, ddel_fxda
    real(kdp) :: wfac, dwfacdp, vfac, dvfacdp, dvfacda
    real(kdp) :: yfac, dyfacdp, dyfacda

    real(kdp), parameter :: A1 = 4.9479d0
    real(kdp), parameter :: K0 = 0.174d0
    real(kdp), parameter :: MU = 10.0d0/81.0d0
    real(kdp), parameter :: K1 = 0.065d0

    real(kdp), parameter :: cfx1 = 0.667d0
    real(kdp), parameter :: cfx2 = 0.800d0
    real(kdp), parameter :: cfdx1 = 1.24d0

    real(kdp), parameter :: B1 = 0.156632d0
    real(kdp), parameter :: B2 = 0.12083d0
    real(kdp), parameter :: B3 = 0.5d0
    real(kdp), parameter :: B4 = MU*MU/K1 - 0.112654d0

    real(kdp), parameter, dimension(8) :: PARAMS = (/ &
    &      -0.023185843322d0,0.234528941479d0,-0.887998041597d0, &
    &      1.451297044490d0,-0.663086601049d0,-0.4445555d0,-0.667d0, &
    &      1.0d0/)
    integer, dimension(8), parameter :: f_x_e = (/7,6,5,4,3,2,1,0/)
    real(kdp), parameter :: D_DAMP2 = 0.361d0

    ALPHA_GE = 20.0d0/27.0d0 + ETA*5.0d0/3.0d0

    ief = 0.0d0
    diefda = 0.0d0
    oma = 1.0d0 - alpha
    if (IINTERP .eq. 0) then  ! scan interpolation function
        if (alpha .lt. 1.0d0) then
            ief = exp(-cfx1*alpha/oma)
            diefda = -cfx1*exp(-cfx1*alpha/oma)/oma**2
        else
            ief = -cfdx1*exp(cfx2/oma)
            diefda = -cfx2*cfdx1*exp(cfx2/oma)/oma**2
        endif

    else if (IINTERP .eq. 1) then ! rscan, r2scan, r4scan interpolation function
        if (alpha .lt. 1.0d-13) then
            ief = exp(-cfx1*alpha/oma)
            diefda = -cfx1*exp(-cfx1*alpha/oma)/oma**2
        else if( alpha .lt. 2.5d0) then
            ief = dot_product(alpha**f_x_e, PARAMS)
            diefda = dot_product(alpha**f_x_e(2:), f_x_e(:7)*PARAMS(:7))
        else if (alpha .ge. 2.5d0) then
            ief = -cfdx1*exp(cfx2/oma)
            diefda = -cfx2*cfdx1*exp(cfx2/oma)/oma**2
        endif

    else
        !call quit('ERROR: Unknown IINERP in SCAN')
        print *,'ERROR: Unknown IINERP in SCAN'
        call stopit
    end if

!       Single orbital enhancement
    h0x = 1.0d0 + K0

!       Slowly varying enhancement
    if (IDELFX .eq. 0) then     ! scan, rscan
        wfac = B4*p**2*exp(-B4*p/MU)
        dwfacdp = B4*p*exp(-B4*p/MU)*(2.0d0 - B4*p/MU)

        vfac = B1*p + B2*oma*exp(-B3*oma**2)
        yfac = MU*p + wfac + vfac**2
        h1x = 1.0d0 + K1 - K1/(1.0d0 + yfac/K1)

        dvfacdp = B1
        dyfacdp = MU + dwfacdp + 2.0d0*vfac*dvfacdp
        dh1xdp = dyfacdp/(1.0d0 + yfac/K1)**2

        dvfacda = -B2*(1.0d0 - 2.0d0*B3*oma**2)*exp(-B3*oma**2)
        dyfacda = 2.0d0*vfac*dvfacda
        dh1xda = dyfacda/(1.0d0 + yfac/K1)**2

    else if (IDELFX .eq. 1 .or. IDELFX .eq. 2) then  ! Second order corrections
        del_f2 = dot_product(f_x_e(:7), PARAMS(:7))
        C2 = -del_f2*(1.0d0 - h0x)

!       Damping
        damp = exp(-p**2/D_DAMP2**4)
        ddampdp = -2.0d0*damp*p/D_DAMP2**4

!       Slowly varying contribution [FD]
        h1x = 1.0d0 + K1 - K1/(1.0d0 + p*(MU + ALPHA_GE*C2*damp)/K1)
        dh1xdp = K1**2*(MU + ALPHA_GE*C2*(damp + p*ddampdp)) &
        &              /(K1 + MU*p + ALPHA_GE*C2*p*damp)**2
        dh1xda = 0.0d0

    else
        !call quit('ERROR: Unknown IDELFX in SCAN')
        print *,'ERROR: Unknown IDELFX in SCAN'
        call stopit
    end if

    gx = 1.0d0 - exp(-A1/p**(1.0d0/4.0d0))
    dgxdp = -A1*exp(-A1/p**(1.0d0/4.0d0))/(4.0d0*p**(5.0d0/4.0d0))

    if (IDELFX .eq. 2) then  ! 4th order corrections for r4scan
        call get_del_fx(p, alpha, del_fx, ddel_fxdp, ddel_fxda, &
        &                        K0, K1, C2, ETA, ALPHA_GE, PARAMS, f_x_e)
    else
        del_fx = 0.0d0
        ddel_fxdp = 0.0d0
        ddel_fxda = 0.0d0
    end if

    fx = (h1x + ief*(h0x - h1x) + del_fx)*gx
    dfxdp = (del_fx + h1x + (h0x - h1x)*ief)*dgxdp &
    &       + gx*(ddel_fxdp + dh1xdp - ief*dh1xdp)
    dfxda = gx*((h0x - h1x)*diefda + ddel_fxda + dh1xda - ief*dh1xda)

end subroutine

subroutine get_del_fx(p, alpha, del_fx, ddel_fxdp, ddel_fxda, &
&                        K0, K1, C2, ETA, ALPHA_GE, PARAMS, f_x_e)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)

    real(kdp), intent(in) :: p, alpha, K0, K1, C2, ETA, ALPHA_GE
    real(kdp), dimension(8), intent(in) :: PARAMS
    integer, dimension(8), intent(in) :: f_x_e
    real(kdp), intent(out) :: del_fx, ddel_fxdp, ddel_fxda
    real(kdp) :: order_1, dorder_1dp, dorder_1da, C_pa, C_aa, C_pp
    real(kdp) :: oma, damp, t1, dt1dp, dt1da, ddampdp, ddampda

    call get_dx_terms(C_aa, C_pa, C_pp, &
    &              K0, K1, C2, ETA, PARAMS, f_x_e)

    oma = 1.0d0 - alpha

    order_1 = C2*(oma - ALPHA_GE*p)
    dorder_1dp = -C2*ALPHA_GE
    dorder_1da = -C2

!       Correcting contribution [FD]
    t1 = order_1 + C_aa*oma**2 + C_pa*p*oma + C_pp*p**2
    dt1dp = dorder_1dp + 2*C_pp*p + C_pa*oma
    dt1da = dorder_1da - 2*C_aa*oma - C_pa*p

    call get_fourth_order_damp(p, alpha, damp, ddampdp, ddampda)

    del_fx = t1*damp
    ddel_fxdp = damp*dt1dp + t1*ddampdp
    ddel_fxda = t1*ddampda + dt1da*damp

end subroutine get_del_fx

subroutine get_dx_terms(C_aa, C_pa, C_pp, K0, K1, C2, ETA, &
&                  PARAMS, f_x_e)
    IMPLICIT NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)

    real(kdp), intent(in) :: K0, K1, C2, ETA
    real(kdp), dimension(8), intent(in) :: PARAMS
    integer, dimension(8), intent(in) :: f_x_e
    real(kdp), intent(out) :: C_aa, C_pa, C_pp
    real(kdp) :: h0x, eta_term, del_f2, del_f4
    real(kdp) :: ALPHA_GE

    real(kdp), parameter :: MU = 10.0d0/81.0d0

    ALPHA_GE = 20.0d0/27.0d0 + ETA*5.0d0/3.0d0

    eta_term = ETA*3.0d0/4.0d0 + 2.0d0/3.0d0
    h0x = 1.0d0 + K0

    del_f2 = dot_product(f_x_e(:7), PARAMS(:7))
    del_f4 = dot_product(f_x_e(:7)*(f_x_e(:7) - 1.0d0), PARAMS(:7))

    C_aa = 73.0d0/5000.0d0 - 0.5d0*del_f4*(h0x - 1.0d0)

    C_pa = 511.0d0/13500.0d0 - 73.0d0/1500.0d0*ETA &
    &            - del_f2*(ALPHA_GE*C2 + MU)

    C_pp = 146.0d0/2025.0d0*eta_term**2 - 73.0d0/405.0d0*eta_term &
    &            + (ALPHA_GE*C2 + MU)**2/K1

end subroutine

subroutine get_fourth_order_damp(p, alpha, damp4, &
&            ddamp4dp, ddamp4da)
    implicit NONE
!    include 'fortrankinds.h'
integer, PARAMETER :: kdp = KIND(1.0d0)

    real(kdp), intent(in) :: p, alpha
    real(kdp), intent(out) :: damp4, ddamp4dp, ddamp4da
    real(kdp) :: t1, t2, dt1da, dt2dp, dt2da, oma
    real(kdp), parameter :: DX_DAMP4_P = 0.232d0
    real(kdp), parameter :: DX_DAMP4_A = 0.232d0

    oma = 1.0d0 - alpha

    t1 = 2.0d0*alpha**2/(1.0d0 + alpha**4)
    dt1da = -4.0d0*alpha*(alpha**4 - 1.0d0)/(1.0d0 + alpha**4)**2

    t2 = exp(-oma**2/DX_DAMP4_A**2 - p**2/DX_DAMP4_P**4)
    dt2dp = -2.d0*t2*p/DX_DAMP4_P**4
    dt2da = 2*oma*t2/DX_DAMP4_A**2

    damp4 = t1*t2
    ddamp4dp = t1*dt2dp
    ddamp4da = t1*dt2da + t2*dt1da

end subroutine
