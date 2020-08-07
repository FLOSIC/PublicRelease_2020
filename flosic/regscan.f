C UTEP Electronic Structure Lab (2020)
C  Regularized SCAN 
C  The original code is from SCAN 
C  YY This implementation is as described in our PCCP article
C  [PCCP (2020) doi:10.1039/D0CP02717K].
C  For the original rSCAN implementation, see regscan_jcp.f
C  6/4/2019 YY
C  In certain case, alpha become zero (mostly due to numerical precision
C  issue). When that happens, alpha is set to zero. 
C
C************************ SUBROUTINE VSCANx *****************************
C
C calculates the first order derivatives of Ex wrt n and |grad(n)|
C Sun et. al. PRL (2015)
C
CWritten by Jianwei Sun
C
C everything in Hartree units
C
C ATTANTION: Every values are passed "as they are", i.e. including
C possibly unphysical numerical errors (e.g. negative charge densities)
C values need to be checked accordingly
C***********************************************************************

! YY. Edit the VSCANx functional into spin unpolarized functional
! to make it work better on NRLMOL

       SUBROUTINE RSCANxUNP(RHO,DRHO,TAU_RHO,
     &            EX_metagga,VXD1,VXDD1,AMUXD1)

CInputs
C RU,RD                        density up,down
C DRU, DRD                     abs. val. gradient of density up/down
C DRT                          abs. val. gradient of total density
C TAUU,TAUD                    kinetic energy density up/down
C TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down

C Outputs
C VXD1 VXD2                    THE DERIVATIVES OF EX WRT n 
C VXDD1,VXDD2                  THE DERIVATIVES OF EX WRT |grad n| 
C AMUXD1, AMUXD2                   THE DERIVATIVES OF EX WRT TAU
C EX_metagga                   The exchange energy of metagga

      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER I

      PARAMETER (CFC1X=0.667d0)  !This looks like 2/3 (?)
      PARAMETER (CFC2X=0.800d0)
      PARAMETER (CFD1X=1.240d0)
      PARAMETER (CFK1 =0.065d0)

C other parameters
      PARAMETER (ONE=1.0d0)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)
      PARAMETER (PI=3.1415927d0)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (AX=-0.7385587663820224058842300326808360d0)

      VXD1=0.0d0
      VXDD1=0.0d0
      AMUXD1=0.0d0

!YY. Added these to prevent NAN
      !if(RHO == 0.0d0) RHO = 1d-30
      !if(DRHO == 0.0d0) DRHO = 1d-30

C spin up
C IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD

      !RHO=TWO*RU
      !DRHO=TWO*DRU
      TAUW_RHO=DRHO**TWO/RHO/8.0d0
      !TAU_RHO=TWO*TAUU

C----------------------------------------------------------------------
C construct LDA exchange energy density AND ITS DERIVATIVE WRT n
! YY. EXUNIF LDA Ex per particle and EXLDA is Ex
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      TAU_UNIF=THREE/10.0d0*(THREE*PI**TWO)**THRD2*RHO**THRD5
      ALPHA=(TAU_RHO-TAUW_RHO)/TAU_UNIF
!YY If alpha becomes < 0, it's most likely numerical instability.
! Assume tau - tauw = 0 then.
! For FLOSIC-rSCAN, this is very problematic. We need to make
! an assumption that when ALPAH<0, Fx=1 and d_Fx/d_a=0.
! Not just setting ALPHA=0.
      !if(ALPHA .lt. 0.0d0) ALPHA = 0.0d0

      DPD=-THRD8*P/RHO
!YY Potential DPDD = 0/0 expression
      DPDD=TWO*P/DRHO
!     DPDD=DRHO/(TWO*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      DPDTAU=0.0d0
      DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)
     $        +DRHO**TWO/RHO**(11.0d0/THREE)
     $        *(10.0d0/9.0d0/(THREE*PI**TWO)**THRD2)
      DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)
      DALPHADTAU=ONE/TAU_UNIF




C calculate the exchange enhancement factor and its derivatives wrt p
C and alpha
      call regmetaggafx(cfc1x,cfc2x,cfd1x,cfk1,
     $                p, alpha, Fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2)


C Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=fx1p*DPD+fx1a*DALPHAD
      DFX_DD=fx1p*DPDD+fx1a*DALPHADD
      DFX_DTAU=fx1a*DALPHADTAU

C   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
! YY. orginal lines return Exchange energy for given mesh point
! we want exchange energy per particle instead in NRLMOL.
      !EX_metagga=0.0d0
      !EX_metagga=EX_metagga+EXLDA*FX
      EX_metagga=EXUNIF*FX
      !write(*,'("SCANUNP",5F16.8)') RHO,sqrt(p),Fx,EXUNIF,EX_metagga

      RETURN
      END
      
      

C************************ SUBROUTINE metaggafx *****************************
C
C calculates the first and second order derivatives of the exchange enhancement factor
C wrt its two ingredients p and alpha
C
CWritten by Jianwei Sun
C
C everything is dimentionless
C
C  INPUTS: P ALPHA
C
C  OUTPUTS: D FX /DP, D FX /D ALPHA, D2 FX /D P2, D2 FX /D P D ALPHA, 
C  D2 FX/D ALPHA2
C***********************************************************************

       SUBROUTINE regmetaggafx(cfc1,cfc2,cfd1,cfk1,
     $            p, a, fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2)

CInputs
c cfc1, cfc2, cfd1         fitting parameters
C p                        square of the dimensionless reduced density gradient
C a                        dimensionless deviation from single orbital

C Outputs
C fx                    exchange enhancement factor
C fx1p                  d fx /dp
C fx1a                  d fx /da
c fx2p2                 d2 fx /dp2
c fx2pa                 d2 fx /dpda
c fx2a2                 d2 fx /da2

       IMPLICIT REAL*8 (A-H,O-Z)

       INTEGER I

       PARAMETER (cfk0=0.1740d0)
       PARAMETER (CFA1=4.9479d0)

       PARAMETER (CFMUAK=10.00d0/81.0d0)
       PARAMETER (CFB1=0.156632d0)
       PARAMETER (CFB2=0.12083d0)
       PARAMETER (CFB3=0.5d0)


c      PARAMETER (CFC1=0.9165d0)
c      PARAMETER (CFC2=1.55d0)
c      PARAMETER (CFD1=2.6d0)

C other parameters
       PARAMETER (ONE=1.0d0)
       PARAMETER (TWO=2.0d0)
       PARAMETER (THREE=3.0d0)
       PARAMETER (FOUR=4.0d0)

       real*8 :: poly(1:9)
       poly(1) =  1.00000000D00
       poly(2) = -0.67700000D00
       poly(3) = -0.44455550D00
       poly(4) = -0.6210866010493596D00
       poly(5) =  1.3968970444898194D00
       poly(6) = -0.8591980415966991D00
       poly(7) =  0.22746334147858077D00
       poly(8) = -0.02252024332234138D00
       poly(9) = -0.0D00
! Degree 8 version
!      poly(4) = -0.2715678271666667D00
!      poly(5) =  0.4881482323948366D00
!      poly(6) =  0.03557001954295549D00
!      poly(7) = -0.19755148756274582D00
!      poly(8) =  0.07590424340301653D00
!      poly(9) = -0.008947680611396027D00


       fx=0.0d0
       fx1p=0.0d0
       fx1a=0.0d0
       fx2p2=0.0d0
       fx2pa=0.0d0
       fx2a2=0.0d0


       P2=P*P
       oma=one-a
       oma2=oma*oma

C Enhancement factor Fx
C Hx0
       HX0=ONE+CFK0
C HX1
       cfb4=cfmuak**two/cfk1-0.112654d0
       if (cfb4 .gt. 0.0d0) then
          wfac=cfb4*p2*dexp(-cfb4*p/cfmuak)
          d_wfac_dp=cfb4*p*dexp(-cfb4*p/cfmuak)*(two-cfb4*p/cfmuak)
       else
          wfac=cfb4*p2*dexp(cfb4*p/cfmuak)
          d_wfac_dp=cfb4*p*dexp(cfb4*p/cfmuak)*(two+cfb4*p/cfmuak)
       endif
       vfac=cfb1*p+cfb2*oma*dexp(-cfb3*oma2)
       yfac=cfmuak*p+wfac+vfac**two
       hx1=one+cfk1-cfk1/(one+yfac/cfk1)
C FA
       FA=0.d0
       if(a .lt. 0.0d0) then
          FA=1.0d0  ! fail safe
       else if (a .lt. 2.5d0) then
          FA= poly(1) + poly(2)*a + poly(3)*a**2.0 + poly(4)*a**3.0 
     &      + poly(5)*a**4.0 + poly(6)*a**5.0 + poly(7)*a**6.0 
     &      + poly(8)*a**7.0 + poly(9)*a**8.0
       else
          FA=-cfd1*dExp(cfc2/oma)
       end if
c gx
       p14=p**(one/four)
       gx=1.d0
       if (p .gt. 0.d0) then
           gx=one-dexp(-cfa1/p14)
       endif
c Fx1
       Fx1 = hx1+FA*(hx0-hx1)
C Fx
       Fx = Fx1*gx
C First order derivatives of Fx
C d Hx0 /dp
       d_hx0_dp=0.d0
C d Hx1 /dp
       d_vfac_dp=cfb1
       d_yfac_dp=cfmuak+d_wfac_dp+two*vfac*d_vfac_dp
       d_hx1_dp=d_yfac_dp/(one+yfac/cfk1)**two
C d Hx1 /da
       d_vfac_da=-cfb2*(one-two*cfb3*oma2)*dexp(-cfb3*oma2)
       d_yfac_da=two*vfac*d_vfac_da
       d_hx1_da=d_yfac_da/(one+yfac/cfk1)**two
C d FA /da
       d_FA_da=0.d0
       if(a. lt. 0.0d0) then
          d_FA_da=0.0d0
       else if (a .lt. 2.5d0) then
        d_FA_da = poly(2) + 2.0d0*poly(3)*a + 3.0d0*poly(4)*a**2.0 
     &          + 4.0d0*poly(5)*a**3.0 + 5.0d0*poly(6)*a**4.0 
     &          + 6.0d0*poly(7)*a**5.0 + 7.0d0*poly(8)*a**6.0
     &          + 8.0d0*poly(9)*a**7.0
       else
         d_FA_da=-cfc2*cfd1*dExp(cfc2/oma)/oma**two
       endif
c d gx / dp
       d_gx_dp=0.d0
       if (p .gt. 0.d0) then
           d_gx_dp=-cfa1/four/p/p14*dexp(-cfa1/p14)
       endif
C d Fx1 /dp
       d_Fx1_dp=d_hx1_dp+FA*(d_hx0_dp-d_hx1_dp)
C d Fx1 /da
       d_Fx1_da=(one-FA)*d_hx1_da+d_FA_da*(hx0-hx1)

c d Fx /dp
       fx1p=d_Fx1_dp*gx+Fx1*d_gx_dp
C d Fx /da
       fx1a=d_Fx1_da*gx


C ********************************* Attention second-order derivatives
C are not ready yet**************************************************
C second order derivatives of Fx
C d2 Hx0 /dpp
       d2_hx0_dpp=-two*cfmuak*cfmuak/cfk0
     $            /(one+(cfmuak*p+cfc)/cfk0)**three
C d2 Hx1 /dpp
       d2_hx1_pbe_dpp=-two*cfmuak*cfmuak/cfk0
     $                /(one+cfmuak*p/cfk0)**three
       d2_vpa_dpp=0.0d0
       d2_wpa_dpp=-two*cfb3*(one-two*cfb3*p*p)*wpa
       d2_hx1_dpp=d2_hx1_pbe_dpp+d2_vpa_dpp*wpa
     $            +d_vpa_dp*d_wpa_dp*two+vpa*d2_wpa_dpp
C d2 Hx1 /dpa
       d2_vpa_dpa=-cfb1
       d2_wpa_dpa=-two*cfb3*p*d_wpa_da
       d2_hx1_dpa=d2_vpa_dpa*wpa+d_vpa_dp*d_wpa_da
     $            +d_vpa_da*d_wpa_dp+vpa*d2_wpa_dpa
C d2 Hx1 /daa
       d2_vpa_daa=two*cfb2
       d2_wpa_daa=two*cfb3*(two*cfb3*oma2-one)*wpa
       d2_hx1_daa=d2_vpa_daa*wpa+two*d_vpa_da*d_wpa_da
     $            +vpa*d2_wpa_daa
C d2 FA /daa
       d2_FA_daa=0.d0
       if(a .lt. 0.0d0) then
        d2_FA_daa = 0.0d0
       else if (a .lt. 2.5d0) then
        d2_FA_daa = 2.0d0*poly(3) + 6.0d0*poly(4)*a 
     &          + 12.0d0*poly(5)*a**2.0 + 20.0d0*poly(6)*a**3.0 
     &          + 30.0d0*poly(7)*a**4.0 + 42.0d0*poly(8)*a**5.0
     &          + 56.0d0*poly(9)*a**6.0
       else
        d2_FA_daa=-cfc2*cfd1*(two+cfc2/oma)*dExp(cfc2/oma)/oma**three
       endif
C d2 Fx /dpp
       fx2p2=d2_hx1_dpp+FA*(d2_hx0_dpp-d2_hx1_dpp)
C d2 Fx /dpa
       fx2pa=(one-FA)*d2_hx1_dpa+d_fa_da*(d_hx0_dp-d_hx1_dp)
c d2 Fx /daa
       fx2a2=(one-FA)*d2_hx1_daa+d2_FA_daa*(hx0-hx1)
     $         -two*d_fa_da*d_hx1_da
      RETURN
      END 







C*****************************VSCANc*********************************
C
C
C RU,RD                        density up,down
C DRU, DRD                     abs. val. gradient of density up/down
C DRT                          abs. val. gradient of total density
C TAUU,TAUD                    kinetic energy density up/down
C TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
C VCD1 VCD2                    THE DERIVATIVES OF EC WRT n 
C VCDD1,VCDD2                  THE DERIVATIVES OF EC WRT |grad n| 
C AMUCD1, AMUCD2                   THE DERIVATIVES OF EC WRT TAU
C
C***********************************************************************

      SUBROUTINE RSCANc(
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, 
     &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2)

      IMPLICIT REAL*8 (A-H,O-Z)



C other parameters
      PARAMETER (ONE=1.0d0)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)

      PARAMETER (cfc1=0.640d0)
      PARAMETER (cfc2=1.500d0)
      PARAMETER (cfd1=0.700d0)


      call regmetaggac(cfc1,cfc2,cfd1,
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD,
     &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2)


      RETURN
      END 


C*****************************VrevTPSSc*********************************
C
C
C RU,RD                        density up,down
C DRU, DRD                     abs. val. gradient of density up/down
C DRT                          abs. val. gradient of total density
C TAUU,TAUD                    kinetic energy density up/down
C TAUWU,TAUWD                  Weizsaecker kinetic energy density up/down
C VCD1 VCD2                    THE DERIVATIVES OF EC WRT n 
C VCDD1,VCDD2                  THE DERIVATIVES OF EC WRT |grad n| 
C AMUCD1, AMUCD2                   THE DERIVATIVES OF EC WRT TAU
C
C***********************************************************************

      SUBROUTINE regmetaggac(cfc1,cfc2,cfd1,
     &   RU,RD,DRU,DRD,DRT,TAUU,TAUD, 
     &   Ec_revTPSS,VCD1,VCDD1,VCD2,VCDD2,AMUCD1,AMUCD2)

      IMPLICIT REAL*8 (A-H,O-Z)



C other parameters
      PARAMETER (ONE=1.0d0)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)

      real*8 :: poly(1:9)
      poly(1) =  1.00000000D00
      poly(2) = -0.6400D00
      poly(3) = -0.4352D00
      poly(4) = -1.5356842886420188D00
      poly(5) =  3.061557572493911D00
      poly(6) = -1.9157083913363522D00
      poly(7) =  0.5168839298923616D00
      poly(8) = -0.051848822407902895D00
      poly(9) = -0.0D00
! degree 8 polynomial version
!     poly(4) = -0.27409D00
!     poly(5) = -0.21858757797533326D00
!     poly(6) =  1.3139729875872022D00
!     poly(7) = -1.0172147250963217D00
!     poly(8) =  0.30341612927368594D00
!     poly(9) = -0.03229681378923527D00


! YY. when DRU or DRD is zero you get 0/0 problem.
! To avoid it, assign DRU (DRD) a tiny number
! ---> I added a conditional statement 
      if(DRU == 0.0d0) DRU = 1.0d-100
      if(DRD == 0.0d0) DRD = 1.0d-100 


      VCD1=0.0d0
      VCD2=0.0d0
      VCDD1=0.0d0
      VCDD2=0.0d0
      AMUCD1=0.0d0
      AMUCD2=0.0d0

      RT=RU+RD
      YA=DRU**2.0d0
      YB=DRD**2.0d0
      Y=DRT**2.0d0
C    YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2.0d0
      TAUW=1.0d0/8.0d0*(Y/RT)
      TAU=TAUU+TAUD
      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.99999999999990d0),0.99999999999990d0)
      DS_ZETA=((ONE+ZETA)**THRD5+(ONE-ZETA)**THRD5)/TWO
      DX_ZETA=((ONE+ZETA)**THRD4+(ONE-ZETA)**THRD4)/TWO
      TAU0=0.30d0*(THREE*PISQ)**THRD2*RT**THRD5*DS_ZETA
      ALPHA=(TAU-TAUW)/TAU0
!If alpha becomes < 0, that is likely due to numerical precision. Make
!it zero. -YY
!Update: Need to explicitly set Fc=1 and d_Fc/d_a=0 when alpha < 0.
      !if(ALPHA .lt. 0.0d0) ALPHA = 0.0d0

C     DERIVATIVES OF ZETA, ALPHA, ETC
C     C D1 AND D2
      D_ZETA_D1=TWO*RD/RT**TWO
      D_ZETA_D2=-TWO*RU/RT**TWO
      D_DS_ZETA_D1=THRD5*((ONE+ZETA)**THRD2
     $             -(ONE-ZETA)**THRD2)*D_ZETA_D1/TWO
      D_DS_ZETA_D2=THRD5*((ONE+ZETA)**THRD2
     $             -(ONE-ZETA)**THRD2)*D_ZETA_D2/TWO
      DTAU0D1=THRD5*TAU0/RT+TAU0*D_DS_ZETA_D1/DS_ZETA
      DTAU0D2=THRD5*TAU0/RT+TAU0*D_DS_ZETA_D2/DS_ZETA

      DTAUWD1=-TAUW/RT
      DTAUWD2=DTAUWD1

      D_ALPHA_D1=(-DTAUWD1-ALPHA*DTAU0D1)/TAU0
      D_ALPHA_D2=(-DTAUWD2-ALPHA*DTAU0D2)/TAU0

C     C DD1 AND DD2
      DYCDD1=YC/DRU
      DYCDD2=YC/DRD
!YY. work around for 0/0 
!     if(DRU .eq. 0.0d0) DYCDD1=-DRU/2.0d0
!     if(DRD .eq. 0.0d0) DYCDD2=-DRD/2.0d0 
!     YC/DRU==(Y-YA-YB)/(2.0d0*DRU)
!           ==(DRT**2-DRU**2-DRD**2)/(2*DRU)
!end work around
      DTAUWDD1=ONE/FOUR/RT*(DRU+DYCDD1)
      DTAUWDD2=ONE/FOUR/RT*(DRD+DYCDD2)

      D_ALPHA_DD1=-DTAUWDD1/TAU0
      D_ALPHA_DD2=-DTAUWDD2/TAU0

C     C DTAU
      D_ALPHA_DTAU=ONE/TAU0

C     FCAKE1, FUNKC_ALPHA, GUNKC_ALPHA AND THEIR DERIVATIVES: D_FCAKE1_D1, D_GUNKC_ALPHA_D1 AND
C     D_FUNKC_ALPHA_D1, _D2, _DD1, _DD2, _DTAU
      ALPHA2=ALPHA*ALPHA
      ALPHA3=ALPHA*ALPHA2
      ALPHA4=ALPHA*ALPHA3

C-----------------------------------------------------------------------------------------------------------------------------
C This part needs to be modified when modeling a metaGGA correlation
C-----------------------------------------------------------------------------------------------------------------------------
C   FUNKC_ALPHA, the interpolation function between alpha=0 and 1
      FUNKC_ALPHA=0.0d0
      IF(ALPHA .LT. 0.0d0) THEN
       FUNKC_ALPHA = 1.0d0
      ELSE IF(ALPHA .LT. 2.5d0) THEN
       FUNKC_ALPHA = poly(1) + poly(2)*ALPHA + poly(3)*ALPHA**2.0 
     &             + poly(4)*ALPHA**3.0 + poly(5)*ALPHA**4.0 
     &             + poly(6)*ALPHA**5.0 + poly(7)*ALPHA**6.0 
     &             + poly(8)*ALPHA**7.0 + poly(9)*ALPHA**8.0
      ELSE
          FUNKC_ALPHA=-cfd1*EXP(-CFC2/(ALPHA-ONE))
      ENDIF

C  D_FUNKC_ALPHA_DALPHA
      D_FUNKC_ALPHA_DALPHA=0.0d0
      IF(ALPHA .LT. 0.0d0) THEN
       D_FUNKC_ALPHA_DALPHA=0.0d0
      ELSE IF(ALPHA .LT. 2.5d0) THEN
       D_FUNKC_ALPHA_DALPHA = poly(2) + 2.0d0*poly(3)*ALPHA 
     &       + 3.0d0*poly(4)*ALPHA**2.0 + 4.0d0*poly(5)*ALPHA**3.0
     &       + 5.0d0*poly(6)*ALPHA**4.0 + 6.0d0*poly(7)*ALPHA**5.0
     &       + 7.0d0*poly(8)*ALPHA**6.0 + 8.0d0*poly(9)*ALPHA**7.0
      ELSE
       D_FUNKC_ALPHA_DALPHA=cfc2*FUNKC_ALPHA
     $                         /(ALPHA-ONE)**TWO
      ENDIF

C epsilon_0(rs,zeta,t) and its derivatives
      CALL CORGGA_0(RU,RD,DRU,DRD,DRT,
     $        EPPGGA0,EPPGGA0_D1,EPPGGA0_D2,EPPGGA0_DD1,EPPGGA0_DD2)
C epsilon_1(rs,zeta,t), which is modified from the modPBE by simplifying
C it t-dependence, and its derivatives
      CALL CORGGA_1(RU,RD,DRU,DRD,DRT,
     $       EPPGGA1,EPPGGA1_D1,EPPGGA1_D2,EPPGGA1_DD1,EPPGGA1_DD2)


      EPP=EPPGGA1+FUNKC_ALPHA*(EPPGGA0-EPPGGA1)
      D_EPP_D1=EPPGGA1_D1
     $          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_D1
     $          +FUNKC_ALPHA*(EPPGGA0_D1-EPPGGA1_D1)
      D_EPP_D2=EPPGGA1_D2
     $          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_D2
     $          +FUNKC_ALPHA*(EPPGGA0_D2-EPPGGA1_D2)
      D_EPP_DD1=EPPGGA1_DD1
     $          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_DD1
     $          +FUNKC_ALPHA*(EPPGGA0_DD1-EPPGGA1_DD1)
      D_EPP_DD2=EPPGGA1_DD2
     $          +(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_DD2
     $          +FUNKC_ALPHA*(EPPGGA0_DD2-EPPGGA1_DD2)
      !YY. D_ALPHA_DD2 may have 0/0
      !write(*,*) D_ALPHA_DD2
      D_EPP_DTAU=(EPPGGA0-EPPGGA1)*D_FUNKC_ALPHA_DALPHA*D_ALPHA_DTAU

C-----------------------------------------------------------------------------------------------------------------------------
C End of modification
C-----------------------------------------------------------------------------------------------------------------------------


      VCD1=EPP+RT*D_EPP_D1
      VCD2=EPP+RT*D_EPP_D2
      VCDD1=RT*D_EPP_DD1
      VCDD2=RT*D_EPP_DD2
      AMUCD1=RT*D_EPP_DTAU
      AMUCD2=AMUCD1

! YY. Since we want correlation energy per particle, don't multiply with
! RT(=RU+RD)
      !EC_REVTPSS=RT*EPP
      EC_REVTPSS=EPP

      RETURN
      END 


C######################################################################
