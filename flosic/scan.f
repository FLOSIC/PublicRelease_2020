C UTEP Electronic Structure Lab (2020)
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

       SUBROUTINE VSCANx(RU,RD,DRU,DRD,DRT,TAUU,TAUD,
     $               EX_metagga,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

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

      PARAMETER (CFC1X=0.667d0)
      PARAMETER (CFC2X=0.800d0)
      PARAMETER (CFD1X=1.240d0)
      PARAMETER (CFK1 =0.065d0)

      !PARAMETER (CFC1X=0.667_q)
      !PARAMETER (CFC2X=0.800_q)
      !PARAMETER (CFD1X=1.240_q)
      !PARAMETER (CFK1=0.065_q)
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
      VXD2=0.0d0
      VXDD1=0.0d0
      VXDD2=0.0d0
      AMUXD1=0.0d0
      AMUXD2=0.0d0

!YY. Added these to prevent NAN
      !if(RU == 0.0d0) RU = 1d-30
      !if(RD == 0.0d0) RD = 1d-30
      !if(DRU == 0.0d0) DRU = 1d-30
      !if(DRD == 0.0d0) DRD = 1d-30

C spin up
C IN EXD1(2*RU), TAUWU AND TAUU SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RU
      DRHO=TWO*DRU
      TAUW_RHO=DRHO**TWO/RHO/8.0d0
      TAU_RHO=TWO*TAUU

C----------------------------------------------------------------------
C construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      TAU_UNIF=THREE/10.0d0*(THREE*PI**TWO)**THRD2*RHO**THRD5
      ALPHA=(TAU_RHO-TAUW_RHO)/TAU_UNIF

      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0.0d0
      DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)
     $        +DRHO**TWO/RHO**(11.0d0/THREE)
     $        *(10.0d0/9.0d0/(THREE*PI**TWO)**THRD2)
      DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)
      DALPHADTAU=ONE/TAU_UNIF




C calculate the exchange enhancement factor and its derivatives wrt p
C and alpha
       call metaggafx(cfc1x,cfc2x,cfd1x,cfk1,
     $                p, alpha, Fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2)

C Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=fx1p*DPD+fx1a*DALPHAD
      DFX_DD=fx1p*DPDD+fx1a*DALPHADD
      DFX_DTAU=fx1a*DALPHADTAU

C   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
      EX_metagga=0.0d0
      EX_metagga=EX_metagga+EXLDA*FX

C spin down
C IN EXD1(2*RD), TAUWD AND TAUD SCALES AS 2 AND TAU0 SCALES AS 2**FTHRD
      RHO=TWO*RD
      DRHO=TWO*DRD
      TAUW_RHO=DRHO**TWO/RHO/8.0d0
      TAU_RHO=TWO*TAUD

C----------------------------------------------------------------------
C construct LDA exchange energy density AND ITS DERIVATIVE WRT n
      EXUNIF=AX*RHO**THRD
      EXLDA=EXUNIF*RHO
      EXDLDA=EXUNIF*THRD4


      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      !write(*,*) 'SCAN',P,DRHO,RHO
      TAU_UNIF=THREE/10.0d0*(THREE*PI**TWO)**THRD2*RHO**THRD5
      ALPHA=(TAU_RHO-TAUW_RHO)/TAU_UNIF

      DPD=-THRD8*P/RHO
      DPDD=TWO*P/DRHO
      DPDTAU=0.0d0
      DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)
     $        +DRHO**TWO/RHO**(11.0d0/THREE)
     $        *(10.0d0/9.0d0/(THREE*PI**TWO)**THRD2)
      DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)
      DALPHADTAU=ONE/TAU_UNIF

C calculate the exchange enhancement factor and its derivatives wrt p
C and alpha
       call metaggafx(cfc1x,cfc2x,cfd1x,cfk1,
     $                p, alpha, Fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2)

C Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=fx1p*DPD+fx1a*DALPHAD
      DFX_DD=fx1p*DPDD+fx1a*DALPHADD
      DFX_DTAU=fx1a*DALPHADTAU


C   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD2=EXDLDA*FX+EXLDA*DFX_D
      VXDD2=EXLDA*DFX_DD
      AMUXD2=EXLDA*DFX_DTAU

      EX_metagga=EX_metagga+EXLDA*FX
      EX_metagga=EX_metagga/TWO
      !write(*,*) 'SCAND', EX_metagga, EXLDA, FX
      RETURN
      END SUBROUTINE VSCANx

! YY. Edit the VSCANx functional into spin unpolarized functional
! to make it work better on NRLMOL

       SUBROUTINE VSCANxUNP(RHO,DRHO,TAU_RHO,
     &            EX_metagga,VXD1,VXDD1,AMUXD1)
!       SUBROUTINE VSCANx(RU,RD,DRU,DRD,DRT,TAUU,TAUD,
!     $               EX_metagga,VXD1,VXDD1,VXD2,VXDD2,AMUXD1,AMUXD2)

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

      !PARAMETER (CFC1X=0.667_q)
      !PARAMETER (CFC2X=0.800_q)
      !PARAMETER (CFD1X=1.240_q)
      !PARAMETER (CFK1=0.065_q)
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
      !VXD2=0.0d0
      VXDD1=0.0d0
      !VXDD2=0.0d0
      AMUXD1=0.0d0
      !AMUXD2=0.0d0

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

      DPD=-THRD8*P/RHO
!YY Potential DPDD = 0/0 expression
      DPDD=TWO*P/DRHO
!      DPDD=DRHO/(TWO*(THREE*PISQ)**THRD2*(RHO)**THRD8)
      DPDTAU=0.0d0
      DALPHAD=-TAU_RHO*(THREE*PI**TWO*RHO)**THRD2/(TWO*TAU_UNIF**TWO)
     $        +DRHO**TWO/RHO**(11.0d0/THREE)
     $        *(10.0d0/9.0d0/(THREE*PI**TWO)**THRD2)
      DALPHADD=-DRHO/(FOUR*RHO*TAU_UNIF)
      DALPHADTAU=ONE/TAU_UNIF




C calculate the exchange enhancement factor and its derivatives wrt p
C and alpha
       call metaggafx(cfc1x,cfc2x,cfd1x,cfk1,
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

       SUBROUTINE metaggafx(cfc1,cfc2,cfd1,cfk1,
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
          if (a .lt. one) then
             FA=dExp(-cfc1*a/oma)
          endif
          if (a .gt. one) then
             FA=-cfd1*dExp(cfc2/oma)
          endif
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
          if (a .lt. one) then
            d_FA_da=-cfc1*dExp(-cfc1*a/oma)/oma**two
          endif
          if (a .gt. one) then
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
     $                 /(one+(cfmuak*p+cfc)/cfk0)**three
C d2 Hx1 /dpp
            d2_hx1_pbe_dpp=-two*cfmuak*cfmuak/cfk0
     $                     /(one+cfmuak*p/cfk0)**three
            d2_vpa_dpp=0.0d0
            d2_wpa_dpp=-two*cfb3*(one-two*cfb3*p*p)*wpa
            d2_hx1_dpp=d2_hx1_pbe_dpp+d2_vpa_dpp*wpa
     $                 +d_vpa_dp*d_wpa_dp*two+vpa*d2_wpa_dpp
C d2 Hx1 /dpa
            d2_vpa_dpa=-cfb1
            d2_wpa_dpa=-two*cfb3*p*d_wpa_da
            d2_hx1_dpa=d2_vpa_dpa*wpa+d_vpa_dp*d_wpa_da
     $                 +d_vpa_da*d_wpa_dp+vpa*d2_wpa_dpa
C d2 Hx1 /daa
            d2_vpa_daa=two*cfb2
            d2_wpa_daa=two*cfb3*(two*cfb3*oma2-one)*wpa
            d2_hx1_daa=d2_vpa_daa*wpa+two*d_vpa_da*d_wpa_da
     $                 +vpa*d2_wpa_daa
C d2 FA /daa
          d2_FA_daa=0.d0
          if (a .lt. one) then
          d2_FA_daa=cfc1*(cfc1/oma-two)*dExp(-cfc1*a/oma)/oma**three
          endif
       if (a .gt. one) then
       d2_FA_daa=-cfc2*cfd1*(two+cfc2/oma)*dExp(cfc2/oma)/oma**three
       endif
C d2 Fx /dpp
            fx2p2=d2_hx1_dpp+FA*(d2_hx0_dpp-d2_hx1_dpp)
C d2 Fx /dpa
            fx2pa=(one-FA)*d2_hx1_dpa+d_fa_da*(d_hx0_dp-d_hx1_dp)
c d2 Fx /daa
            fx2a2=(one-FA)*d2_hx1_daa+d2_FA_daa*(hx0-hx1)
     $              -two*d_fa_da*d_hx1_da
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

      SUBROUTINE VSCANc(
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


      call Vmetaggac(cfc1,cfc2,cfd1,
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

      SUBROUTINE Vmetaggac(cfc1,cfc2,cfd1,
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
!      if(DRU .eq. 0.0d0) DYCDD1=-DRU/2.0d0
!      if(DRD .eq. 0.0d0) DYCDD2=-DRD/2.0d0 
!      YC/DRU==(Y-YA-YB)/(2.0d0*DRU)
!            ==(DRT**2-DRU**2-DRD**2)/(2*DRU)
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
      IF (ALPHA .LT. ONE) THEN
          FUNKC_ALPHA=EXP(CFC1*ALPHA/(ALPHA-ONE))
       ENDIF
      IF (ALPHA .GT. ONE) THEN
          FUNKC_ALPHA=-cfd1*EXP(-CFC2/(ALPHA-ONE))
       ENDIF

C  D_FUNKC_ALPHA_DALPHA
      D_FUNKC_ALPHA_DALPHA=0.0d0
      IF (ALPHA .lt. one) THEN
          D_FUNKC_ALPHA_DALPHA=-cfc1*FUNKC_ALPHA
     $                         /(ALPHA-ONE)**TWO
       ENDIF
      IF (ALPHA .gt. one) THEN
          D_FUNKC_ALPHA_DALPHA=cfc2*FUNKC_ALPHA
     $                         /(ALPHA-ONE)**TWO
      endif

C epsilon_0(rs,zeta,t) and its derivatives
      CALL CORGGA_0(RU,RD,DRU,DRD,DRT,
     $        EPPGGA0,EPPGGA0_D1,EPPGGA0_D2,EPPGGA0_DD1,EPPGGA0_DD2)
C epsilon_1(rs,zeta,t), which is modified from the modPBE by simplifying
C it t-dependence, and its derivatives
      CALL CORGGA_1(RU,RD,DRU,DRD,DRT,
     $       EPPGGA1,EPPGGA1_D1,EPPGGA1_D2,EPPGGA1_DD1,EPPGGA1_DD2)

      !YY. debug
      !write(*,*) 'SCAN',EPPGGA0_DD2,EPPGGA1_DD2

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
C----------------------------------------------------------------------
      SUBROUTINE CORPBE_revtpss(RS,ZET,EC,VCUP,VCDN,g,sk, 
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
C----------------------------------------------------------------------
C  Official PBE correlation code. K. Burke, May 14, 1996.
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : lgga=flag to do gga (0=>LSD only)
C       : lmetagga=flag to do metagga (revTPSS)
C       : lpot=flag to do potential (0=>energy only)
C  output: ec=lsd correlation energy from [a]
C        : vcup=lsd up correlation potential
C        : vcdn=lsd dn correlation potential
C        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
C        : dvcup=nonlocal correction to vcup
C        : dvcdn=nonlocal correction to vcdn
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C References:
C [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
C     {\sl Generalized gradient approximation made simple}, sub.
C     to Phys. Rev.Lett. May 1996.
C [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
C     construction of a generalized gradient approximation:  The PW91
C     density functional}, submitted to Phys. Rev. B, Feb. 1996.
C [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      logical lgga
C thrd*=various multiples of 1/3
C numbers for use in LSD energy spin-interpolation formula, [c](9).
C      GAM= 2^(4/3)-2
C      FZZ=f''(0)= 8/(9*GAM)
C numbers for construction of PBE
C      gamma=(1-log(2))/pi^2
C      bet=coefficient in gradient expansion for correlation, [a](4).
C      eta=small number to stop d phi/ dzeta from blowing up at 
C          |zeta|=1.      
      parameter(thrd=1.0d0/3.0d0,thrdm=-thrd,thrd2=2.0d0*thrd)
      parameter(sixthm=thrdm/2.0d0)
      parameter(thrd4=4.0d0*thrd)
      parameter(GAM=0.51984209978974632953442121455650d0)
      parameter(fzz=8.0d0/(9.0d0*GAM))
      parameter(gamma=0.031090690869654895034940863712730d0)
      parameter(bet_mb=0.066724550603149220d0)
      parameter(eta=1.d-120)
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C find LSD energy contributions, using [c](10) and Table I[c].
C EU=unpolarized LSD correlation energy
C EURS=dEU/drs
C EP=fully polarized LSD correlation energy
C EPRS=dEP/drs
C ALFM=-spin stiffness, [c](3).
C ALFRSM=-dalpha/drs
C F=spin-scaling factor from [c](9).
C construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor22(0.03109070d0,0.213700d0,7.59570d0,3.58760d0, 
     &    1.63820d0,0.492940d0,rtrs,EU,EURS)
      CALL gcor22(0.015545350d0,0.205480d0,14.11890d0,6.19770d0, 
     &    3.36620d0,0.625170d0,rtRS,EP,EPRS)
      CALL gcor22(0.01688690d0,0.111250d0,10.3570d0,3.62310d0,
     &    0.880260d0,0.496710d0,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1.0d0+ZET)**THRD4+(1.0d0-ZET)**THRD4-2.0d0)/GAM
      EC = EU*(1.0d0-F*Z4)+EP*F*Z4-ALFM*F*(1.0d0-Z4)/FZZ
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C LSD potential from [c](A1)
C ECRS = dEc/drs [c](A2)
C ECZET=dEc/dzeta [c](A3)
C FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1.0d0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.0d0-Z4)/FZZ
      FZ = THRD4*((1.0d0+ZET)**THRD-(1.0d0-ZET)**THRD)/GAM
      ECZET = 4.0d0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU 
     &        -(1.0d0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.0d0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
Cwrite(*,*)'rs,VCUP,VCDN',rs,VCUP,VCDN
      if (.not.lgga) return
C----------------------------------------------------------------------
C PBE correlation energy
C G=phi(zeta), given after [a](3)
C DELT=bet/gamma
C B=A of [a](8)
      bet = bet_mb*(1.0d0 + 0.10d0*RS)/(1.0d0 + 0.17780d0*RS)

      delt=bet/gamma
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1.0d0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.0d0+B*T2
      Q5 = 1.0d0+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1.0d0+DELT*Q4*T2/Q5)
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3.0d0
      GZ=(((1.0d0+zet)**2+eta)**sixthm- 
     &((1.0d0-zet)**2+eta)**sixthm)/3.0d0
      FAC = DELT/B+1.0d0
      BG = -3.0d0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.0d0+2.0d0*B*T2
      hB = -BET*G3*B*T6*(2.0d0+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hZ = 3.0d0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2.0d0*BET*G3*Q9/Q8
      COMM = H+HRS-7.00d0*T2*HT/6.0d0
      PREF = HZ-GZ*T2*HT/G
      COMM = COMM-PREF*ZET
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      ecdd=0.50d0/(sk*g)*t*ht
      RETURN
      END 

C***********************************************************************
C
C Written by Jianwei Sun and Yoon-Suk Kim 06/18/2009
C
C RU,RD                        density up,down
C DRU, DRD                     abs. val. gradient of density up/down
C DRT                          abs. val. gradient of total density
C EPPGGA                       GGA ENERGY PER PARTICLE
C EPPGGA_D1,EPPGGA_D2          THE DERIVATIVE OF EPPGGA WRT n
C EPPGGA_DD1,EPPGGA_DD2        THE DERIVATIVE OF EPPGGA WRT |GRAD n|
C
C***********************************************************************
      SUBROUTINE CORGGA_1(
     &   RU,RD,DRU,DRD,DRT, 
     &   EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD1,EPPGGA_DD2)


      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ONE=1.0d0)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (FIVE=5.0d0)
      PARAMETER (SIX=6.0d0)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)
      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)
  
      PARAMETER (GAMMA=0.031090690869654895034940863712730d0)
      PARAMETER (BETA_mb=0.066724550603149220d0)

      EPPGGA  =0.0d0
      EPPGGA_D1=0.0d0
      EPPGGA_D2=0.0d0
      EPPGGA_DD1=0.0d0
      EPPGGA_DD2=0.0d0

      RT=RU+RD
      YA=DRU**2.0d0
      YB=DRD**2.0d0
      Y=DRT**2.0d0
C     YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2.0d0

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.99999999999990d0),0.99999999999990d0)
      DZETAD1=TWO*RD/RT**TWO
      DZETAD2=-TWO*RU/RT**TWO


      DTHRD=exp(log(RT)*THRD)
      RS=(0.750d0/PI)**THRD/DTHRD
      DRSD1=-THRD/RT*RS
      DRSD2=-THRD/RT*RS


      PHI = (exp((TWO*THRD)*log(1.0d0+ZETA)) 
     $       +exp((TWO*THRD)*log(1.0d0-ZETA)))/2.0d0
      D_PHI_DZETA=THRD*((ONE+ZETA)**(-THRD)-(ONE-ZETA)**(-THRD))
      D_PHI_D1=D_PHI_DZETA*DZETAD1
      D_PHI_D2=D_PHI_DZETA*DZETAD2

      AFIX_T=SQRT(PI/FOUR)*(9.0d0*PI/FOUR)**(ONE/SIX)
      S=DRT/(TWO*(THREE*PISQ)**THRD*RT**THRD4)
      DSD1=-THRD4*S/RT
      DSD2=-THRD4*S/RT
C |GRAD N|/|GRAD N(UP OR DOWN)|
      GNGNU=ONE/(TWO*DRT)*(TWO*DRU+TWO*YC/DRU)
      GNGND=ONE/(TWO*DRT)*(TWO*DRD+TWO*YC/DRD)
!YY. work around for DRU(DRD) == 0 
!      if(DRU .eq. 0.0d0) GNGNU=DRU/(TWO*DRT)
!      if(DRD .eq. 0.0d0) GNGND=DRD/(TWO*DRT)
!end work around
      DSDD1=S/DRT*GNGNU
      DSDD2=S/DRT*GNGND

      T=AFIX_T*S/SQRT(RS)/PHI
      T2 = T*T
      T4 = T2*T2

      DTD1=AFIX_T*(PHI*RS*DSD1-0.50d0*S*PHI*DRSD1-S*RS*D_PHI_D1)
     $           /(PHI**TWO*RS**(THREE/TWO))
      DTD2=AFIX_T*(PHI*RS*DSD2-0.50d0*S*PHI*DRSD2-S*RS*D_PHI_D2)
     $           /(PHI**TWO*RS**(THREE/TWO))
      DTDD1=AFIX_T*DSDD1/(SQRT(RS)*PHI)
      DTDD2=AFIX_T*DSDD2/(SQRT(RS)*PHI)

      FK=(3.0d0*PI*PI)**THRD*DTHRD
      SK = SQRT(4.00d0*FK/PI)
C Only the local part (EC,VCUPLDA,VCDNLDA) are used from the following call to CORPBE
      CALL CORPBE_revtpss(RS,ZETA,EC,VCUPLDA,VCDNLDA,PHI,SK, 
     $           T,H,DVCUP,DVCDN,ECQ,.FALSE.)
C OBTAIN THE CONTRIBUTION FROM LDA PART AND D_EC_DRS AND D_EC_DZETA
C      EPPGGA=EC
C      EPPGGA_D1=(VCUPLDA-EC)/RT
C      EPPGGA_D2=(VCDNLDA-EC)/RT
      D_EC_DZETA=(VCUPLDA-VCDNLDA)/TWO
      D_EC_DRS=(EC-ZETA*D_EC_DZETA-(VCUPLDA+VCDNLDA)/TWO)*THREE/RS

C THE RS DEPENDENCE OF BETA
      AFACTOR=0.10d0
      BFACTOR=0.17780d0
      BETA_NUM=ONE + AFACTOR*RS
      BETA_DEN=ONE+ BFACTOR*RS
      D_BETA_NUM=AFACTOR
      D_BETA_DEN=BFACTOR
      BETA = BETA_MB*BETA_NUM/BETA_DEN
      D_BETA_DRS=BETA_MB*(BETA_DEN*D_BETA_NUM-BETA_NUM*D_BETA_DEN)
     $                  /BETA_DEN**TWO
      
 
      PHI3=PHI**THREE
      PON=-EC/(PHI3*gamma)
      W=DEXP(PON)-ONE
      D_W_DRS=-(W+ONE)*D_EC_DRS/(GAMMA*PHI3)
      D_W_DZETA=-(W+ONE)/(GAMMA*PHI3)
     $                  *(D_EC_DZETA-THREE*EC*D_PHI_DZETA/PHI)
      D_W_DT=0.0d0


      A=BETA/(GAMMA*W)
      D_A_DRS=(W*D_BETA_DRS-BETA*D_W_DRS)/(GAMMA*W**TWO)
      D_A_DZETA=-BETA*D_W_DZETA/(GAMMA*W**TWO)
      D_A_DT=0.0d0

      V=A*T2
      D_V_DRS=T2*D_A_DRS
      D_V_DZETA=T2*D_A_DZETA
      D_V_DT=T2*D_A_DT+TWO*A*T
        
C      FUNKG=ONE/(ONE+V+V**TWO)
C      D_FUNKG_DV=-(ONE+TWO*V)*FUNKG**TWO

      FUNKG=ONE/(ONE+4.0d0*V)**0.250d0
      D_FUNKG_DV=-ONE/(ONE+4.0d0*V)**1.250d0


      
      HCORE=ONE+W*(ONE-FUNKG)
      AH=GAMMA*PHI3
      H=AH*DLOG(HCORE)
   
      DH1=ONE-FUNKG
      DH2=W*D_FUNKG_DV
      D_H_DRS=(DH1*D_W_DRS-DH2*D_V_DRS)*AH/HCORE
      D_H_DZETA=THREE*H*D_PHI_DZETA/PHI
     $          +(DH1*D_W_DZETA-DH2*D_V_DZETA)*AH/HCORE
      D_H_DT=(DH1*D_W_DT-DH2*D_V_DT)*AH/HCORE

C OUTPUT EPPGGA AND ITS DERIVATIVES EPPGGA_D1,EPPGGA_D2, EPPGGA_DD1, EPPGGA_DD2
      EPPGGA=EC+H
      EPPGGA_D1=D_EC_DRS*DRSD1+D_EC_DZETA*DZETAD1+D_H_DRS*DRSD1
     $                        +D_H_DZETA*DZETAD1+D_H_DT*DTD1
      EPPGGA_D2=D_EC_DRS*DRSD2+D_EC_DZETA*DZETAD2+D_H_DRS*DRSD2
     $                        +D_H_DZETA*DZETAD2+D_H_DT*DTD2
      EPPGGA_DD1=D_H_DT*DTDD1
      EPPGGA_DD2=D_H_DT*DTDD2

      RETURN
      END 



C***********************************************************************
C
C Written by Jianwei Sun 04/01/2013
C
C RU,RD                        density up,down
C DRU, DRD                     abs. val. gradient of density up/down
C DRT                          abs. val. gradient of total density
C EPPGGA                       GGA ENERGY PER PARTICLE
C EPPGGA_D1,EPPGGA_D2          THE DERIVATIVE OF EPPGGA WRT n
C EPPGGA_DD1,EPPGGA_DD2        THE DERIVATIVE OF EPPGGA WRT |GRAD n|
C
C***********************************************************************
      SUBROUTINE CORGGA_0(RU,RD,DRU,DRD,DRT, 
     &   EPPGGA,EPPGGA_D1,EPPGGA_D2,EPPGGA_DD1,EPPGGA_DD2)

      IMPLICIT REAL*8 (A-H,O-Z)

C the coefficients for the correlatin part
      PARAMETER (CFB1=0.0285764d0)
      PARAMETER (CFB2=0.0889d0) 
      PARAMETER (CFB3=0.125541d0)
      PARAMETER (cfkaiLD=0.12802585262625815d0)

      PARAMETER (ONE=1.0d0)
      PARAMETER (TWO=2.0d0)
      PARAMETER (THREE=3.0d0)
      PARAMETER (FOUR=4.0d0)
      PARAMETER (FIVE=5.0d0)
      PARAMETER (SIX=6.0d0)
      PARAMETER (THRD=1.0d0/3.0d0)
      PARAMETER (THRD2=2.0d0*THRD)
      PARAMETER (THRD4=4.0d0*THRD)
      PARAMETER (THRD5=1.0d0+THRD2)
      PARAMETER (THRD8=1.0d0+THRD5)
      PARAMETER (PI =3.141592653589793238d0)
      PARAMETER (PISQ=PI*PI)

      PARAMETER (GAMMA=0.031090690869654895034940863712730d0)
      PARAMETER (BETA_mb=0.066724550603149220d0)

      EPPGGA  =0.0d0
      EPPGGA_D1=0.0d0
      EPPGGA_D2=0.0d0
      EPPGGA_DD1=0.0d0
      EPPGGA_DD2=0.0d0
      AX_LDA=-THREE/(FOUR*PI)*(9.00d0*PI/FOUR)**THRD

      RT=RU+RD
      YA=DRU**2.0d0
      YB=DRD**2.0d0
      Y=DRT**2.0d0
C     YC IS DEL RU DOT DEL RD
      YC=(Y-YA-YB)/2.0d0

      ZETA=(RU-RD)/RT
      ZETA=MIN(MAX(ZETA,-0.99999999999990d0),0.99999999999990d0)
      DZETAD1=TWO*RD/RT**TWO
      DZETAD2=-TWO*RU/RT**TWO


      DTHRD=exp(log(RT)*THRD)
      RS=(0.750d0/PI)**THRD/DTHRD
      DRSD1=-THRD/RT*RS
      DRSD2=-THRD/RT*RS


      PHI = (exp((TWO*THRD)*log(1.0d0+ZETA)) 
     $           +exp((TWO*THRD)*log(1.0d0-ZETA)))/2.0d0
      D_PHI_DZETA=THRD*((ONE+ZETA)**(-THRD)-(ONE-ZETA)**(-THRD))
      D_PHI_D1=D_PHI_DZETA*DZETAD1
      D_PHI_D2=D_PHI_DZETA*DZETAD2

      AFIX_T=SQRT(PI/FOUR)*(9.0d0*PI/FOUR)**(ONE/SIX)
      S=DRT/(TWO*(THREE*PISQ)**THRD*RT**THRD4)
      DSD1=-THRD4*S/RT
      DSD2=-THRD4*S/RT
C |GRAD N|/|GRAD N(UP OR DOWN)|
      if(abs(DRT) < 1.0d-100) DRT=1.0d-100
      !if(DRU.eq.0.0d0) DRU=1.0d-50
      !if(DRD.eq.0.0d0) DRD=1.0d-50
      GNGNU=ONE/(TWO*DRT)*(TWO*DRU+TWO*YC/DRU)
      GNGND=ONE/(TWO*DRT)*(TWO*DRD+TWO*YC/DRD)
!YY. work around for DRU(DRD) == 0
!# eqn:    (2*DRU+(DRT**2-DRU**2-DRD**2)/DRU)/(2*DRT)
!#         =(DRT**2-DRD**2)/(2*DRU*DRT)
!      if(DRU .eq. 0.0d0) GNGNU=DRU/(TWO*DRT)
!      if(DRD .eq. 0.0d0) GNGND=DRD/(TWO*DRT)
!end work around

      DSDD1=S/DRT*GNGNU
      DSDD2=S/DRT*GNGND
      !write(*,*) 'SCAN',YC,DRD 0/0 problem if DRU or DRD is zero  
      ! In that case, you have
      ! YC/DRD = (DRT**2-DRU**2-DRD**2)/(2*DRD) = (0 - DRD**2)/(2*DRD) = - DRD/2
C--------------------------------------------------------------------------------------
C EC0 and its derivatives
C--------------------------------------------------------------------------------------

C EC0LDA and its derivatives
      RSHALF=SQRT(RS)
      FACTOR1=ONE+CFB2*RSHALF+CFB3*RS
      EC0LDA=-CFB1/FACTOR1

      D_EC0LDA_DRS=(CFB3+CFB2/RSHALF/TWO)*EC0LDA**TWO/cfb1
      
C gc(zeta) and its derivatives
      dx_zeta=((one+zeta)**THRD4+(one-zeta)**THRD4)/two
      d_dx_zeta_dzeta=two*((one+zeta)**THRD-(one-zeta)**THRD)/three
      TWO13=TWO**THRD
c      gc_zeta=(two13-dx_zeta)/(two13-one)
c      d_gc_zeta_dzeta=-d_dx_zeta_dzeta/(two13-one)

      gc_zeta=(one-2.363d0*(dx_zeta-one))*(one-zeta**12.0d0)
      d_gc_zeta_dzeta=-(one-2.363d0*(dx_zeta-one))*12.d0*zeta**11.d0
      d_gc_zeta_dzeta=d_gc_zeta_dzeta
     $                -2.363d0*d_dx_zeta_dzeta*(one-zeta**12.0d0)

C h0 and its derivatives
      w0=dexp(-EC0LDA/cfb1)-one
      d_w0_drs=-(w0+one)*D_EC0LDA_DRS/cfb1
      
      gfunkinf=one/(one+four*cfkaiLD*s*s)**(one/four)
      d_gfunkinf_ds=-two*cfkaiLD*s*gfunkinf**five

      hcore0=one+w0*(one-gfunkinf)
      h0=cfb1*dlog(hcore0)
      d_h0_drs=cfb1*d_w0_drs*(one-gfunkinf)/hcore0
      d_h0_ds=-cfb1*w0*d_gfunkinf_ds/hcore0


C OUTPUT EPPGGA AND ITS DERIVATIVES EPPGGA_D1,EPPGGA_D2, EPPGGA_DD1, EPPGGA_DD2
      EPPGGA=(EC0LDA+h0)*gc_zeta

      d_EPPGGA_dzeta=(EC0LDA+h0)*d_gc_zeta_dzeta
      d_EPPGGA_drs=(D_EC0LDA_DRS+d_h0_drs)*gc_zeta
      d_EPPGGA_ds=d_h0_ds*gc_zeta

      EPPGGA_D1=d_EPPGGA_dzeta*DZETAD1+d_EPPGGA_drs*DRSD1
     $          +d_EPPGGA_ds*DSD1
      EPPGGA_D2=d_EPPGGA_dzeta*DZETAD2+d_EPPGGA_drs*DRSD2
     $          +d_EPPGGA_ds*DSD2



      EPPGGA_DD1=d_EPPGGA_ds*DSDD1
      EPPGGA_DD2=d_EPPGGA_ds*DSDD2
      
      RETURN
      END 




C----------------------------------------------------------------------
C slimmed down version of GCOR used in PW91 routines, to interpolate
C LSD correlation energy, as given by (10) of
C J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
C K. Burke, May 11, 1996.
C----------------------------------------------------------------------
      SUBROUTINE GCOR22(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
      IMPLICIT REAL*8 (A-H,O-Z)
      Q0 = -2.0d0*A*(1.0d0+A1*rtrs*rtrs)
      Q1 = 2.0d0*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1.0d0+1.0d0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2.0d0*B2+rtrs*(3.0d0*B3+4.0d0*B4*rtrs))
      GGRS = -2.0d0*A*A1*Q2-Q0*Q3/(Q1*(1.0d0+Q1))
      RETURN
      END


