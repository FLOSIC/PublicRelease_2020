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

! modSCANx routine - Currently it does not returns potential.

! YY. Edit the VSCANx functional into spin unpolarized functional
! to make it work better on NRLMOL

       SUBROUTINE MODVSCANxUNP(RHO,DRHO,TAU_RHO,
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

      !open(888,file='modscanfitparam',form='formatted',status='old')
      !rewind(888)
      !  read(888,*) FPARAMA,FPARAMB,FPARAMC
      !close(888)

      VXD1=0.0d0
      !VXD2=0.0d0
      VXDD1=0.0d0
      !VXDD2=0.0d0
      AMUXD1=0.0d0
      !AMUXD2=0.0d0

!YY. Added these to prevent NAN
      !if(RHO == 0.0d0) RHO = 1d-30
      !if(DRHO == 0.0d0) DRHO = 1d-30
      !DRHO = DRHO /1.259921049894873d0

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

!In modSCAN, alpha=0, s_{i\sigma} is divided by 2^(1/3) in Fx^{SCAN}
      P=(DRHO)**TWO/(FOUR*(THREE*PISQ)**THRD2*(RHO)**THRD8)
!     P=(DRHO/1.259921049894873d0)**TWO/(FOUR*(THREE*PISQ)**THRD2*
!    &  (RHO)**THRD8)
      TAU_UNIF=THREE/10.0d0*(THREE*PI**TWO)**THRD2*RHO**THRD5
!      ALPHA=(TAU_RHO-TAUW_RHO)/TAU_UNIF
      ALPHA=0.0d0
!Beta
      BETA =(TAU_RHO-TAUW_RHO)/(TAU_RHO+TAU_UNIF)

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
     $                p, alpha, Fx, fx1p, fx1a, fx2p2, fx2pa, fx2a2,
     &                d_Fx1_da,FA,d_hx1_da,d_FA_da,hx0,hx1)


C Calculate fx DERIVATIVES WRT n, |grad n|, and tau
      DFX_D=fx1p*DPD+fx1a*DALPHAD  !<- YY. Testing this now - fx1a seems to be responsible for kink
      DFX_DD=fx1p*DPDD+fx1a*DALPHADD
      DFX_DTAU=fx1a*DALPHADTAU

C   OUTPUT THE VXD1,VXDD1 AND AMUXD1
      VXD1=EXDLDA*FX+EXLDA*DFX_D   !<- YY. DFX_D causes kink
      VXDD1=EXLDA*DFX_DD
      AMUXD1=EXLDA*DFX_DTAU
! YY. orginal lines return Exchange energy for given mesh point
! we want exchange energy per particle instead in NRLMOL.
      !EX_metagga=0.0d0
      !EX_metagga=EX_metagga+EXLDA*FX
      EX_metagga=EXUNIF*FX !*1.259921049894873d0
!     &    *(1.0d0+FPARAMA*BETA+FPARAMB*BETA*BETA+FPARAMC*BETA*BETA*BETA)
      !write(*,'("SCANUNP",5F16.8)') RHO,sqrt(p),Fx,EXUNIF,EX_metagga

      RETURN
      END
      
      

