C UTEP Electronic Structure Lab (2020)
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE PBEEX(MODE,RHO,S,U,V,LGGA,LPOT,EXL,EXN,VX)
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT MODE: RUN PBE(MODE=1),REVPBE(MODE=2) OR RPBE(MODE=3)
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (LOCAL: EXL, NONLOCAL: EXN) 
C           AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(THRD=1.D0/3.D0,THRD4=4.D0/3.D0)
Cdvp  parameter(pi=3.14159265358979323846264338327950d0)
      PARAMETER(AX=-0.738558766382022405884230032680836D0)
C      PARAMETER(UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK)
      PARAMETER(UM=0.2195149727645171D0)
      SAVE

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      IF(MODE.EQ.2) THEN
C MODE=2 IS FOR REVPBE
        UK=1.2450D0
      ELSE
C MODE=1 PBE AND MODE=3 RPBE USE THIS
        UK=0.8040D0
      ENDIF
      UL=UM/UK
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      EXUNIF = AX*RHO**THRD
      IF(LGGA.EQ.0)THEN
        EX=EXUNIF
        VX=EX*THRD4
        RETURN
      ENDIF
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct PBE enhancement factor
      S2 = S*S
      IF(MODE.EQ.3)THEN
C RPBE USES THIS
        P0=EXP(-UL*S2)
        FXPBE=1D0+UK*(1.0D0-P0)
      ELSE
C PBE AND REVPBE USE THIS
        P0=1.D0+UL*S2
        FXPBE = 1D0+UK-UK/P0
      ENDIF
      EXL = EXUNIF
      EXN = EXUNIF*(FXPBE-1.0D0)
      IF(LPOT.EQ.0)RETURN
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c  Fss=d Fs/ds
      IF(MODE.EQ.3)THEN
C RPBE USES THIS
        FS=2.D0*UM*P0
        FSS=-2.D0*UL*S*FS
      ELSE
C PBE AND REVPBE USE THIS
        FS=2.D0*UK*UL/(P0*P0)
        FSS=-4.D0*UL*S*FS/P0
      ENDIF
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c calculate potential from [b](24) 
      VX = EXUNIF*(THRD4*FXPBE-(U-THRD4*S2*S)*FSS-V*FS)
      RETURN
      END
