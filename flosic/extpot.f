C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE EXTPOT(R,POT,DERIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXT=10)
      COMMON/EXPOT/C(MAXT),M(MAXT),N(3,MAXT),A(3,MAXT),G(3,MAXT),NTERMS
      LOGICAL FIRST
      DIMENSION R(3),DERIV(3)
      DATA FIRST/.TRUE./
      DATA EPS/ 1.0D-7/
      IF(FIRST)THEN
        NTERMS=0
        CALL READEXT
        FIRST=.FALSE.
      END IF
      POT=0.0D0
      DO J=1,3
      DERIV(J)=0.0D0
      END DO
      DO I=1,NTERMS
          XX=R(1)-A(1,I)
          YY=R(2)-A(2,I)
          ZZ=R(3)-A(3,I)
          RR=SQRT(XX*XX+YY*YY+ZZ*ZZ)
                RP=1.0D0
                DO J=1,M(I)
                RP=RP/RR
                END DO
          POLX=1.0D0
          POLY=1.0D0
          POLZ=1.0D0
          DERX=1.0D0
          DERY=1.0D0
          DERZ=1.0D0
          DO J=1,N(1,I)
          POLX=POLX*XX
             IF(J.GT.1)DERX=DERX*XX
          END DO
          DO J=1,N(2,I)
          POLY=POLY*YY
             IF(J.GT.1)DERY=DERY*YY
          END DO
          DO J=1,N(3,I)
          POLZ=POLZ*ZZ
             IF(J.GT.1)DERZ=DERZ*ZZ
          END DO
          DERX=N(1,I)*DERX
          DERY=N(2,I)*DERY
          DERZ=N(3,I)*DERZ
          X2=XX*XX
          Y2=YY*YY
          Z2=ZZ*ZZ
          ENV=EXP(-G(1,I)*X2-G(2,I)*Y2-G(3,I)*Z2)
          POT=POT+C(I)*POLX*POLY*POLZ*ENV*RP
          DERIV(1)=DERIV(1)+C(I)*DERX*POLY*POLZ*ENV*RP
          DERIV(2)=DERIV(2)+C(I)*POLX*DERY*POLZ*ENV*RP
          DERIV(3)=DERIV(3)+C(I)*POLX*POLY*DERZ*ENV*RP
          DERIV(1)=DERIV(1)-C(I)*POLX*POLY*POLZ*ENV*2.0*G(1,I)*XX*RP
          DERIV(2)=DERIV(2)-C(I)*POLX*POLY*POLZ*ENV*2.0*G(2,I)*YY*RP
          DERIV(3)=DERIV(3)-C(I)*POLX*POLY*POLZ*ENV*2.0*G(3,I)*ZZ*RP
          IF(RR.GT.EPS)THEN
          DERIV(1)=DERIV(1)-C(I)*POLX*POLY*POLZ*ENV*M(I)*XX*RP/RR/RR
          DERIV(2)=DERIV(2)-C(I)*POLX*POLY*POLZ*ENV*M(I)*YY*RP/RR/RR
          DERIV(3)=DERIV(3)-C(I)*POLX*POLY*POLZ*ENV*M(I)*ZZ*RP/RR/RR
          END IF
      END DO
      END 
