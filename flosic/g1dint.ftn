C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       SUBROUTINE G1DINT(X1,X2,NPOW,ALP,GIRES)
C
C G1DINT CALCULATES THE INTEGRAL
C
C   X2  N            2
C  INT X * EXP(-ALP*X ) 
C   X1
C
C FOR ALL N=0,NPOW. X1 AND X2 MUST BE >= 0, X2 MUST BE >= X1.
C WRITTEN BY DV POREZAG
C
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (MAXPOW=20)
        PARAMETER (MAXORD=30)
        DIMENSION GIRES(0:NPOW)
        DIMENSION U(2),U2(2),ERF1(2),ERF2(2),UEXP(2)
        DIMENSION RCINT(MAXPOW+2*MAXORD+2)
        SAVE
        DATA NIRC/0/
        DATA RTPIRC/0.56418958354775628695D0/
        DATA RTPIHF/0.88622692545275801364D0/
        DATA TOLER/1.0D-14/       
C
C CHECKS 
C     
        IF (NPOW .LT. 0) RETURN
        IF (ALP .LE. 0.0D0) THEN
         write(6,*)'G1DINT: EXPONENTS MUST BE > 0'
         CALL STOPIT
        END IF
        IF (NPOW .GT. MAXPOW) THEN
         write(6,*)'G1DINT: MAXPOW MUST BE AT LEAST: NPOW'
         CALL STOPIT
        END IF
        IF (X1 .LT. 0.0D0) THEN
         write(6,*)'G1DINT: LOWER INTEGRAL BOUNDARY MUST BE >= 0'
         CALL STOPIT
        END IF
        IF (X2 .LE. X1) THEN
         DO IPOW=0,NPOW
          GIRES(IPOW)= 0.0D0
         END DO
         RETURN
        END IF
C
C SETUP
C
        NEED=NPOW+2*MAXORD+2
        IF (NIRC .LT. NEED) THEN
         DO I=NIRC+1,NEED
          RCINT(I)= 1.0D0/I
         END DO
         NIRC=NEED
        END IF
        ALPS= SQRT(ALP)
        U(1)= ALPS*X1
        U(2)= ALPS*X2
        U2(1)= U(1)*U(1)
        U2(2)= U(2)*U(2)
C
C IF U2 IS SMALL, USE EXPANSION OF CORRECT SOLUTION ABOUT U2=0
C
        IF (U2(2) .LE. 1.0D0) THEN
         DO IPOW= 0,NPOW
          GIRES(IPOW)= 0.0D0
         END DO
         U2P1= 1.0D0
         U2P2= 1.0D0
         FAC= 1.0D0
         DO IORD=0,MAXORD
          IOR2P1= 2*IORD+1
          XP1= X1
          XP2= X2
          FCT= FAC*RCINT(IOR2P1)
          OLD= GIRES(NPOW)
          DO IPOW= 0,NPOW
           GIRES(IPOW)=GIRES(IPOW)+FCT*(XP2*U2P2-XP1*U2P1)
           XP1= XP1*X1
           XP2= XP2*X2
           FCT= FCT*(IPOW+IOR2P1)*RCINT(IPOW+IOR2P1+1)
          END DO
          IF (ABS(GIRES(NPOW)-OLD) .LT. TOLER*ABS(OLD)) GOTO 10
          FAC= -FAC*RCINT(IORD+1)
          U2P1= U2P1*U2(1)
          U2P2= U2P2*U2(2)
         END DO
   10    CONTINUE
        ELSE
C
C DETERMINE INTEGRALS THE REGULAR WAY: ERROR FUNCTION FOR IPOW=0,
C ANALYTICAL SOLUTION FOR IPOW=1, RECURSION FOR IPOW > 1
C
         DO 60 ILIM=1,2
          XS= U(ILIM)
          XS2= U2(ILIM)
          UEXP(ILIM)= EXP(-XS2)
          IF (XS2 .LE. 4.0D0) THEN
           FAC= 2*RTPIRC*UEXP(ILIM)
           SUM= 0.0D0
           I=0
   20      CONTINUE
            I=I+1
            OLD= SUM
            SUM= SUM+FAC
            IF (FAC .LE. TOLER*OLD) GOTO 30
            FAC= FAC*XS2/(0.5D0+I)
            GOTO 20
   30      CONTINUE
           ERF1(ILIM)=SUM*XS
           ERF2(ILIM)=1.0D0-ERF1(ILIM)
          ELSE
           GOLD= 0.0D0
           A0= 1.0D0
           B0= 0.0D0
           A1= XS2
           B1= 1.0D0
           FAC= 1.0D0
           I=0
   40      CONTINUE
            I=I+1
            AI= I
            AIA= AI-0.5D0
            A0= (A1+A0*AIA)*FAC
            B0= (B1+B0*AIA)*FAC
            AIF= AI*FAC
            A1= XS2*A0+AIF*A1
            B1= XS2*B0+AIF*B1
            IF (A1 .NE. 0.0D0) THEN
             FAC= 1.0D0/A1
             G= B1*FAC
             IF (ABS(G-GOLD) .LE. TOLER*GOLD) GOTO 50
             GOLD= G
            END IF
            GOTO 40
   50      CONTINUE
           ERF2(ILIM)= RTPIRC*XS*UEXP(ILIM)*G
           ERF1(ILIM)= 1.0D0-ERF2(ILIM)
          END IF
   60    CONTINUE
         IF (ERF1(1) .LT. 0.5D0) THEN
          GIRES(0)= RTPIHF*(ERF1(2)-ERF1(1))
         ELSE
          GIRES(0)= RTPIHF*(ERF2(1)-ERF2(2))
         END IF
         IF (NPOW .GT. 0) GIRES(1)= 0.5D0*(UEXP(1)-UEXP(2))
C
C RECURSION
C
         U1PW= UEXP(1)
         U2PW= UEXP(2)
         DO IPOW=2,NPOW
          U1PW= U1PW*U(1) 
          U2PW= U2PW*U(2) 
          GIRES(IPOW)= 0.5D0*(U1PW-U2PW+(IPOW-1)*GIRES(IPOW-2))
         END DO
         FAC= 1.0D0/ALPS
         FCT= FAC
         DO IPOW=0,NPOW
          GIRES(IPOW)= GIRES(IPOW)*FCT
          FCT= FCT*FAC
         END DO
        END IF
        RETURN
       END
