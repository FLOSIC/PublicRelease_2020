C UTEP Electronic Structure Lab (2020)
C
C *********************************************************************
C
C GTGRAD CALCULATES THE NEGATIVE GRADIENT OF THE ATOMIC BASIS 
C FUNCTIONS OF ONE ATOM FOR ONE POINT
C
       SUBROUTINE GTGRAD(X,Y,Z,I_ATOM,GRAD)
C ORIGINAL VERSION BY KOBLAR A. JACKSON (1988)
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:52 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I_ATOM, I, I_BARE, IANG, IBEG, IEND, IX, J, LI,
     & LMAX1, MU, NB
       REAL*8 :: X , Y, Z, GRAD, ALF, ALF2, ANG, DANG, GAUS, GEPS,
     & PANG, REXPON, RSQR
       SAVE
       DIMENSION GRAD(3,MAX_CON,10,3),ANG(10),PANG(3,10)
     &  ,DANG(3,10),ALF2(3)
       DIMENSION IBEG(3), IEND(3)
C
       DATA IBEG /1,2,5/
       DATA IEND /1,4,10/
C
       LMAX1 = LSYMMAX(I_ATOM)+1
       GEPS = EXP(-CUTEXP)
C
C GET ANGULAR FUNCTION ARRAY
C
       ANG(1) = 1.0D0
       ANG(2) = X  
       ANG(3) = Y  
       ANG(4) = Z
       ANG(5) = X*X
       ANG(6) = Y*Y
       ANG(7) = Z*Z
       ANG(8) = X*Y
       ANG(9) = X*Z
       ANG(10)= Y*Z
       RSQR=ANG(5)+ANG(6)+ANG(7)
C
C GET DERIVATIVE ARRAY FOR POLYNOMIAL FUNCTIONS:  -D/DX(ANG)
C
       DO 3 I=1,3
        DO 2 J=1,10
         PANG(I,J)=0.0D0
    2   CONTINUE
    3  CONTINUE
C
C X-DERIVATIVE
C
       PANG(1,2) = -1.0D0
       PANG(1,5) = -2*X
       PANG(1,8) = -Y
       PANG(1,9) = -Z
C
C  Y DERIVATIVE
C
       PANG(2,3) = -1.0D0
       PANG(2,6) = -2*Y
       PANG(2,8) = -X
       PANG(2,10)= -Z
C
C  Z DERIVATIVE
C
       PANG(3,4) = -1.0D0
       PANG(3,7) = -2*Z
       PANG(3,9) = -X
       PANG(3,10)= -Y
C
C INITIALIZE GRAD
C
       DO 20 LI=1,LMAX1
        DO 15 IANG=1,IEND(LI)
         DO 10 NB=1,N_CON(LI,I_ATOM)
          DO 5 IX=1,3
           GRAD(IX,NB,IANG,LI)=0.0D0
    5     CONTINUE
   10    CONTINUE
   15   CONTINUE
   20  CONTINUE
C
C  GET GRADIENT TERM FOR DIFFERENT ANGULAR FUNCTIONS
C
       DO 70 I_BARE=1,N_BARE(I_ATOM)
        ALF=BFALP(I_BARE,I_ATOM)
        REXPON=ALF*RSQR
        IF (REXPON.GT.CUTEXP) GOTO 70
        ALF2(1)=2*ALF*X
        ALF2(2)=2*ALF*Y
        ALF2(3)=2*ALF*Z
        GAUS=EXP(-REXPON)
C
C X,Y,Z DERIVATIVE
C
        DO 30 I=1,3
         DO 25 J=1,IEND(LMAX1)
          DANG(I,J) = ALF2(I)*ANG(J)+PANG(I,J)
   25    CONTINUE
   30   CONTINUE
C
C CALCULATING DERIVATIVE
C
        DO 60 LI=1,LMAX1
         DO 55 NB=1,N_CON(LI,I_ATOM)
          IF (ABS(BFCON(I_BARE,NB,LI,I_ATOM)*GAUS).GT.GEPS) THEN
           DO 50 MU=IBEG(LI),IEND(LI)
            DO 45 IX=1,3
             GRAD(IX,NB,MU,LI)=GRAD(IX,NB,MU,LI) 
     &       +BFCON(I_BARE,NB,LI,I_ATOM)*GAUS*DANG(IX,MU)
   45       CONTINUE
   50      CONTINUE
          END IF
   55    CONTINUE
   60   CONTINUE
   70  CONTINUE
       RETURN
       END
