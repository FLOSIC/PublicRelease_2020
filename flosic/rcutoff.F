C UTEP Electronic Structure Lab (2020)
C
C ******************************************************************
C
C RCUTOFF CALCULATES THE ALLOWED CUTOFF RADIUS FOR A TEST FUNCTION 
C F(R)= R**N * EXP(-AMIN*R**2) SO THAT 
C INTEGRAL(RCUTOFF->INFINITY)/INTEGRAL(0->INFINITY) = EPS
C
C N,AMIN,EPS: SEE ABOVE FORMULA
C
       FUNCTION RCUTOFF(N,AMIN,EPS)
        IMPLICIT REAL*8 (A-H,O-Z)
        SAVE
        DATA ACCU /5.0D-16/
        DATA ZERO /1.0D-30/
C
        IF (AMIN .LT. ZERO) THEN
         write(6,*)'RCUTOFF: AMIN IS TOO SMALL'
         CALL STOPIT
        END IF
C
C MACHINE ACCURACY
C
        CALL MACHEPS(SMALL)
        SMALL=10*SMALL
        SMALL=MAX(SMALL,ACCU)
        M=MAX(N,1)
        IF (MOD(M,2).NE.1) M=M+1
        ALRC=1.0D0/AMIN
        AL2RC=0.5D0*ALRC 
C
C INTEGRAL(0->INFINITY)
C
        TOTI=AL2RC
        DO I=3,M,2
         TOTI=(I-1)*AL2RC*TOTI
        END DO
        TOTI=TOTI*EPS
C
C ITERATION
C
        R=SQRT((1.0D0*(M-1)-LOG(EPS))*ALRC)
        ITER=0
  110   CONTINUE    
         ROLD=R
         ITER=ITER+1
         IF (ITER .GT. 100*M) THEN
          write(6,*)'RCUTOFF: ITERATION FAILURE'
          CALL STOPIT
         END IF
         RSQR=R*R
         FAC=EXP(-AMIN*RSQR)
         ACTI=AL2RC*FAC
         RPOW=FAC
         DO I=3,M,2
          RPOW=RPOW*RSQR
          ACTI=AL2RC*(RPOW+(I-1)*ACTI)
         END DO
         DRVI= -RPOW*R
         R=ROLD-(ACTI-TOTI)/DRVI
         IF (R .LT. 0.0D0)  R=ROLD-5.0D0
         IF (R .GT. 2*ROLD) R=ROLD+5.0D0
         IF (ABS(R-ROLD) .LE. SMALL*ABS(ROLD)) GOTO 120
         GOTO 110
  120   CONTINUE
        RCUTOFF= R
        RETURN
       END
