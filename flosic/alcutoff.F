C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       FUNCTION ALCUTOFF(N,RSQR,EPS)
C
C ALCUTOFF DETERMINES THE VALUE ALP SO THAT
C
C  INF  N            2          INF  N            2
C  INT R * EXP(-ALP*R ) - EPS * INT R * EXP(-ALP*R ) = 0
C   R                            0
C
C INPUT: N, RSQR= R*R, EPS (SEE ABOVE FORMULA)
C
        IMPLICIT REAL*8 (A-H,O-Z)
        SAVE
        DATA TOLER/1.0D-14/       
C
C SIMPLE CASES
C
        IF ((N .LT. 0) .OR. (EPS .GE. 1.0D0)) THEN
         ALCUTOFF= 0.0D0
         RETURN
        ELSE IF ((RSQR .LE. 0.0D0) .OR. (EPS .LE. 0.0D0)) THEN
         ALCUTOFF= 1.0D30
         RETURN
        END IF
C
C INITIALIZATION: SOLUTION FOR M=1
C THEN ITERATE
C
        M=N
        IF (MOD(M,2) .NE. 1) M=M+1
        ALP= -LOG(EPS)/RSQR
   10   CONTINUE
         OLD= ALP
         ALRCH= 0.5D0/ALP
         FAC=  ALRCH*EXP(-ALP*RSQR)
         SUM1= FAC
         SUM2= ALRCH
         DO I=3,M,2
          FAC= FAC*RSQR
          FF= (I-1)*ALRCH
          SUM1= FAC+FF*SUM1
          SUM2=     FF*SUM2
         END DO
         FAC= FAC*RSQR
         FF= (M+1)*ALRCH
         DR1= FAC+FF*SUM1
         DR2=     FF*SUM2
         ANEW= ALP+(SUM1-EPS*SUM2)/(DR1-EPS*DR2)
         IF (ANEW .GT. 2.0D0*ALP) ANEW= 2.0D0*ALP
         IF (ANEW .LT. 0.5D0*ALP) ANEW= 0.5D0*ALP
         ALP= ANEW
         IF (ABS(ALP-OLD) .LT. TOLER*OLD) GOTO 20
         GOTO 10
   20   CONTINUE
        ALCUTOFF= ALP
        RETURN
       END
