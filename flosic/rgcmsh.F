C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE RGCMSH(IMSH,NPTS,RST,RMID,RRET,WRET)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RRET(*),WRET(*)
        SAVE
C
C WRITTEN BY DV POREZAG 
C THIS SUBROUTINE RETURNS A ONE DIMENSIONAL MESH SUITABLE FOR
C INTEGRATING FUNCTIONS OF THE FORM:
C
C  INF
C  INT  DR * R**2 * F(R)
C  RST
C
C METHOD: 
C IF (IMSH .EQ. 1) WE USE A GAUSS-CHEBYSHEV QUADRATURE OF
C THE SECOND KIND:
C
C       2     /                    2              3 \
C  X= ---- * | 3*ASIN(Z) + SQRT(1-Z ) * (5*Z - 2*Z ) |
C     3*PI    \                                     /
C
C RESULTING IN POINTS AND WEIGHTS:
C
C      N+1-2*I    2    /    2    2 I*PI \        I*PI       I*PI
C  X = ------- + -- * | 1 + - SIN  ----  | * SIN ---- * COS ----
C   I    N+1     PI    \    3       N+1 /         N+1        N+1
C
C         16        4 I*PI
C  W = ------- * SIN  ----       WHERE (1 <= I <= N)
C   I  3*(N+1)         N+1
C
C IF (IMSH .EQ. 2) WE USE GAUSS-LEGENDRE POINTS INSTEAD.
C THE MESH IS TRANSFORMED VIA:
C
C     RMID                               / RMID * LOG(2) \
C R= ------ * LOG(B*(1-X))  WHERE B= EXP|  -------------  |
C    LOG(B)                              \   RST - RMID  /
C
C ********************************************************************
C
        IF (NPTS .LT. 1) RETURN
        IF (RST .LT. 0.0D0) THEN
         write(6,*)'RGCMSH: RST MUST BE >= 0'
         CALL STOPIT
        END IF
        IF (RMID .LE. RST) THEN
         write(6,*)'RGCMSH: RMID MUST BE > RST'
         CALL STOPIT
        END IF
C
C TRANSFORMED GAUSS-CHEBYSHEV QUADRATURE OF THE SECOND KIND
C
        IF (IMSH .EQ. 1) THEN
         PI= 4*ATAN(1.0D0)
         N2= (NPTS+1)/2
         REC= 1.0D0/(NPTS+1)
         PREC= PI*REC
         FAC= 2.0D0/PI
         TWTH= 2.0D0/3.0D0
         WFAC= 16.0D0/3.0D0*REC
         DO I= 1,N2
          X= (NPTS+1-I)*PREC
          S= SIN(X)
          C= COS(X)
          S2= S*S
          RRET(I)= (2*I-1-NPTS)*REC+FAC*(1.0D0+TWTH*S2)*S*C
          WRET(I)= WFAC*S2*S2
          RRET(NPTS-I+1)= -RRET(I)
          WRET(NPTS-I+1)=  WRET(I)
         END DO
        ELSE
         CALL GAUSSP(-1.0D0,1.0D0,NPTS,RRET,WRET)
        ENDIF
C
C RADIAL TRANSFORMATION
C
        B= EXP(RMID*LOG(2.0D0)/(RST-RMID))
        FAC= RMID/LOG(B)
        DO I= 1,NPTS
         OMX= 1.0D0-RRET(I)
         RRET(I)= FAC*LOG(B*OMX)
         WRET(I)= -FAC*WRET(I)/OMX
         WRET(I)= WRET(I)*RRET(I)**2
        END DO
        RETURN
       END
