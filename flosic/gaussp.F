C UTEP Electronic Structure Lab (2020)
C
C ***********************************************************
C
C GAUSSP CREATES A N-POINT GAUSS-LEGENDRE MESH SUITABLE FOR THE 
C INTEGRATION OF POLYNOMIALS OF DEGREE 2*N-1 IN THE INTERVAL (Y1,Y2)
C
       SUBROUTINE GAUSSP(Y1,Y2,N,Y,WY)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),WY(N)
        PARAMETER (EPS=1.0D-14)
        PARAMETER (MAX=1000)
        DIMENSION X(MAX),W(MAX)
        SAVE
        DATA N_SAV/-1/
        DATA PI/3.14159265358979323844D0/
C
        IF (N. LE. 0) RETURN
        IF (N .GT. MAX) THEN
         write(6,*)'GAUSSP: MAX MUST BE AT LEAST: ',N
         CALL STOPIT
        END IF
        IF (N .EQ. 1) THEN
         WY(1)=Y2-Y1
         Y(1)=0.5D0*(Y2+Y1)
         RETURN
        END IF
C
C FIRST TRY PRESTORED MESHES
C IF NONE OF THEM MATCHES, DEFINE MESH VIA ITERATION OF CHEBYSHEV POINTS
C POINTS/WEIGHTS WILL BE DEFINED FOR INTERVAL (-1:1)
C
        IF (N .NE. N_SAV) THEN
         N_SAV=N
         M=N
         CALL GAUSSPP(M,X,W)
         IF (M .NE. N) THEN
          M=(N+1)/2
          DO I=1,M
           Y(I)=COS(PI*(I-0.25D0)/(N+0.5D0))
          END DO
          DO 100 I=1,M
           Z=Y(I)
  20       CONTINUE
            P1=1.0D0
            P2=0.0D0
            DO J=1,N
             P3=P2
             P2=P1
             P1=((2*J-1)*Z*P2-(J-1)*P3)/J
            END DO
            PP=N*(Z*P1-P2)/(Z*Z-1.0D0)
            Z1=Z
            Z=Z1-P1/PP
            IF (ABS(Z-Z1).GT.EPS) GOTO 20
           CONTINUE
           X(I)    = -Z
           X(N+1-I)=  Z
           W(I)=2.0D0/((1.0D0-Z*Z)*PP*PP)
           W(N+1-I)=W(I)
 100      CONTINUE
         END IF
        END IF
C
C SCALE WEIGHTS AND POINTS:
C
        DO I=1,N
         FACT= 0.5D0*(Y2-Y1)
         Y (I)= Y1+FACT*(X(I)+1.0D0)
         WY(I)= FACT*W(I)
        END DO
        RETURN
        END
