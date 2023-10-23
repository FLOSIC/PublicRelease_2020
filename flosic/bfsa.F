C UTEP Electronic Structure Lab (2020)
C
C ***********************************************************
C
       SUBROUTINE BFSA(MAXPTS,LMAX,NPTS,ANGLE,WTS)
C WRITTEN BY MARK R PEDERSON  (1992-1994)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(MAXANG=500)
        DIMENSION X(MAXANG),Z(MAXANG),WTX(MAXANG),WTZ(MAXANG)
        DIMENSION ANGLE(3,MAXPTS),WTS(MAXPTS)
        SAVE
C
C SETUP SOME STUFF AND CHECK INPUT
C
        IF (LMAX.LT.0) RETURN
        N1= (LMAX+2)/2
        LMAX1=LMAX+1
        NPTS=N1*LMAX1
        IF (LMAX1 .GT. MAXANG) THEN
         write(6,*)'BFSA: MAXANG MUST BE AT LEAST ',LMAX1
         CALL STOPIT
        END IF
        IF (NPTS.GT.MAXPTS) THEN
         write(6,*)'BFSA: SIZE OF ARRAY ANGLE EXCEEDED'
         CALL STOPIT
        END IF
C
C THETA MESH
C
        DN= -1.0D0
        UP= 1.0D0
        CALL GAUSSP(DN,UP,N1,Z,WTZ)
        DO I=1,N1
         Z(I)= ACOS(Z(I))
        END DO
C
C PHI MESH
C
        PI=4*ATAN(1.0D0)
        ALMAX=LMAX1
        DPHI=(2*PI)/ALMAX
        PHI= -0.5D0*DPHI
        DO I=1,LMAX1
         PHI=PHI+DPHI
         X(I)=PHI
         WTX(I)=DPHI
        END DO
C
C SPHERE MESH
C
        NPTS=0
        DO I=1,LMAX1
         DO J=1,N1
          NPTS=NPTS+1
          ANGLE(1,NPTS)=SIN(Z(J))*SIN(X(I))
          ANGLE(2,NPTS)=SIN(Z(J))*COS(X(I))
          ANGLE(3,NPTS)=COS(Z(J))
          WTS(NPTS)=WTX(I)*WTZ(J)
         END DO
        END DO
        RETURN
       END
