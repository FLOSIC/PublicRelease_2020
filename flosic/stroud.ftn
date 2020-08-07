C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       SUBROUTINE STROUD(MAXPTS,LMAX,NPTS,ANGLE,DOMEGA)
C ORIGINAL VERSION BY MARK R PEDERSON (1989)
        IMPLICIT REAL*8 (A-H,O-Z)
C
C SOME REALLY GOOD MESHES OUT OF:
C
C A.H. STROUD "APPROXIMATE CALCULATION OF MULTIPLE INTEGRALS"
C              PRENTICE-HALL SERIES IN AUTOMATIC COMPUTATION
C SEE PAGES 294-299
C
C U3:  5-1 12 PT MESH OF DEGREE  5 
C U3:  9-1 32 PT MESH OF DEGREE  9 
C U3: 14-1 72 PT MESH OF DEGREE 14
C
        DIMENSION ANGLE(3,MAXPTS),DOMEGA(MAXPTS)
        DIMENSION XS14(3,5)
        SAVE
        DATA XS14/ -0.15110827466456159D0, 0.15524060002531609D0,
     &              0.97625132288344653D0,
     &              0.31583835294225182D0, 0.25704938737035939D0,
     &              0.91333003195079960D0,
     &              0.34630711173365857D0, 0.66627779003666489D0,
     &              0.66041296994119471D0,
     &             -0.10180878734504331D0, 0.81738606470122647D0,
     &              0.56702292021714103D0,
     &             -0.40922840266630550D0, 0.50154771175897465D0,
     &              0.76222175728623802D0/
C
C CHECK SIZES
C
        PI=4*ATAN(1.0D0)
        NPTS=0
        IF(LMAX.LE.5)THEN 
         NPTS=12
        ELSE IF(LMAX.LE.9)THEN
         NPTS=32
        ELSE IF(LMAX.LE.14)THEN
         NPTS=72
        END IF
        IF (NPTS.GT.MAXPTS) THEN
         write(6,*)'STROUD: SIZE OF ARRAY ANGLE EXCEEDED'
c         call errorfile('STROUD: SIZE OF ARRAY ANGLE EXCEEDED')
         CALL STOPIT
        END IF
        NPTS=0
C
C MESH OF DEGREE 5 
C
        IF(LMAX.LE.5)THEN 
         NPTS=0
         WTMP=(1.0D0/3.0D0)*PI
         ANG1= 0.85065080835203988D0
         ANG2= 0.52573111211913359D0
         DO 30 ISW=1,3
          K1= 1+MOD(ISW-1,3)
          K2= 1+MOD(ISW  ,3)
          K3= 1+MOD(ISW+1,3)
          DO 20 I1=1,2
           ANG1= -ANG1
           DO 10 I2=1,2
            NPTS=NPTS+1
            ANG2= -ANG2
            ANGLE(K1,NPTS)= ANG1
            ANGLE(K2,NPTS)= ANG2
            ANGLE(K3,NPTS)= 0.0D0
            DOMEGA(NPTS)=WTMP
   10      CONTINUE
   20     CONTINUE
   30    CONTINUE
C
C MESH OF DEGREE 9
C
        ELSE IF(LMAX.LE.9)THEN 
         NPTS=0
         ANG1=SQRT(0.1D0*(5.0D0+SQRT(5.0D0)))
         ANG2=SQRT(0.1D0*(5.0D0-SQRT(5.0D0)))
         WTMP=(5.0D0/42.0D0)*PI
         DO 140 IMESH=1,2
          DO 130 ISW=1,3
           K1= 1+MOD(ISW-1,3)
           K2= 1+MOD(ISW  ,3)
           K3= 1+MOD(ISW+1,3)
           DO 120 I1=1,2
            ANG1= -ANG1
            DO 110 I2=1,2
             NPTS=NPTS+1
             ANG2= -ANG2
             ANGLE(K1,NPTS)= ANG1
             ANGLE(K2,NPTS)= ANG2
             ANGLE(K3,NPTS)= 0.0D0
             DOMEGA(NPTS)=WTMP
  110       CONTINUE
  120      CONTINUE
  130     CONTINUE
          ANG1=SQRT((3.0D0-SQRT(5.0D0))/6.0D0)
          ANG2=SQRT((3.0D0+SQRT(5.0D0))/6.0D0)
          WTMP=(9.0D0/70.0D0)*PI
  140    CONTINUE
C
         XTMP=1.0D0/SQRT(3.0D0)
         SIGNX=XTMP
         SIGNY=XTMP
         SIGNZ=XTMP
         DO 190 IX=1,2
          SIGNX= -SIGNX
          DO 180 IY=1,2
           SIGNY= -SIGNY
           DO 170 IZ=1,2
            NPTS=NPTS+1
            SIGNZ= -SIGNZ
            ANGLE(1,NPTS)=SIGNX
            ANGLE(2,NPTS)=SIGNY
            ANGLE(3,NPTS)=SIGNZ
            DOMEGA(NPTS)=WTMP
  170      CONTINUE
  180     CONTINUE
  190    CONTINUE
C
C MESH OF DEGREE 14
C
        ELSE IF(LMAX.LE.14)THEN 
         NPTS=0
         X1=0.52573111211913360D0
         X2=0.85065080835203990D0
         X3=0.0D0
         WTMP=(25.0D0/504.0D0)*PI
         DO 230 J=1,2
          X1= -X1
          DO 220 K=1,2
           X2= -X2
           DO 210 L=1,3
            XT=X1
            X1=X3
            X3=X2
            X2=XT
            NPTS=NPTS+1
            ANGLE(1,NPTS)=X1
            ANGLE(2,NPTS)=X2
            ANGLE(3,NPTS)=X3
            DOMEGA (NPTS)=WTMP
  210      CONTINUE
  220     CONTINUE
  230    CONTINUE
C
         WTMP=(143.0D0/2520.0D0)*PI
         DO 260 N=1,5
          X1=XS14(1,N)
          X2=XS14(2,N)
          X3=XS14(3,N)
          DO 250 I=1,3
           XT= X1
           X1= X3
           X3= -X2
           X2= -XT
           DO 240 J=1,3
            XT=X1
            X1=X3
            X3=X2
            X2=XT
            NPTS=NPTS+1
            ANGLE(1,NPTS)=X1
            ANGLE(2,NPTS)=X2
            ANGLE(3,NPTS)=X3
            DOMEGA (NPTS)=WTMP
  240      CONTINUE
           X2= -X2
           X3= -X3
           NPTS=NPTS+1
           ANGLE(1,NPTS)=X1
           ANGLE(2,NPTS)=X2
           ANGLE(3,NPTS)=X3
           DOMEGA (NPTS)=WTMP
  250     CONTINUE
  260    CONTINUE
        ELSE
         NPTS=0
        END IF
        RETURN
       END
