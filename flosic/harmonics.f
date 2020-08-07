C UTEP Electronic Structure Lab (2020)
C
C **********************************************************************
C
       SUBROUTINE HARMONICS(MXPTS,NPTS,LM,RXYZ,YLM,NYLM)
C
C SPHERICAL HARMONICS GENERATOR BASED ON RECURRENCE FORMULAE
C WRITTEN BY DIRK POREZAG
C
C MXPTS: USED TO DEFINE DIMENSION OF RXYZ, YLM
C NPTS:  ACTUAL NUMBER OF POINTS
C LM:    LARGEST ANGULAR MOMENTUM REQUIRED
C RXYZ:  COORDINATES OF POINTS (normalized)
C YLM:   SPHERICAL HARMONICS AS CALCULATED ON POINTS
C NYLM:  TOTAL NUMBER OF SPHERICAL HARMONICS PER POINT
C 
C ATTENTION: THIS ROUTINE EXPECTS NORMALIZED POINTS IN ARRAY RXYZ
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RXYZ(3,MXPTS),YLM(MXPTS,(LM+1)**2)
        DIMENSION PREFACRC(50)
        LOGICAL FIRST
        SAVE
        DATA FIRST/.TRUE./

        IF (LM .GT. 50) THEN
         PRINT *,'HARMONICS: THIS ROUTINE IS NOT STABLE FOR L > 50'
         CALL STOPIT
        END IF
        NYLM=0
        IF (LM .LT. 0) RETURN
        NYLM=(LM+1)**2
        IF (NPTS .LE. 0) RETURN
C
C SET UP 1/I ARRAY
C
        IF (FIRST) THEN
         FIRST= .FALSE.
         DO I=1,50 
          PREFACRC(I)= 1.0D0/I
         END DO
        END IF
C
C L=0
C
        DO IPTS=1,NPTS
         YLM(IPTS,1)= 1.0D0
        END DO
        IF (LM .LT. 1) GOTO 200
C
C L=1
C
        DO IPTS=1,NPTS
         YLM(IPTS,2)= RXYZ(3,IPTS)
         YLM(IPTS,3)= RXYZ(1,IPTS)
         YLM(IPTS,4)= RXYZ(2,IPTS)
        END DO
C
C HIGHER L: START WITH RECURSION FOR Y_LM (M <= L-1)
C
        DO 100 L=2,LM
         IFC2L1=2*L-1
         NOW=L**2
         IL1=(L-1)**2
         IL2=(L-2)**2
         DO M=0,L-1
          FACRC=PREFACRC(L-M)
          IFCLM1=L+M-1
          IF (M .EQ. L-1) IFCLM1=0 
          AA= FACRC*IFC2L1
          BB= FACRC*IFCLM1
          NRUN=2
          IF (M .EQ. 0) NRUN=1
          DO IRUN=1,NRUN
           NOW=NOW+1
           IL1=IL1+1
           IL2=IL2+1
           DO IPTS=1,NPTS
            YLM(IPTS,NOW)= AA*YLM(IPTS,IL1)*YLM(IPTS,2)
     &                    -BB*YLM(IPTS,IL2)
           END DO
          END DO
         END DO
C
C RECURSION FOR Y_LL
C
         LAST=L**2
         NOW=(L+1)**2
         DO IPTS=1,NPTS
          YLM(IPTS,NOW)=   IFC2L1*(YLM(IPTS,LAST  )*YLM(IPTS,3)
     &                            +YLM(IPTS,LAST-1)*YLM(IPTS,4))
          YLM(IPTS,NOW-1)= IFC2L1*(YLM(IPTS,LAST-1)*YLM(IPTS,3)
     &                            -YLM(IPTS,LAST  )*YLM(IPTS,4))
         END DO
  100   CONTINUE
C
C PREFACTORS
C
  200   FOURPI=16*ATAN(1.0D0)
        ONOFPI=1.0D0/FOURPI
        DO L=0,LM
         FACFC= (2*L+1)*ONOFPI
         NOW=L**2
         DO M=0,L
          NRUN=2
          IF (M .EQ. 0) NRUN=1
          VFAC=SQRT(NRUN*FACFC)
          DO IRUN=1,NRUN
           NOW=NOW+1
           DO IPTS=1,NPTS
            YLM(IPTS,NOW)= VFAC*YLM(IPTS,NOW)
           END DO
          END DO
          FACFC= FACFC/(MAX(L-M,1)*(L+M+1))
         END DO
        END DO
        RETURN
       END
