C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       SUBROUTINE TABDRV(NDEG,NDRV,RATIO,NTAB,XTAB,FTAB,NSIZ,FOUT)
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NDEG, NDRV, NTAB, NSIZ, MAXDEG, I, IND1, IND2, IPTS,
     & ISTT, J, K, MDEG, NLEFT, NPOL
       REAL*8 :: SYMBOL , RATIO, XTAB, FTAB, FOUT, DIST, FAC, FC, FPOL,
     & PRD, SUM1, SUM2, XC, XPOL
        SAVE
        PARAMETER (MAXDEG=30)
C
C GIVEN A TABLE XTAB/FTAB WHERE FTAB= F(XTAB), TABDRV CALCULATES:
C   (IF NDRV.EQ.0): F(X)
C   (IF NDRV.EQ.1): F(X) AND dF(X)/dX 
C   (IF NDRV.EQ.2): F(X), dF(X)/dX, AND d2F(X)/dX**2
C FOR EVERY POINT XTAB(I), I=1,NTAB USING A SIMPLE POLYNOMIAL 
C INTERPOLATION SCHEME OF DEGREE NDEG.
C IF THE DISTANCE BETWEEN TWO TABULATED POINTS X1 AND X2 IS SMALLER
C THAN RATIO*MAX(ABS(X-X1),ABS(X-X2)) WHERE X IS THE POINT FOR WHICH 
C THE INTERPOLATION IS DONE, X2 WILL BE IGNORED (THIS IS NECESSARY IN 
C ORDER TO OBTAIN A STABLE SCHEME)
C
        DIMENSION XTAB(NTAB),FTAB(NTAB),FOUT(NSIZ,NTAB)
        DIMENSION XPOL(MAXDEG),FPOL(MAXDEG)
C
        IF ((NDRV .LT. 0) .OR. (NDRV .GT. 2)) THEN
         PRINT *,'TABDRV: NDRV MUST BE 0, 1, OR 2'
         CALL STOPIT
        END IF
        IF (NSIZ .LT. NDRV) THEN
         PRINT *,'TABDRV: NSIZ MUST BE >= NDRV'
         CALL STOPIT
        END IF
        IF (NDEG .GT. MAXDEG) THEN
         PRINT *,'TABDRV: MAXDEG MUST BE AT LEAST: ',NDEG
         CALL STOPIT
        END IF
        MDEG=MIN(NDEG,NTAB-1)
        IF (MDEG .LT. 2) THEN
         PRINT *,'TABDRV: MDEG MUST BE >= 2'
         CALL STOPIT
        END IF
C
C LOOP OVER ALL POINTS
C LOOK FOR NDEG CLOSEST POINTS
C
        DO 100 IPTS= 1,NTAB 
         XC= XTAB(IPTS)
         FC= FTAB(IPTS)
         FOUT(1,IPTS)= FC
         IF (NDRV .EQ. 0) GOTO 90
         NLEFT= NDEG/2
         IND1= IPTS-NLEFT
         IND2= IND1+NDEG
         IF (IND1 .LT. 1) THEN
          IND1= 1
          IND2= NDEG+1
         END IF
         IF (IND2 .GT. NTAB) THEN
          IND2= NTAB
          IND1= NTAB-NDEG
         END IF
C
C THROW AWAY UNUSABLE POINTS
C
         ISTT= IND1
         IF (IPTS .EQ. ISTT) ISTT= IND1+1
         NPOL= 1
         XPOL(1)= XTAB(ISTT)
         FPOL(1)= FTAB(ISTT)
         DO I= ISTT+1,IND2
          DIST= MAX(ABS(XC-XPOL(NPOL)),ABS(XC-XTAB(I)))
          IF ((I .NE. IPTS) .AND. 
     &        (XTAB(I)-XPOL(NPOL) .GT. RATIO*DIST)) THEN
           NPOL= NPOL+1
           XPOL(NPOL)= XTAB(I)
           FPOL(NPOL)= FTAB(I)
          END IF
         END DO
C
C CALCULATE DERIVATIVES
C          
         FOUT(2,IPTS)= 0.0D0
         IF (NDRV .EQ. 1) THEN
          DO I= 1,NPOL
           PRD= 1.0D0
           DO K= 1,NPOL
            IF (K .NE. I) PRD= PRD*(XC-XPOL(K))/(XPOL(I)-XPOL(K))
           END DO
           FOUT(2,IPTS)= FOUT(2,IPTS)+(FC-FPOL(I)*PRD)/(XC-XPOL(I))
          END DO
         ELSE
          FOUT(3,IPTS)= 0.0D0
          DO I= 1,NPOL
           SUM1= 0.0D0
           SUM2= 0.0D0
           DO J= 1,NPOL
            IF (J .NE. I) THEN
             SUM1= SUM1+1.0D0/(XC-XPOL(J))
             PRD= 1.0D0
             DO K= 1,NPOL
              IF ((K .NE. I) .AND. (K .NE. J)) THEN
                PRD= PRD*(XC-XPOL(K))/(XPOL(I)-XPOL(K))
              END IF
             END DO
             SUM2= SUM2+PRD/(XPOL(I)-XPOL(J))
            END IF
           END DO
           PRD= 1.0D0
           DO K= 1,NPOL
            IF (K .NE. I) PRD= PRD*(XC-XPOL(K))/(XPOL(I)-XPOL(K))
           END DO
           FAC= 1.0D0/(XC-XPOL(I))
           FOUT(2,IPTS)= FOUT(2,IPTS)+(FC-FPOL(I)*PRD)*FAC
           FOUT(3,IPTS)= FOUT(3,IPTS)+(FC*SUM1-2*FPOL(I)*SUM2)*FAC
          END DO
         END IF
   90    CONTINUE
  100   CONTINUE
        RETURN
       END
