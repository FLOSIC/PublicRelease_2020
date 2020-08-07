C UTEP Electronic Structure Lab (2020)
C
C ************************************************************
C
       SUBROUTINE FINTPOL(NDEG,NPTS,X,RATIO,NTAB,NFSIZ,NFUSE,
     &                    XTAB,FTAB,FVAL)
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:45 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NDEG, NPTS, NTAB, NFSIZ, NFUSE, MAXSIZ, MAXDEG, I,
     & IND1, IND2, IPTS, IUSE, J, MDEG, NPOL
       REAL*8 :: SYMBOL , X, RATIO, XTAB, FTAB, FVAL, FAC, FPOL, XC,
     & XPOL
        SAVE
        PARAMETER (MAXDEG=20)
        PARAMETER (MAXSIZ=3)
C
C GIVEN A TABLE XTAB/FTAB WHERE FTAB= F(XTAB), FINTPOL CALCULATES F(X)
C BY SIMPLE POLYNOMIAL INTERPOLATION.
C XTAB MUST BE SORTED IN ASCENDING ORDER.
C IF X IS OUTSIDE OF THE INTERVAL SPANNED BY XTAB, ZERO WILL BE RETURNED.
C IF THE DISTANCE BETWEEN TWO TABULATED POINTS X1 AND X2 IS SMALLER
C THAN RATIO*(X2-X1), X2 WILL BE IGNORED (THIS IS NECESSARY IN ORDER TO 
C OBTAIN A STABLE SCHEME)
C
        DIMENSION X(NPTS),FVAL(NFSIZ,NPTS),XTAB(NTAB),FTAB(NFSIZ,NTAB)
        DIMENSION XPOL(MAXDEG),FPOL(MAXSIZ,MAXDEG)
C
        IF (NDEG .LT. 1) THEN
         PRINT *,'FINTPOL: NDEG MUST BE >= 1'
         CALL STOPIT
        END IF
        IF (NDEG .GT. MAXDEG) THEN
         PRINT *,'FINTPOL: MAXDEG MUST BE AT LEAST: ',NDEG
         CALL STOPIT
        END IF
        IF (NFUSE .GT. MAXSIZ) THEN
         PRINT *,'FINTPOL: MAXSIZ MUST BE AT LEAST: ',NFUSE
         CALL STOPIT
        END IF
        MDEG=MIN(NDEG,NTAB)
C
C LOOP OVER ALL POINTS
C FIRST, BRACKET X(IPTS)
C
        DO 100 IPTS=1,NPTS
         XC= X(IPTS)
         J= 0
         IF (XC .LT. XTAB(   1)) J=1
         IF (XC .GT. XTAB(NTAB)) J=NTAB
         IF (J .NE. 0) THEN
          DO IUSE=1,NFUSE
           FVAL(IUSE,IPTS)= FTAB(IUSE,J)
          END DO
          GOTO 90
         END IF
         DO IUSE=1,NFUSE
          FVAL(IUSE,IPTS)= 0.0D0
         END DO
         IND1=1
         IND2=NTAB
   10    CONTINUE
          IF (IND2-IND1 .LE. 1) GOTO 20
          I=(IND1+IND2)/2
          IF (XC .GT. XTAB(I)) THEN
           IND1=I
          ELSE
           IND2=I
          END IF
          GOTO 10
   20    CONTINUE
         IF (ABS(XC-XTAB(IND1)) .GT. ABS(XC-XTAB(IND2))) THEN
          IND1=IND2
         ELSE
          IND2=IND1
         END IF
         DO I=2,MDEG
          IF (IND1 .EQ.1) THEN
           IND2=IND2+1 
          ELSE IF (IND2 .EQ. NTAB) THEN
           IND1=IND1-1
          ELSE
           IF (ABS(XC-XTAB(IND1-1)) .GT. ABS(XC-XTAB(IND2+1))) THEN
            IND2=IND2+1
           ELSE
            IND1=IND1-1
           END IF
          END IF
         END DO
C
C THROW AWAY UNNECESSARY POINTS
C
         NPOL=1
         FAC= XTAB(IND1)
         XPOL(1)= FAC
         DO IUSE=1,NFUSE
          FPOL(IUSE,1)= FTAB(IUSE,IND1)
         END DO
         DO I=IND1+1,IND2
          IF ((XTAB(I)-FAC) .GT. RATIO*ABS(XC-FAC)) THEN
           NPOL=NPOL+1
           FAC= XTAB(I)
           XPOL(NPOL)=FAC
           DO IUSE=1,NFUSE
            FPOL(IUSE,NPOL)=FTAB(IUSE,I)
           END DO
          END IF
         END DO
C
C CALCULATE FUNCTION VALUE 
C          
         DO I=1,NPOL
          FAC= 1.0D0  
          DO J=1,NPOL
           IF (J .NE. I) FAC= FAC*(XC-XPOL(J))/(XPOL(I)-XPOL(J))
          END DO
          DO IUSE=1,NFUSE
           FVAL(IUSE,IPTS)= FVAL(IUSE,IPTS)+FAC*FPOL(IUSE,I)
          END DO
         END DO
   90    CONTINUE
  100   CONTINUE
        RETURN
       END
