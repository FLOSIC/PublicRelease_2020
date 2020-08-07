C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE ATMFIT(ALMIN,NALP,NPTS,RPTS,WPTS,FPTS,
     &                   MATSZ,FCOF,AMAT,AVEC,IVEC)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RPTS(NPTS),WPTS(NPTS),FPTS(NPTS)
        DIMENSION FCOF(MATSZ),AMAT(MATSZ,MATSZ),AVEC(MATSZ),IVEC(MATSZ)
        SAVE
C
C LOW-ACCURACY LEAST-SQUARES FIT OF FUNCTION F(R) GIVEN ON A MESH 
C OF POINTS FPTS(RPTS) WITH WEIGHT WPTS TO THE FUNCTIONAL FORM:
C
C       NALP                        2                   I-1
C F(R)= SUM  FCOF(I) * EXP(-ALP(I)*R )   WHERE ALP(I)= 2   * ALMIN
C       I=1
C
        IF (NALP .LT. 1) RETURN
        IF (NALP .GT. MATSZ) THEN
         PRINT *,'ATMFIT: MATSZ MUST BE AT LEAST: ',NALP
         CALL STOPIT
        END IF
        DO I= 1,NALP
         FCOF(I)= 0.0D0
         DO J= I,NALP
          AMAT(J,I)= 0.0D0
         END DO
        END DO
C
C SETUP AMAT, AVEC
C
        DO IPTS= 1,NPTS
         RR= RPTS(IPTS)
         WR= WPTS(IPTS)
         FR= FPTS(IPTS)
         R2= RR*RR
         REX1= EXP(-ALMIN*R2)
         DO I= 1,NALP
          FAC1= WR*REX1
          FCOF(I)= FCOF(I)+FAC1*FR
          REX2= REX1
          DO J= I,NALP
           AMAT(J,I)= AMAT(J,I)+FAC1*REX2
           REX2= REX2*REX2
          END DO
          REX1= REX1*REX1
         END DO
        END DO
        DO I= 1,NALP
         DO J= I,NALP
          AMAT(I,J)= AMAT(J,I)
         END DO
        END DO
C
C GET SOLUTION
C
        CALL GAUSSPIV(MATSZ,1,NALP,1,AMAT,FCOF,AVEC,IVEC,IER)
        IF (IER .NE. 0) THEN
         PRINT *,'ATMFIT: ERROR IN GAUSSPIV'
         CALL STOPIT
        END IF
        RETURN
       END 
