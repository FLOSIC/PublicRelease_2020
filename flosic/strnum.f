C UTEP Electronic Structure Lab (2020)
      CHARACTER*(*) FUNCTION STRNUM(STRING,MATRIX,MDIM,NDIM,M,N)
C
C     EXTEND A STRING WITH A NUMBER.
C
C     ------------------------------------------------------------------
C
C     LOCAL DIMENSIONS:
C
C     MDIM: Row dimension of MATRIX in the calling subroutine.
C     NDIM: Column dimension of MATRIX in the calling subroutine.
C
C     LOCAL VARIABLES:
C
C     COL   : Column index.
C     M     : Number of rows.
C     MATRIX: Matrix of already defined strings.
C     N     : Number of columns.
C     ROW   : Row index.
C     STRING: String to be extended.
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER MDIM,NDIM
      CHARACTER*(*) MATRIX(MDIM,NDIM)
C
      CHARACTER*(*) STRING
C
      CHARACTER*9 LABEL
      CHARACTER*10 MATSTR
      CHARACTER*80 PRTSTR,STRCOMP
      INTEGER COL,LENLAB,LENSTR,M,N,NUMBER,STREXT,ROW,APPEAR
C
C     ------------------------------------------------------------------
C
C     --- INITIALIZATION ---
C
      APPEAR = 0
      NUMBER = 0
      WRITE (LABEL,'(I9)') NUMBER
C
      IF (STRING(1:1).EQ.'-') THEN
        STRNUM = STRCOMP(STRING(2:))
      ELSE
        STRNUM = STRCOMP(STRING)
      END IF
C
C     --- LOOP OVER ALREADY DEFINED STRINGS ---
C
  500 CONTINUE
C
      LENSTR = STREXT(STRNUM) + 1
C
      DO ROW=1,M
        DO COL=1,N
C
          IF (MATRIX(ROW,COL)(1:1).EQ.'-') THEN
            MATSTR = MATRIX(ROW,COL)(2:LENSTR+1)
          ELSE
            MATSTR = MATRIX(ROW,COL)(1:LENSTR)    
          END IF
C
          IF (INDEX(MATSTR,STRNUM).NE.0) APPEAR = APPEAR + 1
C
          IF (APPEAR.GT.1) THEN
C
            APPEAR = 1
            NUMBER = NUMBER + 1
            WRITE (LABEL,'(I9)') NUMBER
C
            IF (STRING(1:1).EQ.'-') THEN
              STRNUM = STRCOMP(STRING(2:))
            ELSE
              STRNUM = STRCOMP(STRING)
            END IF
C
            LABEL = STRCOMP(LABEL)
            LENSTR = STREXT(STRNUM)
            LENLAB = STREXT(LABEL)
            STRNUM = STRNUM(1:LENSTR)//LABEL(1:LENLAB)
C
            IF (STREXT(STRNUM).GT.8) THEN
              WRITE(*,*)'MESSAGE IN STRNUM: INPUT ERROR DETECTED'
              WRITE (PRTSTR,5000) STRNUM(1:LENSTR) 
 5000         FORMAT ('EQUIVALENT VARIABLE < ',A,' > IS TO LARGE')
              WRITE(*,*)'ERROR IN STRNUM: ',STRCOMP(PRTSTR)
            END IF 
C
            GO TO 500
C
          END IF
C
        END DO                  
      END DO
C
C     ------------------------------------------------------------------
C
      END
