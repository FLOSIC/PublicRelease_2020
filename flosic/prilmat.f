C UTEP Electronic Structure Lab (2020)
      SUBROUTINE PRILMAT(MATRIX,DMAT,M1,M2,N1,N2,BEFORE,AFTER,TAPE,
     $                   OPTION,STORED,TEXT)
C
C     PRINT MATRIX WITH ATOMIC ORBITAL LABELS ON TAPE.
C
C     BY ULISES REVELES, DEC. 2013.
C
C     ------------------------------------------------------------------
C
C     --- GLOBAL VARIABLES ---
C
C     ATOMORB_LAB: Full atomic orbital labels.
C
C     --- LOCAL DIMENSIONS ---
C
C     DMAT: Dimension of MATRIX in the calling subroutine.
C
C     --- LOCAL VARIABLES ---
C
C     AFTER : Number of digits after the decimal point (AFTER.GE.0).
C     BEFORE: Number of digits before the decimal point (BEFORE.GE.1).
C     MATRIX: Matrix to be printed.
C     M1    : First row of MATRIX (M.LE.DMAT) to be printed.
C     M2    : Last row of MATRIX (M.LE.DMAT) to be printed.
C     N1    : First column of MATRIX (N.LE.DMAT) to be printed.
C     N2    : Last column of MATRIX (N.LE.DMAT) to be printed.
C     TAPE  : Number of output tape (TAPE.GE.0).
C     TEXT  : Text for matrix header.
C     OPTION: a) OPTION = NONE, print no labels.
C             b) OPTION = BASIS, print orbital labels.
C             c) OPTION = TBBASIS, print tight-binding orbital labels.
C     STORED: a) STORED = TOTAL, complete matrix form.
C             b) STORED = LOWDIA, lower diagonal matrix form.
C             c) STORED = UPPDIA, upper diagonal matrix form.
C
C     ------------------------------------------------------------------
C
      use xmol,only: ATOMORB_LAB
C
      INTEGER DMAT
      REAL*8 MATRIX(DMAT,*)
C
      CHARACTER*(*) OPTION,STORED,TEXT
      INTEGER AFTER,BEFORE,COLUMN,FROM,I,J,K,L,LEFT,M1,M2,N1,N2,RIGHT,
     $        STREXT,TAPE,TEXLEN,TO,WIDTH
C
      CHARACTER FORM1*24,FORM2*24,FORM3*25,FORM4*38
C
C     ------------------------------------------------------------------
C
C     --- DEFINITION OF VARIABLE FORMATS ---
C
      FORM1 = '(/,T2,7X,??(??X,I5,??X))'
      FORM2 = '(T2,I5,2X,??(1X,F??.??))'
      FORM3 = '(/,T2,17X,??(??X,I5,??X))'
      FORM4 = '(T2,I4,1X,I3,2X,A2,A5, ??(1X,F??.??))'
C
C     ------------------------------------------------------------------
C
C     --- WRITE HEADER ---
C
      IF ((TAPE.GT.0).AND.(TAPE.LT.100)) THEN
        TEXLEN = STREXT(TEXT)
        IF (TEXLEN.GT.0) THEN
          WRITE (*,5000) TEXT(:TEXLEN)
C         WRITE (TAPE,5000) TEXT(:TEXLEN)
 5000     FORMAT (/,/,T2,A,/)
        ELSE
          WRITE (*,'(3(/))')
C         WRITE (TAPE,'(3(/))')
        END IF
      ELSE
        GO TO 99970
      END IF
C
C     --- CHECK PARAMETERS M1, M2, N1, N2, AND BEFORE ---
C
      IF ((M1.LT.1).OR.(M1.GT.DMAT).OR.(M2.LT.1).OR.(M2.GT.DMAT)) THEN
        GO TO 99960
      END IF
C
      IF ((N1.LT.1).OR.(N2.LT.1)) THEN
        GO TO 99950
      END IF
C
      IF ((BEFORE.LT.1).OR.(BEFORE.GT.15)) THEN
        GO TO 99940
      END IF
C
C     --- WRITE THE VARIABLE FORMAT ---
C
      WIDTH = BEFORE + AFTER + 3
      IF (OPTION.EQ.'NONE') THEN
        COLUMN = INT(72/WIDTH)
      ELSE IF (OPTION.EQ.'BASIS') THEN
        COLUMN = INT(56/WIDTH)
      ELSE IF (OPTION.EQ.'TBBASIS') THEN
        COLUMN = INT(56/WIDTH)
      END IF
      WRITE (FORM1(10:11),5010) COLUMN
 5010 FORMAT (I2)
      WRITE (FORM3(11:12),5010) COLUMN
      RIGHT = MAX((AFTER-2),1)
      LEFT =  WIDTH - RIGHT - 5
      WRITE (FORM1(13:14),5010) LEFT
      WRITE (FORM3(14:15),5010) LEFT
      WRITE (FORM1(20:21),5010) RIGHT
      WRITE (FORM3(21:22),5010) RIGHT
      IF ((AFTER.GE.0).AND.(AFTER.LT.16)) THEN
        WRITE (FORM2(11:12),5010) COLUMN
        WRITE (FORM4(24:25),5010) COLUMN
        WRITE (FORM2(18:19),5010) WIDTH - 1
        WRITE (FORM4(31:32),5010) WIDTH - 1
        WRITE (FORM2(21:22),5010) AFTER
        WRITE (FORM4(34:35),5010) AFTER
      ELSE
        GO TO 99930
      END IF
C
C     --- WRITE MATRIX IN FORMAT M x N ---
C
      DO 20 I=N1,N2,COLUMN
        FROM = I
        TO = MIN((FROM+COLUMN-1),N2)
        IF (OPTION.EQ.'NONE') THEN
          WRITE (*,FORM1) (J,J=FROM,TO)
        ELSE IF (OPTION.EQ.'BASIS') THEN
          WRITE (*,FORM3) (J,J=FROM,TO)
        ELSE IF (OPTION.EQ.'TBBASIS') THEN
          WRITE (*,FORM3) (J,J=FROM,TO)
        END IF
C
        DO 10 J=M1,M2
C
          IF (OPTION.EQ.'NONE') THEN 
C
            IF (STORED.EQ.'LOWDIA') THEN
              WRITE (*,FORM2) J,(MATRIX(MAX(J,K),MIN(J,K)),K=FROM,TO)
            ELSE IF (STORED.EQ.'UPPDIA') THEN
              WRITE (*,FORM2) J,(MATRIX(MIN(J,K),MAX(J,K)),K=FROM,TO)
            ELSE IF (STORED.EQ.'TOTAL') THEN
              WRITE (*,FORM2) J,(MATRIX(J,K),K=FROM,TO)
            ELSE
              GO TO 99920
            END IF
C
          ELSE IF (OPTION.EQ.'BASIS') THEN
            IF ((J.NE.1).AND.(ATOMORB_LAB(J)%ORBLAB.EQ.'1S ')) 
     %        WRITE (*,*)
            IF (STORED.EQ.'LOWDIA') THEN
              WRITE (*,FORM4) J,ATOMORB_LAB(J)%ATNUM,
     %                          ATOMORB_LAB(J)%ATLAB,
     %                          ATOMORB_LAB(J)%ORBLAB,
     %                          (MATRIX(J,K),K=FROM,TO)
            ELSE IF (STORED.EQ.'UPPDIA') THEN
              WRITE (*,FORM4) J,ATOMORB_LAB(J)%ATNUM,
     %                          ATOMORB_LAB(J)%ATLAB,
     %                          ATOMORB_LAB(J)%ORBLAB,
     %                          (MATRIX(J,K),K=FROM,TO)
            ELSE IF (STORED.EQ.'TOTAL') THEN
              WRITE (*,FORM4) J,ATOMORB_LAB(J)%ATNUM,
     %                          ATOMORB_LAB(J)%ATLAB,
     %                          ATOMORB_LAB(J)%ORBLAB,
     %                          (MATRIX(J,K),K=FROM,TO)
            ELSE
              GO TO 99920
            END IF
C
          END IF
   10   CONTINUE
        WRITE (*,'(/)')
   20 CONTINUE
C
      GO TO 99990
C
C     ------------------------------------------------------------------
C
C     *** Handling of output errors ***
C
99920 CONTINUE
      WRITE (*,99925) STORED
99925 FORMAT (/,T2,'*** ERROR IN PRILMAT! INADMISSIBLE ',
     $             'CLASSIFICATION OF MATRIX STORAGE ***',/,
     $        /,T2,'           STORED = ', A)
      GO TO 99980
C
99930 CONTINUE
      WRITE (*,99935) AFTER
99935 FORMAT (/,T2,'*** ERROR IN PRILMAT! INADMISSIBLE VALUE ',
     $             'FOR LETTERS AFTER DECIMAL POINT ***',/,
     $        /,T2,'            AFTER = ', I5)
      GO TO 99980
C
99940 CONTINUE
      WRITE (*,99945) BEFORE
99945 FORMAT (/,T2,'*** ERROR IN PRILMAT! INADMISSIBLE VALUE ',
     $             'FOR LETTERS BEFORE DECIMAL POINT ***',/,
     $        /, T2,'           BEFORE = ', I5)
      GO TO 99980
C
99950 CONTINUE
      WRITE (*,99955) N1,N2
99955 FORMAT (/,T2,'*** ERROR IN PRILMAT! INADMISSIBLE VALUE ',
     $             'FOR COLUMN NUMBER ***',/,
     $        /,T2,' First column: N1 = ', I5,
     $        /,T2,'  Last column: N2 = ', I5)
      GO TO 99980
C
99960 CONTINUE
      WRITE (*,99965) M1,M2
99965 FORMAT (/,T2,'*** ERROR IN PRILMAT! INADMISSIBLE VALUE ',
     $             'FOR ROW NUMBER ***',/,
     $        /,T2,'    First row: M1 = ', I5,
     $        /,T2,'     Last row: M2 = ', I5)
      GO TO 99980
C
99970 CONTINUE
      WRITE (*,99975)
99975 FORMAT (/,T2,'*** ERROR IN PRILMAT! INADMISSIBLE ',
     $             'CLASSIFICATION OF OUTPUT TAPE ***',/,
     $        /,T2,'             TAPE = ', I5)
      STOP
C
99980 CONTINUE
      WRITE (*,99985)
99985 FORMAT (/,/,T2,'****************************',
     $          /,T2,'*** Abnormal termination ***',
     $          /,T2,'****************************')
      STOP
C
C     *** Normal termination ***
C
99990 CONTINUE
C
C     ------------------------------------------------------------------
C
      END
