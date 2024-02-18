C UTEP Electronic Structure Lab (2020)
      SUBROUTINE MOLDENMOS
C
C     WRITE BASIS SETS AND MO COEFFICIENTS INTO THE MOLDEN FILE AND
C     ALLOWS MOLDEN TO GRAPHICALLY DISPLAY RESULTS.
C
C     SUBHENDU PAUL (2011)
C     ULISES REVELES (JUNE 2013)
C
C     ------------------------------------------------------------------
C
C     COMMON VARIABLES:
C
      use common2,only : ISPN, NSPN, RIDT
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: ICOUNT, IL, ITEST, J, JAT, JATOMS, JJ, JSET, L,
     & M_NUC, MSITES, NALP, NBASF, NCOUNT, NFIND, NSETS, NSUM, NUPP
       REAL*8 :: ALPSET , CHARGE_E, CHARGE_N, CONSET, REQV
C
C     LOCAL VARIABLES:
C
      COMMON/TMP1/REQV(3,MX_GRP)
C
      CHARACTER*1 NSPD(3)
      CHARACTER*132 INLINE
C
      DIMENSION ALPSET(MAX_BARE),CONSET(MAX_BARE,MAX_CON,3)
      DIMENSION NBASF(2,3)
C
      REAL*8 FJL
C
      DATA NSPD/'s','p','d'/
C
C     ------------------------------------------------------------------
C
C     --- OPEN AUXILIARY FILES ---
C
      OPEN(20,file='ISYMGEN')
C
C     --- OPEN AND APPEND MOLDEN FILE ---
C
      OPEN(99,FILE='CLUSTER.MOLDEN',FORM='FORMATTED',STATUS='OLD')
C
C     --- SCAN MOLDEN FILE FOR THE BASIS SET HEADER, APPEND OR ---
C     --- WRITE OVER THE NEW DATA AS NEEDED ---
C
      ITEST = 0
      DO WHILE(ITEST.EQ.0)
        READ(99,'(132a)',END=100) INLINE
        DO J=1,20
          IF(INLINE(J:J+3).EQ.'[6D]') THEN
            ITEST=1
            EXIT
          END IF
        END DO
      END DO
C
  100 CONTINUE
      BACKSPACE(99)
C
      WRITE(99,*)'[6D]'
      WRITE(99,*)'[10F]'
      WRITE(99,*)'[GTO]'
C
C     --- READ BASIS FUNCTIONS ---
C
      ICOUNT = 0
      NCOUNT = 0
      READ(20,*)NSETS
      DO JSET = 1, NSETS
        READ(20,*)charge_e,charge_n
        READ(20,*)
        READ(20,*)NFIND
        DO JATOMS = 1,NFIND
          READ(20,*) 
        END DO
        READ(20,*)
        READ(20,*)NALP
        READ(20,*)(NBASF(1,L),L=1,3)
        READ(20,*)(NBASF(2,L),L=1,3)
        READ(20,*)(ALPSET(J),J=1,NALP)
        READ(20,*)
        DO L = 1, 3
          DO IL = 1, NBASF(1,L)+NBASF(2,L)
            READ(20,*)(CONSET(J,IL,L),J=1,NALP)
            READ(20,*)
          END DO
        END DO
C
C     --- GET NUMBER OF EQUIVALENT ATOMS ---
C
        CALL GASITES(1,RIDT(1,JSET),M_NUC,REQV,MSITES)
        NUPP = MAX(NFIND,M_NUC)
C
        DO JAT = 1, NUPP
          ICOUNT = ICOUNT + 1
C
C     ---  NOW WRITE INTO MOLDEN FILE
C
          WRITE(99,*)NCOUNT+JAT,'0'
          DO L = 1, 3
            DO IL = 1, NBASF(1,L) !+NBASF(2,L)
C
              NSUM = 0
              DO JJ = 1, NALP
                IF(CONSET(JJ,IL,L).NE.0) THEN
                  NSUM = NSUM + 1
                END IF
              END DO 
              WRITE(99,'(1X,A1,3X,I2,2X,A4)')NSPD(L),NSUM,'1.00'
C
              DO J = 1, NALP 
                IF(CONSET(J,IL,L).NE.0) THEN
C
C     --- NORMALIZE BASIS SETS ---
C
                  FJL = 1.0/(3.14**(0.75))
                  FJL = FJL*ALPSET(J)**(0.75+0.5*(L-1))
                  FJL = FJL*(SQRT(2.0)**(L-1)) 
C
                  WRITE (99,'(2D18.10)')ALPSET(J),CONSET(J,IL,L)/FJL
                END IF
              END DO
            END DO
          END DO
          WRITE(99,*)' '
        END DO
        NCOUNT = ICOUNT
      END DO
C
      READ(20,*)
      READ(20,*)
C
C     --- DONE, CLOSE UP SHOP AND CONTINUE ---
C
      CLOSE(20)
      CLOSE(99)
C
C     --- NOW WRITE MO's ENERGIES AND OCCUPATIONS ---
C
      CALL WRITEUS(0) 
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END
