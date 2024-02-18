! UTEP Electronic Structure Lab (2020)
!###############################################################################
subroutine wfxmos
!This file will create a gto temp file compatible with MOLDEN2AIM
!Modiefied version of MOLDENMOS subroutine written by
!SUBHENDU PAUL (2011)
use common2,only : ISPN, NSPN, RIDT
! Conversion to implicit none.  Raja Zope Sun Aug 20 11:17:36 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
  INTEGER :: ICOUNT, IL, J, JAT, JATOMS, JJ, JSET, L, M_NUC, &
     & MSITES, NALP, NBASF, NCOUNT, NFIND, NSETS, NSUM, NUPP
  REAL*8 :: ALPSET , CHARGE_E, CHARGE_N, CONSET, REQV
!LOCAL VARIABLES:
COMMON/TMP1/REQV(3,MX_GRP)
CHARACTER*1 NSPD(3)
CHARACTER*132 INLINE
DIMENSION ALPSET(MAX_BARE),CONSET(MAX_BARE,MAX_CON,3)
DIMENSION NBASF(2,3)
REAL*8 FJL
DATA NSPD/'S','P','D'/

!------------------------------------------------------------------
!--- OPEN AUXILIARY FILES ---
OPEN(20,file='ISYMGEN')
!--- OPEN AND WRITE TMP FILE ---
OPEN(991,FILE='GTO.TMP',FORM='FORMATTED',STATUS='UNKNOWN')
REWIND(991)

100 CONTINUE
BACKSPACE(991)

!--- READ BASIS FUNCTIONS ---
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

!--- GET NUMBER OF EQUIVALENT ATOMS ---
  CALL GASITES(1,RIDT(1,JSET),M_NUC,REQV,MSITES)
  NUPP = MAX(NFIND,M_NUC)

  DO JAT = 1, NUPP
    ICOUNT = ICOUNT + 1
!---  NOW WRITE INTO MOLDEN FILE
    WRITE(991,*)NCOUNT+JAT,'0'
    DO L = 1, 3
      DO IL = 1, NBASF(1,L) !+NBASF(2,L)
        NSUM = 0
        DO JJ = 1, NALP
          IF(CONSET(JJ,IL,L).NE.0) THEN
            NSUM = NSUM + 1
          END IF
        END DO
        WRITE(991,'(1X,A1,3X,I2,1X,A3)')NSPD(L),NSUM,'1.0'
        DO J = 1, NALP
          IF(CONSET(J,IL,L).NE.0) THEN
!--- NORMALIZE BASIS SETS ---
! YY. You shouldn't normalize basis sets unless you diagonalize
! Hamiltonian again.
          !  FJL = 1.0/(3.14**(0.75))
          !  FJL = FJL*ALPSET(J)**(0.75+0.5*(L-1))
          !  FJL = FJL*(SQRT(2.0)**(L-1))
          !  FJL = 1.0d0
!            WRITE (991,'(2D18.10)')ALPSET(J),CONSET(J,IL,L)/FJL
            WRITE (991,'(2D20.10)')ALPSET(J),CONSET(J,IL,L)
          END IF
        END DO
      END DO
    END DO
    WRITE(991,*)'E 0'  !YY. Writing E to indicate the end of the section
  END DO
  NCOUNT = ICOUNT
END DO

READ(20,*)
READ(20,*)
!--- DONE, CLOSE UP SHOP AND CONTINUE ---
CLOSE(20)
CLOSE(991)
RETURN
end subroutine wfxmos
