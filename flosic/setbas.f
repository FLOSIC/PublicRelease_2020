C UTEP Electronic Structure Lab (2020)
C
C *******************************************************************
C
      SUBROUTINE SETBAS(IZNUC,NAMPSP,ALP,CON,NALP,NBASF)
C
C DIRK POREZAG, AUGUST 1997
C SETUP BASIS SET DEPENDING ON NUCLEAR CHARGE AND PSP
C
C       use common2,only : ZNUC
       use global_inputs
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:02 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IZNUC, NALP, NBASF
       REAL*8 :: ALP , CON
       SAVE
       CHARACTER*3  NAMPSP
       CHARACTER*10 BASIS_TXT
       DIMENSION ALP(MAX_BARE),CON(MAX_BARE,MAX_CON,3),NBASF(2,3)
C
       CALL CHECK_INPUTS
C DEFAULT BASIS SET IN NRLMOL
       IF(CALC_BASIS.EQ.1) THEN
         IF (NAMPSP .EQ. 'ALL') THEN
           CALL SALLBAS(IZNUC,ALP,CON,NALP,NBASF)
C   rajendra- use ALL electron basis for ECP
         ELSE IF (NAMPSP .EQ. 'ECP') THEN
         CALL SALLBAS(IZNUC,ALP,CON,NALP,NBASF)
C
         ELSE IF (NAMPSP .EQ. 'BHS') THEN
           CALL SBHSBAS(IZNUC,ALP,CON,NALP,NBASF)
         ELSE
           PRINT *,'SETBAS: PSEUDOPOTENTIAL TYPE ',NAMPSP,' IS ',
     &          'NOT RECOGNIZED'
           CALL STOPIT
         END IF
C RUN OTHER BASIS SET
       ELSE
         IF(NAMPSP .EQ. 'BHS')THEN
           PRINT *,'SETBAS: PSEUDOPOTENTIAL TYPE ',NAMPSP,' IS ',
     &          'NOT SUPPORTED FOR EXTERNAL BASIS SETS'
           CALL STOPIT
         !This ECP condition is needed.
         ELSEIF(NAMPSP .EQ. 'ALL'.OR.NAMPSP.EQ.'ECP')THEN
           CALL SALLBAS3(CALC_BASIS,IZNUC,ALP,CON,NALP,NBASF)
         ENDIF
       ENDIF
       RETURN
      END
