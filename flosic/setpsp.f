C UTEP Electronic Structure Lab (2020)
C
C *******************************************************************
C
      SUBROUTINE SETPSP(IUNIT,SYMPSP,IZNUC)
C
C DIRK POREZAG, AUGUST 1997
C SETUP PSEUDOPOTENTIAL (PSP) DATA 
C
       use common2,only : ZNUC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:02 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IUNIT, IZNUC
       SAVE
       CHARACTER*7 SYMPSP
       CHARACTER*3 NAMPSP
C
       NAMPSP=SYMPSP(1:3)
       IF (NAMPSP .EQ. 'ALL') RETURN
       IF (NAMPSP .EQ. 'ECP') RETURN
       IF (NAMPSP .EQ. 'BHS') THEN
        CALL SBHSPSP(IUNIT,SYMPSP,IZNUC)
       ELSE
        PRINT *,'SETPSP: PSEUDOPOTENTIAL TYPE ',NAMPSP,' IS ',
     &          'NOT RECOGNIZED'
        CALL STOPIT
       END IF
       WRITE(IUNIT,'(A3)') '***'
       RETURN
      END
