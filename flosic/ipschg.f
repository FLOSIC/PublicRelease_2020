C UTEP Electronic Structure Lab (2020)
          INTEGER  FUNCTION IPSCHG(NAMPSP,IZNUC)
C
C DEFINES THE VALENCE CHARGE FOR PSEUDOPOTENTIAL CALCULATIONS
C
       use common2,only : ZNUC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MXZBHS, IZNUC, NZVBHS, i,
     &            Ncore, LMax, NZVECP, MXZECP !<<<< ECP
       REAL*8 :: SYMBOL 
       SAVE
       PARAMETER (MXZBHS=94)
       PARAMETER (MXZECP=36)
       CHARACTER*3 NAMPSP
       DIMENSION NZVBHS(MXZBHS)
       DIMENSION NZVECP(MXZECP)
       DATA (NZVBHS(I),I=1,MXZBHS)/
     &       1,2,
     &       1,2,3,4,5,6,7,8,
     &       1,2,3,4,5,6,7,8,
     &       1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
     &       1,2,3,4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
     &       1,2,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
     &         4,5,6,7,8,9,10,11,12,3,4,5,6,7,8,
     &       1,2,3,4,5,6,7,8/

!<<<< ECP
       DATA (NZVECP(I),I=1,MXZECP)/
     &       1,2,
     &       1,2,3,4,5,6,7,8,
     &       1,2,3,4,5,6,7,8,
     &       9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26/
!>>>>

       IF (IZNUC .LE. 0) THEN
        IPSCHG=0
        RETURN
       END IF
       IF (NAMPSP .EQ. 'ALL') THEN
        IPSCHG=IZNUC
C<<<< rajendra 
       ELSE IF (NAMPSP .EQ. 'ECP') THEN
        IF (IZNUC .LE. MXZECP) THEN
         IPSCHG=NZVECP(IZNUC)
        END IF
C>>>> rajendra 
       ELSE IF (NAMPSP .EQ. 'BHS') THEN
        IF (IZNUC .LE. MXZBHS) THEN
         IPSCHG=NZVBHS(IZNUC)
        ELSE
         PRINT *,'IPSCHG: BHS PSP NOT DEFINED FOR Z= ',IZNUC
         CALL STOPIT
        END IF
       ELSE
        PRINT *,'IPSCHG: UNKNOWN PSP TYPE: ',NAMPSP
        CALL STOPIT
       END IF
       RETURN
      END
