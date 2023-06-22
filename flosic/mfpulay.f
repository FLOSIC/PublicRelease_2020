C UTEP Electronic Structure Lab (2020)

      SUBROUTINE FPULAY
C ORIGINAL VERSION BY KOBLAR A JACKSON (1990)
C
C CALCULATE PULAY CORRECTIONS BY NUMERICAL INTEGRATION 
C
       use debug1
       use common2,only : RIDT, NIDENT, FHELLF, FNONL, FTOT, FRC1, FRC2,
     &                    IFUIDT, ZNUC, EDISP, over1, ek, dftV, allH
       use common3,only : RMAT
       use SICMAT,only  : DERSIC
       use FRM,only     : LFRM
       use SICFLAG,only : LSICF
       use global_inputs
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:28:41 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IDVERSION, IID, IPAR, IX
       INTEGER :: IWF, JWF, IORBX,NWFTOT !Variables used for Force 0.0
       REAL*8 :: EDISPD3 , FDIFF, FDISPD3, FPUL, SMALL, TIME1, TIME2

       dimension FDISPD3(3,MAX_IDENT)

       DIMENSION FPUL(3)
       DATA SMALL/1.0D-3/

       CALL CHECK_INPUTS
       EDISP=0
C
C ZERO FRC ARRAYS
C
c      DO IID=1,NIDENT 
c       DO IX=1,3
c        FRC1(IX,IID)=0.0D0
c        FRC2(IX,IID)=0.0D0
c       END DO
c      END DO
         FRC1 = 0.0d0
         FRC2 = 0.0d0
C Force 0.0
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       IF(LSICF) THEN
c
c  print LFRM to check occupations
c
c       print *, 'mfpulay lfrm', (lfrm(ix),ix=1,2)
c       print *, 'mfpulay nwfs', nwf, (nwfs(ix),ix=1,2)
        NWFTOT = LFRM(1) + LFRM(2)
c
c zero test arrays
c
c       do iwf = 1, NWFTOT
c        do jwf = iwf, NWFTOT
c         over1(iwf,jwf) = 0.0d0
c         over1(jwf,iwf) = 0.0d0
c         ek(iwf,jwf) = 0.0d0
c         ek(jwf,iwf) = 0.0d0
c         dftV(iwf,jwf)=0.0d0
c         dftV(jwf,iwf)=0.0d0
c         allH(iwf,jwf)=0.0d0
c         allH(jwf,iwf)=0.0d0
c        end do
c       end do
       ENDIF
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C CALCULATE PULAY CORRECTIONS NUMERICALLY
C
       CALL GTTIME(TIME1)
       CALL NUMFORCE
       CALL GTTIME(TIME2)
       CALL TIMOUT('PULAY CORRECTIONS TO FORCES:       ',TIME2-TIME1)
 1010  FORMAT(1X,A,7X,F12.3)
CForce 0.0
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
c Add SIC part to FRC1
c  At beginning of this subroutine
c    FRC1 = sum_i <del(phi_i)/del(R_nu)|HDFT|phi_i>
c  Adding: sum_i <del(phi_i)/del(R_nu)|VSIC_i|phi_i>
c  Notes 
c    i) the sum is over all occupied orbitals, spin up and down
c    ii)the - sign is needed because of the def of VSIC_i earlier
c  This should be fixed to make the sign more transparent here.
c
       IF(LSICF) THEN
        do IID =1,NIDENT
         do iorbx = 1,NWFTOT              !KAJ 1-24-2018
          do ix = 1,3
           FRC1(IX,IID) = FRC1(IX,IID) - DERSIC(IX,iorbx,IID)
          end do
         end do
        end do
       END IF
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C COMPUTE TOTAL FORCE ON IDENTITY MEMBERS
C SYMMETRIZE PULAY FORCE BEFORE ADDING TO HELLMANN-FEYNMAN TERM
C
       DO 20 IID=1,NIDENT
        DO IX=1,3
         FPUL(IX)= -2*(FRC1(IX,IID)-FRC2(IX,IID))
        END DO
        CALL FRCSYM(RIDT(1,IID),FPUL,FDIFF)
        IF (FDIFF .GT. SMALL) THEN
         PRINT 1020,IID,FDIFF 
 1020    FORMAT(' WARNING: PULAY FORCE CORRECTION OF ATOM ',I3,
     &          ' VIOLATES SYMMETRY BY ',D12.4)
        END IF
        DO IX=1,3
         FTOT(IX,IID)= FPUL(IX)
        END DO
   20  CONTINUE
C--------------------------------------------
C     Compute vanderWall forces using Grimme's DFTD3 method ! raja 
       IF(DFTD31) THEN
        write(6,*) '----------------------------------------'
        idversion = 3
        write(6,*) 'Calling  D_DFTd3 version 3'
        call dftd3_grimme(idversion,edispd3,FDISPD3)
         
        EDISP=EDISPD3
        write(6,*) '----------------------------------------'
       ENDIF
C--------------------------------------------
C
C OBTAIN TOTAL FORCE: HELLMANN-FEYNMAN PART ALREADY IN FHELLF
C
       PRINT '(A)',' '
       PRINT '(A)','HERE ARE THE FORCES:'
       PRINT '(A)','===================='
       IPAR=0
       DO IID=1,NIDENT
        DO IX=1,3
         IPAR=IPAR+1
         FTOT(IX,IID)= FTOT(IX,IID)+FHELLF(IX,IID)+FNONL(IX,IID) 
         if (DFTD31) then
           FTOT(IX,IID) =  FTOT(IX,IID) + FDISPD3(IX,IID)
         endif
        ENDDO

        PRINT 1030,'ATOM ',IID,', POSITION:',(RIDT(IX,IID),IX=1,3)
        PRINT 1040,'HELLMANN-FEYNMAN:      ',(FHELLF(IX,IID),IX=1,3)
        PRINT 1040,'NONLOCAL FORCE:        ',(FNONL(IX,IID),IX=1,3)
        PRINT 1040,'PULAY CORRECTION:      ',
     &             (FTOT(IX,IID)-FHELLF(IX,IID)-FNONL(IX,IID),IX=1,3)
        PRINT 1040,'PULAY CORRECTION1:     ',
     &             (FRC1(IX,IID),IX=1,3)
        PRINT 1040,'PULAY CORRECTION2:     ',
     &             (FRC2(IX,IID),IX=1,3)
        if (DFTD31) then
           PRINT 1040,'DISPERSION CORRECTION: ',
     &             (FDISPD3(IX,IID),IX=1,3)
        endif 
        PRINT 1040,'TOTAL:             ',(FTOT(IX,IID),IX=1,3)
       END DO 
       PRINT '(A)',' '
 11    format(3x,3(f13.6,1x))
 12    format(3x,4(f13.6,1x))
 13    FORMAT(I3,3(1X,F13.6,1x))
 1030  FORMAT(A,I3,A,3(1X,F18.12))
 1040  FORMAT(A,3(1X,D18.12))
cForce 0.0
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
c  print the test integrals:  nwftot = total of up and down orbitals
c   over1 = <phi_i|phi_j>
c   ek = <phi_i|-del^2/2|phi_j>
c   dftV = <phi_i|VDFT|phi_j>
c   allH = <phi_i|HDFT|phi_j>
c
!      IF(LSICF) THEN
!       print *, 'over1'
!       do iwf = 1,nwftot
!        print 83, (over1(jwf,iwf),jwf=1,nwftot)
!       end do
!       print *, 'ek'
!       do iwf = 1, nwftot
!        print 83, (ek(jwf,iwf),jwf=1,nwftot)
!       end do
!       print *, 'dftV'
!       do iwf = 1, nwftot
!        print 83, (dftV(jwf,iwf),jwf=1,nwftot)
!       end do
!       print *, 'allH'
!       do iwf = 1, nwftot
!        print 83, (allH(jwf,iwf),jwf=1,nwftot)
!       end do
 83      format(1x,10(1x,f10.4))
!       print *, 'return from mfpulay'
!      END IF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       RETURN
      END
