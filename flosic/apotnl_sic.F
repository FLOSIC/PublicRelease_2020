C UTEP Electronic Structure Lab (2020)
C> @file apotnl_sic.ftn

C *****************************************************
C>  
C> @author APOTNL (DIRK POREZAG, 1994, 1998) (MARK PEDERSON 1987-1994)
C> @brief APOTNL CALCULATES THE POTENTIAL ON THE MESHPOINTS WITHIN LSDA OR GGA
C> @note ATTENTION: FIRST TWO ARRAYS OF COMMON BLOCK TMP1 MUST BE IDENTICAL IN
C> DENSOLD AND APOTNL SINCE THEY ARE USED TO PASS DENSITY AND COULOMB POT
C> 
C> @details Detail description of apotnl
C> @return  returns some values LaTeX test \f$ \sum_i^n x_5  \f$
C> @param[out] TOTQNUM Total charge in apotnl
C> @param[in] LSPX Spin index (=0 for DFT)
C> @param[in] IORBX FO index (=0 for DFT)
C> @param[out] RESULTS contains SIC energies and potentials for every FO.
C> @param[out] RESULTS(1) JRBSIC:FO Index
C> @param[out] RESULTS(2) XCH1: Charge
C> @param[out] RESULTS(3) XCH2: Charge (for RHOG down. Should be 0).
C> @param[out] RESULTS(4) ASY: The outermost mesh distance times coulomb
C>                             potential at that point.
C> @param[out] RESULTS(5) XEXX: LDA Exchange energy (The value is not used)
C> @param[out] RESULTS(6) -XCOU/2: Coulomb energy
C> @param[out] RESULTS(7) XDEL: Integral of rhog*(-Vxc - Vcoul) Or Vsic
C> @param[out] RESULTS(8) XENG: Exc-sic - Ecoul/2
C> @param[out] RESULTS(9) SIC_ENG: SIC energy (XENG *2/N)
C> @param[out] RESULTS(10) ERGXL: Local exchange
C> @param[out] RESULTS(11) ERGCL: Local correlation
C> @param[out] RESULTS(12) ERGXN: Nonlocal exchange
C> @param[out] RESULTS(13) ERGCN: Nonlocal correlation
       SUBROUTINE APOTNL_SIC(TOTQNUM,LSPX,IORBX)
       use debug1
       use mpidat1,only  : IRANK,MYGROUP
       use global_inputs,only : SOLVENT1,EXCITED1,SPNORB1,PCM1,ITTOT,
     &                          LIBXC1,CALCTYPE1
       use pot_dens,only : COULOMB,RHOG
       use mixpot1,only  : POTIN,POT=>POTOUT
       use mesh1,only    : WMSH,RMSH,NMSH
       use common2,only  : RCNT, IFUCNT, NIDENT, NCNT, ZELC, IGGA,
     &   NSPN, DIPOLE, EFIELD, ELOCAL, ECOUL,  ERGFLD,ESOLC,
     &   ERGXL,ERGXN,ERGCL,ERGCN
       use common3,only  : RMAT, NGRP
       use common5,only  : CONVERGENCE, NWFS
       use common7,only  : MODDEN, GAUSS_CUT
       use solvent,only  : POTSOL
!SIC modules
       use SICFLAG,only : LSICF,MESH_SAV
       use FRM,only     : BFRM,RESULTS,LFRM
       use HMATSIC,only : OVTM,HMTM
       use FOCENT,only  : CFRM
       use SICMAT,only  : SIC,DERSIC,ZPOT
!      use DIRECTIONS !NSPN_SKIP
       use LOCORB,only : IRBSIC
       use scaledpzsic,only : AVGSICON,AVGSIC,scaledsic
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:46 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LSPX, IORBX, I, IAT, IMSH, IPTS, ITER, IWF, IERR,
     & J, JCALC, JMS2, JMS3, JRBSIC, JWF, KNT, LSPN, LSXX, MMSH, MXXD,
     & NSPN_SAV, IX, IID
       REAL(8) :: TOTQNUM , ADDIT, APT1, ASY, AVG, COU_TIME_TOT, DEL2,
     & DEL3, DERIV, DINT, EXCSIC, FRMTIM, PI, R2, R4, R6, RAD, RM0,
     & RM2, RM4, RM6, RMX, SIC_COU, SIC_ENG, SIC_EXO, SUM1, TIME0,
     & TIME1, TIME2, TIME3, TIME4, TIME_GETFRMFRC, TIME_GETSICMAT,
     & TMEL1, TMEL2, XALPHA, XCGV, XCH1, XCH2, XCOU, XDEL, XENG, XEXC,
     & XEXX, YCGV, ZCGV,ZCOUL,DINT2,DINT3,ZEXCH,ZSUM,FACT
!      INCLUDE 'commons.inc'
       SAVE
C
C FOOL THE COMPILER FOR MXSPN=1 TO SUPRESS WARNING MESSAGES
C THAT ARE REALLY IRRELEVANT
C
C       COMMON/LOCORB/TMAT(MAX_OCC,MAX_OCC,MXSPN),MORB(2),IRBSIC,ZSIC
!       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
!       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
!       COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!     &   LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)
!       COMMON/HMATSIC/ OVTM(MAX_OCC,MAX_OCC,2),HMTM(MAX_OCC,MAX_OCC,2)
!       COMMON/FOCENT/ CFRM(3,MAX_OCC,MXSPN)
!       COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)
!       COMMON/DIRECTIONS/NSPN_SKIP
       LOGICAL FIRST,AVERAGE,EXIST,FIRST_MSH
       CHARACTER(7)  NAMES
       CHARACTER(30) FILENAME
       CHARACTER(3)  ORBTXT
       CHARACTER(3)  ITTOTTXT
       CHARACTER(12) ZPOTSTR
       DIMENSION NAMES(3),DERIV(3)
C
       DATA FIRST/.TRUE./
       DATA FIRST_MSH/.TRUE./
       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/
       DATA AVERAGE/.FALSE./
       DATA COU_TIME_TOT/0.0D0/
C      DATA XALPHA/0.8225/! (Kohn Sham=1, 0.8225=Best Ar for PW91)
       DATA XALPHA/1.0000/! (Kohn Sham=1, 0.8225=Best Ar for PW91)
       DATA TIME_GETSICMAT,TIME_GETFRMFRC/2*0.0D0/
C
C      OPEN(21,FILE='XALPHA')
C      REWIND(21)
C      READ(21,*)XALPHA
C      CLOSE(21)
C
C DIRECTIONS FOR CALCULATING DENSITY:
C MODDEN=1 USE DENSITY MATRIX
C MODDEN=2 USE WAVEFUNCTIONS
C
       MODDEN=1
       WRITE(6+IRANK,'(A)')
     &       'CALCULATING FERMI ORBITAL POTENTIAL ON MESH POINTS'
       WRITE(6+IRANK,*)'CURRENTLY ORBITAL,SPIN',IORBX,LSPX
       CALL GTTIME(APT1)
       ITER=0
C
C READ AVERAGE
C
       AVERAGE=.TRUE.
!
! Cannot write files in parallel, so skip this in group mode
! Global manager will take care of this in the encapsulating global call
!
#ifndef GROUP
       OPEN(99,FILE='AVRGDAT',FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND 99
       READ(99,*,END=30)AVG,AVERAGE
   30  CONTINUE
       REWIND(99)
       WRITE(99,*)AVG,AVERAGE,LSICF,' AVG, AVERAGE,LSICF'
       CLOSE(99)
#endif
!      CALL TRACER('AFTER AVRGDAT')
!      IF(LSICF)THEN
!        AVG=-AVG
!      END IF
       IF(FIRST)THEN
         FIRST=.FALSE.
C
C SET GAUSS_CUT (IN THIS VERSION ONLY NEEDED BY COULOMB1)
C
         DO IAT=1,NIDENT
           GAUSS_CUT(IAT)=1.0D30
         END DO
       END IF
C
C START UP POTENTIAL CALCULATION: CHECK ARRAY DIMENSIONS
C
       IF(NMSH.GT.MAX_PTS) THEN
         PRINT *,'APOTNL: MAX_PTS MUST BE AT LEAST', NMSH
         CALL STOPIT
       END IF
C
C CALCULATE ELECTRONIC COULOMB POTENTIAL AND DENSITY
C
C      IF(IORBX.EQ.1)THEN
C      OPEN(17,FILE='TMATSND',FORM='FORMATTED')
C      DO K=1,2
C      WRITE(17,*)'K:',K
C      WRITE(17,17)((TMAT(J,I,K),J=1,24),I=1,24)
C17    FORMAT(' ',10F12.3)
C      END DO
C      CLOSE(17)
C      END IF
#ifdef MPI
#ifdef GROUP
!      CALL TRACER('STARTING SENNDATA')
       CALL SENDDATA_GRP(101)
       CALL SENDDATA_GRP(102)
       CALL SENDDATA_GRP(103)
       CALL SENDDATA_GRP(104)
       CALL SENDDATA_GRP(105)
C      CALL SENDDATA_GRP(106)
       CALL SENDDATA_GRP(208)
       CALL SENDDATA_GRP(209)
!      CALL SENDDATA_GRP(211)
!      CALL TRACER('DONE SENNDATA')
#else
       CALL SENDDATA(101)
       CALL SENDDATA(102)
       CALL SENDDATA(103)
       CALL SENDDATA(104)
       CALL SENDDATA(105)
C      CALL SENDDATA(106)
       CALL SENDDATA(208)
       CALL SENDDATA(209)
!      CALL SENDDATA(211)
#endif
#endif
       IF(IORBX.EQ.0) then
        print *,"APOTNL_SIC is trying to run DFT calculation"
        call STOPIT
       ENDIF
!##########################################
!##          SIC Section                 ##
!##########################################
C      CALL GTTIME(TIME0)
C      CALL FRMORB2(0,0)
C      CALL GTTIME(TIME1)
C      WRITE(INITOUT,*)'TIME IN FRMORB2:',TIME1-TIME0
C      WRITE(6,*)'IRANK:CALLING FRMIORB',IRANK,IORBX
C      CALL GTTIME(TIME0)
C      CALL FRMIORB(1,1)
C      CALL GTTIME(TIME1)
C      WRITE(INITOUT,*)'TIME IN FRMIORB:',TIME1-TIME0
C      PRINT*,'BACK FROM FRMIORB'
C      CALL STOPIT
C      call flush(6)
C      CALL GTTIME(TIME1)

       IORBX=ABS(IORBX)
       IF(IORBX.EQ.1.AND.LSPX.EQ.1)THEN
         SIC_ENG=0.0D0
         SIC_COU=0.0D0
         SIC_EXO=0.0D0
C        OPEN(41,FILE='SICPOT',FORM='UNFORMATTED')
C        REWIND(41)
       END IF
       KNT=0
C LB taken off due to parallelization
C      IF(LSPX.EQ.1.and.IORBX.EQ.1)THEN
C        OPEN(76,FILE='SIC_ENERGY')
C        REWIND(76)
C      END IF
C LB
       FRMTIM=0.0D0
       LSPN=LSPX
       IRBSIC=IORBX+(LSPN-1)*NWFS(1)
       JRBSIC=IORBX
       KNT=IRBSIC
       CALL GTTIME(TMEL1)
C ORBITAL LOOP
       DO IMSH=1,NMSH
         COULOMB(IMSH)=0.0D0
         RHOG(IMSH,1:10,1)=0.0D0 ! 1-> 1:10
       END DO
       CALL GTTIME(TIME1)
       SUM1=0.0d0
       MXXD=0
       IF(MESH_SAV)THEN
         IF(ITTOT>=4) THEN
           CALL TRACER('READING REDUCED MESH',IORBX)
           WRITE(6+IRANK,*)'READING FILE:',FILENAME,IORBX
           INQUIRE(FILE=FILENAME,EXIST=EXIST)
           IF(.NOT.EXIST) THEN
             CALL TRACER('MESH FILE DOES NOT EXIST',IORBX)
             CALL STOPIT
           ENDIF
           OPEN(99,FILE=FILENAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
           READ(99)MMSH
           READ(99)((RMSH(J,I),J=1,3),I=1,MMSH)
           READ(99)(WMSH(I),I=1,MMSH)
           CLOSE(99)
           NMSH=MMSH
         ENDIF
       ENDIF
#ifdef GROUP
       CALL SENDDATA_GRP(102)
#else
       CALL SENDDATA(102)
#endif
       CALL COUPOT_SIC(MXXD,LSPN,IRBSIC)
       CALL GTTIME(TIME2)
       COU_TIME_TOT=COU_TIME_TOT+TIME2-TIME1
       WRITE(6+IRANK,*)'TIME FOR COUPOT2:',TIME2-TIME1,COU_TIME_TOT,
     &     COU_TIME_TOT/TIME2,LSPN,NMSH,IRBSIC
       CALL FLUSH(6+IRANK)

!Scaling
!Read RSIC factor here
       IF(scaledsic) call readrsicfactor(IORBX,LSPX)

       JMS2=0
       JMS3=0
       DEL2=0.0D0
       DEL3=0.0D0
       ADDIT=0.0D0
       DO IMSH=1,NMSH
         ADDIT=ADDIT+ABS(COULOMB(IMSH))
       END DO
       IF(ADDIT.LT.0.001)THEN
         PRINT*,'COUPOT IS ZERO'
         CALL STOPIT
       END IF
! LB DETERMINE IF WE ARE SAVING THE MESH FOR THIS ORIBTAL
       IF(MESH_SAV)THEN
         WRITE(ORBTXT,'(I0.3)')IORBX
         FILENAME='SCRATCH/VMOLD'//TRIM(ORBTXT)
         IF(ITTOT<=3)THEN
! CALCULATE IT ONCE AND DUMP IT TO FILE
           CALL TRACER('CALCULATING REDUCED MESH',IORBX)
           MMSH=0
           DO IMSH=1,NMSH
             DINT=WMSH(IMSH)*COULOMB(IMSH)*RHOG(IMSH,1,1)
             IF(ABS(DINT).GE.1.0D-11) THEN
               MMSH=MMSH+1
               DO J=1,3
                 RMSH(J,MMSH)=RMSH(J,IMSH)
               END DO
               WMSH(  MMSH)=WMSH(  IMSH)
               COULOMB(MMSH)=COULOMB(IMSH)
               RHOG(MMSH,1,1)=RHOG(IMSH,1,1)
             END IF
           END DO
           OPEN(99,FILE=FILENAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
           WRITE(6+IRANK,*)'WRITING FILE:',FILENAME,IORBX
           WRITE(99)MMSH
           WRITE(99)((RMSH(J,I),J=1,3),I=1,MMSH)
           WRITE(99)(WMSH(I),I=1,MMSH)
           
           CLOSE(99)
!         ELSE
! READ MESH FROM ORBITAL FILE
!           CALL TRACER('READING REDUCED MESH',IORBX)
!           WRITE(6+IRANK,*)'READING FILE:',FILENAME,IORBX
!           INQUIRE(FILE=FILENAME,EXIST=EXIST)
!           IF(.NOT.EXIST) THEN
!             CALL TRACER('MESH FILE DOES NOT EXIST',IORBX)
!             CALL STOPIT
!           ENDIF
!           OPEN(99,FILE=FILENAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
!           READ(99)MMSH
!           READ(99)((RMSH(J,I),J=1,3),I=1,MMSH)
!           READ(99)(WMSH(I),I=1,MMSH)
!           CLOSE(99)
         ENDIF
!      ELSE
! JUST CALCULATE YOUR REDUCED MESH FOR THIS ORBITAL
!        WRITE(ORBTXT,'(I0.3)')IORBX
!        WRITE(ITTOTTXT,'(I0.3)')ITTOT
!        MMSH=0
!        ZCOUL=0.0d0
!        ZEXCH=0.0d0
!        ZSUM=0.0d0
!        ZEXCH =0.0d0        
!        PI=4.0D0*ATAN(1.0D0)
!        FACT=0.75*((6/PI)**(1./3.))
!Mesh screening based on Coulomb energy.                            
!This screening does not work with GGA,meta-GGA.
!        DO IMSH=1,NMSH
!          DINT=WMSH(IMSH)*COULOMB(IMSH)*RHOG(IMSH,1,1)
!!         DINT2=WMSH(IMSH)*RHOG(IMSH,1,1)
!!         DINT3=WMSH(IMSH)*ABS(RHOG(IMSH,1,1))**(4./3.)
!          IF(ABS(DINT).GE.1.0D-11) THEN
!            MMSH=MMSH+1
!            DO J=1,3
!              RMSH(J,MMSH)=RMSH(J,IMSH)
!            END DO
!            WMSH(  MMSH)=WMSH(  IMSH)
!            COULOMB(MMSH)=COULOMB(IMSH)
!            RHOG(MMSH,:,1)=RHOG(IMSH,:,1)
!            JMS2=JMS2+1
!            DEL2=DEL2+DINT
!          END IF
!          DEL3=DEL3+DINT
!          ZCOUL=ZCOUL+DINT
!          ZSUM=ZSUM+DINT2
!          ZEXCH=ZEXCH+DINT3
!        END DO
!        ZEXCH=ZEXCH*FACT
!        ZCOUL=ZCOUL*0.50
       ENDIF
! LB END OF REDUCED MESH
!YY Rebroadcast mesh. It is necessary since our getvlxc is 
!   parallel whereas that of Perfect ver. is a serial.
!#ifdef MPI
!      CALL SENDDATA(102)
!#endif
!      WRITE(6+IRANK,2311)DEL2,DEL3,JMS2,NMSH,
!    &           FLOAT(JMS2)/NMSH
!      WRITE(6+IRANK,*)'INDEX:',IORBX,'ZCOUL:',ZCOUL,
!    &                  'ZSUM:',ZSUM,'ZEXCH:',ZEXCH
!      WRITE(6+IRANK,2312) IORBX,ZCOUL,ZSUM,ZEXCH
!2312  FORMAT('INDEX:',I0.2,' ZCOUL:',F10.6,' ZSUM:',F10.6,
!     &         ' ZEXCH:',F10.6)
 2311  FORMAT(2F12.4,2I10,F12.4,' SPARSITY SAVINGS')
!      NMSH=MMSH  !<<< Needed if mesh screening is done
!      WRITE(6+IRANK,*)'READUCED MESH',IORBX,NMSH
       XCGV=0.0D0
       YCGV=0.0D0
       ZCGV=0.0D0
       DO IMSH=1,NMSH
         XCGV=XCGV+RMSH(1,IMSH)*RHOG(IMSH,1,1)*WMSH(IMSH)
         YCGV=YCGV+RMSH(2,IMSH)*RHOG(IMSH,1,1)*WMSH(IMSH)
         ZCGV=ZCGV+RMSH(3,IMSH)*RHOG(IMSH,1,1)*WMSH(IMSH)
       END DO
       CFRM(1,JRBSIC,LSPN)=XCGV
       CFRM(2,JRBSIC,LSPN)=YCGV
       CFRM(3,JRBSIC,LSPN)=ZCGV
       RM0=0.0D0
       RM2=0.0D0
       RM4=0.0D0
       RM6=0.0D0
       DO IMSH=1,NMSH
         R2=(RMSH(1,IMSH)-XCGV)**2
     &     +(RMSH(2,IMSH)-YCGV)**2
     &     +(RMSH(3,IMSH)-ZCGV)**2
         R4=R2*R2
         R6=R4*R2
         RM0=RM0+RHOG(IMSH,1,1)*WMSH(IMSH)
         RM2=RM2+RHOG(IMSH,1,1)*WMSH(IMSH)*R2
         RM4=RM4+RHOG(IMSH,1,1)*WMSH(IMSH)*R4
         RM6=RM6+RHOG(IMSH,1,1)*WMSH(IMSH)*R6
       END DO
       RM2=RM2**(1./2.)
       RM4=RM4**(1./4.)
       RM6=RM6**(1./6.)
       WRITE(6+IRANK,2209)IRBSIC,XCGV,YCGV,ZCGV,RM0,RM2,RM4,RM6
 2209  FORMAT(' FLO-MOMENTS:',I5,9F12.4)
       RMX=0.0D0
       DO IMSH=1,NMSH
         RAD=SQRT(RMSH(1,IMSH)**2
     &      +RMSH(2,IMSH)**2
     &      +RMSH(3,IMSH)**2)
         COULOMB(IMSH)=COULOMB(IMSH)*FLOAT(NSPN)/2.
         IF(RAD.GT.RMX)THEN
           RMX=RAD
           ASY=COULOMB(IMSH)*RMX
         END IF
         RHOG(IMSH,:,2)=0.0D0 !!TEMPORARY FIX   1->1:10
       END DO
       NSPN_SAV=NSPN
       NSPN=2  !This needs to be broadcasted for getvlxc parallelization.
       DO IPTS=1,NMSH
         POT(IPTS)=0.0D0
         POT(IPTS+NMSH)=0.0D0
         POTIN(IPTS)=0.0D0
         POTIN(IPTS+NMSH)=0.0D0  ! CJK added line
       END DO

       ERGXL=0.0D0
       ERGXN=0.0D0
       ERGCL=0.0D0
       ERGCN=0.0D0

!Broadcast NSPN=2: No need to rebroadcast afterwards since SICLAGM
!doesn't use NSPN. Potentially move this to getvlxc.
!#ifdef MPI
!#ifdef GROUP
!       CALL SENDDATA_GRP(101)
!#else
!       CALL SENDDATA(101)
!#endif
!#endif

       CALL GETVLXC(NMSH,RHOG,POT,POTIN,.false.,MXXD,LSPN,IRBSIC)
!      CALL GETVLXC(MAX_PTS,RHOG,POT,POTIN)
CMP    CALL GETVLXC(NMSH   ,RHOG,POT,POTIN)
!      ERGXL=ERGXL*XALPHA
       PI=4.0D0*ATAN(1.0D0)
!YY XALPHA is set to 1 in parameter
!      DO IPTS=1,NMSH
!        POT(IPTS)=POT(IPTS)+(1.-XALPHA)*
!     &  ((6/PI)**(1./3.))*ABS(RHOG(IPTS,1,1))**(1./3.)
!      END DO
 
       NSPN=NSPN_SAV
C CALCULATE LOCAL CHARGE:
       XCH1=0.0D0
       XCH2=0.0D0
       XCOU=0.0D0
       XEXX=0.0D0
       XDEL=0.0D0
       XENG=0.0D0
C      PRINT*,'SIC WRITE:',LSPN,JRBSIC,KNT
C      WRITE(41)(POT(IMSH)+COULOMB(IMSH),IMSH=1,NMSH)

! Scaling PZSIC potential
       if(scaledsic) call scalepotential(IORBX,LSPX)

! Adding coulomb potential to POT.
! Later, subtract coulomb potential from SIC potential (= XC potential)
! for calculating RESULT array.
       DO IMSH=1,NMSH
         POT(IMSH)=POT(IMSH)+COULOMB(IMSH)
       END DO
       CALL GTTIME(TIME0)
       IF(JRBSIC.NE.IORBX)THEN
         PRINT*,'JRBSIC DIFFERENT THAN IORBX',JRBSIC,IORBX
         CALL STOPIT
       END IF
       CALL GTTIME(TIME2)
!      WRITE(6+IRANK,*)' CALLING SICLAGM', JRBSIC, LSPN

! Save POT for building SIC matrix later
! Draft 1
!      OPEN(888,FILE='SICP.dat',ACCESS='APPEND')
!      DO IMSH=1,NMSH
!       WRITE(888,*) POT(IMSH)
!      END DO
!      CLOSE(888)

! Draft 2
!      ZPOT(:,IORBX,LSPX)=POT(1:NMSH)

! Draft 3 - write POT NPBLOCK at a time
! File name string
       WRITE(ZPOTSTR,'(A,I4.4)')'ZPOT',IRBSIC
!Clean up
       INQUIRE(FILE=ZPOTSTR,EXIST=EXIST)
       IF(EXIST) THEN
        OPEN(206+IRBSIC,FILE=ZPOTSTR,FORM='UNFORMATTED',STATUS='OLD')
        CLOSE(206+IRBSIC,STATUS='DELETE')
       END IF
!Average SIC potential
       if(AVGSICON) then
        DO IMSH=1,NMSH
         !LSIC components
         IF(LSPN.EQ.1) THEN
          AVGSIC(IMSH)= AVGSIC(IMSH)+RHOG(IMSH,1,1)*POT(IMSH)
         ELSEIF(LSPN.EQ.2) THEN
          AVGSIC(IMSH+NMSH)= AVGSIC(IMSH+NMSH)+RHOG(IMSH,1,1)*POT(IMSH)
         END IF
        END DO
       else !Not AVGSIC potential
!Write file
        OPEN(206+IRBSIC,FILE=ZPOTSTR,FORM='UNFORMATTED',STATUS='NEW')
        print *,"ZPOT  created",IRBSIC
        DO IMSH=1,NMSH
         WRITE(206+IRBSIC) POT(IMSH)
        END DO
        CLOSE(206+IRBSIC)
       endif

!      CALL SICLAGM(LSPN,JRBSIC) !>>>> Moving outside of apotnl_sic
       CALL GTTIME(TIME3)

!      WRITE(6+IRANK,*)"SIC mat"
!      DO IWF=1,10
!       WRITE(6+IRANK,'(20G11.3)')(SIC(JWF,IWF,1),JWF=1,20)
!      END DO

!YY. FOD force is moved outside of SIC loop.
!      IF(CONVERGENCE)THEN
!<<<<<<<<<<<<<<<
!Force 0.0
!       !Note: coupot_sic has check_input call that sets the correct 
!       !calctype1 variable in group calculation.
!       IF(CALCTYPE1.NE.2) THEN
c        WRITE(6+IRANK,*) 'calling siclag_der jrbsic', jrbsic
c
c   Initialize DERSIC
c
!        CALL SICLAG_DER(LSPN,JRBSIC)
!        WRITE(6+IRANK,*) 'back from siclag_der jrbsic', jrbsic
!        iwf = JRBSIC + LFRM(1)*(LSPN - 1)   !iwf is the wave function index
!        WRITE(6+IRANK,*) 'iwf', iwf
!        do iid = 1,10   !iid is the atom index
!         WRITE(6+IRANK,'(3G11.3)')(DERSIC(IX,IWF,IID),ix=1,3)
!        end do
!       END IF
!
!>>>>>>>>>>>>>>
! FRMORB2 call uses the SIC array from the previous iteration.
! This is okay as long as orbitals are converged.
!        LSXX=-LSPN
!        WRITE(6+IRANK,*)'CALLING FRMORB2'
!        CALL FRMORB2( LSXX,JRBSIC)
!        CALL GTTIME(TIME4)
!        TIME_GETFRMFRC=TIME_GETFRMFRC+TIME4-TIME3
!        WRITE(6+IRANK,*)'TIME IN GETFRMFRC:',TIME_GETFRMFRC
!      END IF
!      TIME_GETSICMAT=TIME_GETSICMAT+TIME3-TIME2
!      WRITE(6+IRANK,*)'TIME IN GETSICMAT:',TIME_GETSICMAT
       DO IMSH=1,NMSH
         POT(IMSH)=POT(IMSH)-COULOMB(IMSH)
       END DO
       DO IMSH=1,NMSH
         XCH1=XCH1+WMSH(IMSH)*RHOG(IMSH,1,1)
         XCH2=XCH2+WMSH(IMSH)*RHOG(IMSH,1,2)
         XCOU=XCOU+WMSH(IMSH)*RHOG(IMSH,1,1)*COULOMB(IMSH)
         XEXX=XEXX+WMSH(IMSH)*ABS(RHOG(IMSH,1,1))**(4./3.)
         XDEL=XDEL-WMSH(IMSH)*RHOG(IMSH,1,1)*POT  (IMSH)
       END DO
       EXCSIC=ERGXL+ERGXN+ERGCL+ERGCN
       XEXC=XDEL
       XDEL=XDEL-XCOU
       XENG=-EXCSIC-XCOU/2.0d0
!      WRITE(6+IRANK,*)'PRELIM',IORBX,XDEL,XCOU,EXCSIC,XENG
       PI=4.0D0*ATAN(1.0D0)
       XEXX=XEXX*0.75*((6/PI)**(1./3.))
       SIC_ENG=XENG*(2.0d0/NSPN)
C      WRITE(INITOUT,*)'EXCSIC XCOU',EXCSIC,XCOU
       SIC_COU=SIC_COU-(2/NSPN)*XCOU/2.
       SIC_EXO=SIC_EXO+(2/NSPN)*XEXX
       CALL GTTIME(TMEL2)
C      WRITE(INITOUT,*)'FULL TIME:',TMEL2-TIME0

C LB taken off due to parallelization
C      WRITE(INITOUT,800)
C     &JRBSIC,XCH1,XCH2,ASY,XEXX,-XCOU/2.,XDEL,XENG,SIC_ENG,
C     &ERGXL,ERGCL,ERGXN,ERGCN
C              WRITE(76,800)
C     &JRBSIC,XCH1,XCH2,ASY,XEXX,-XCOU/2.,XDEL,XENG,SIC_ENG,
C     &ERGXL,ERGCL,ERGXN,ERGCN
       RESULTS(1, IORBX,LSPN) = JRBSIC
!> JRBSIC: Orbital index number. 1,2,...,#FODs.
       RESULTS(2, IORBX,LSPN) = XCH1
!> XCH1: Charge for spin up. This should be 1.
       RESULTS(3, IORBX,LSPN) = XCH2
!> XCH2: Charge for spin down. This should be 0.
       RESULTS(4, IORBX,LSPN) = ASY
       RESULTS(5, IORBX,LSPN) = XEXX
!> XEXX = (3/4)*(6/pi)^(1/3) * sum(WMSH*|RHOG|^(4/3))
!>      -> LDA Energy
       RESULTS(6, IORBX,LSPN) = -XCOU/2.0d0
!> -XCOU/2 = sum(-WMSH*RHOG*COULOMB)/2 !E_coul 
!>         -> Coulomb energy (nspn=2) or Coulomb energy/2 (nspn=1).
       RESULTS(7, IORBX,LSPN) = XDEL
!> XDEL = sum(-WMSH*RHOG*(POT+COULOMB)).
!> where POT=POT-COULOMB=XCPOT, so it's sum(-WMSH*RHOG*SICPOT).
       RESULTS(8, IORBX,LSPN) = XENG
!> XENG = -Exc -sum(WMSH*RHOG*COULOMB)/2 !E_sic per spin    
       RESULTS(9, IORBX,LSPN) = SIC_ENG
!> SIC_ENG = (-Exc -sum(WMSH*RHOG*COULOMB)/2) * (2/NSPN) !E_sic 
!>         -> This will be the value used in E(SIC).
       RESULTS(10, IORBX,LSPN) = ERGXL
!> Ex local from GETVLXC.
       RESULTS(11, IORBX,LSPN) = ERGCL
!> Ec local from GETVLXC.
       RESULTS(12, IORBX,LSPN) = ERGXN
!> Ex nonloc from GETVLXC.
       RESULTS(13, IORBX,LSPN) = ERGCN
!> Ec nonloc from GETVLXC.
C      DO JSPX=1,NSPN
C        DO IORB=1,NWFS(JSPX)
C          WRITE(INITOUT,800)IORB,(RESULTS(I,IORB,JSPX),I=1,13)
C        END DO
C      END DO
C      PRINT 800,
C    & JRBSIC,XCH1,XCH2,ASY,XEXX,-XCOU/2.,XDEL,XENG,SIC_ENG,
C    & ERGXL,ERGCL,ERGXN,ERGCN
C      ZSIC=SIC_ENG
 800   FORMAT('SIC COR (U,D,A,X,C,XOSIC,DLAM,DENG,TOT):',I5,15F11.4)
 1233  CONTINUE
C LB taken off due to parallelization
C      IF(LSPN.EQ.NSPN.AND.IORBX.EQ.MORB(NSPN)) THEN
C        CLOSE(41)
C        CLOSE(76)
C      END IF
C      CALL GTTIME(TIME1)
C      CALL FRMORB2(2,2)
C      CALL fmat_con(IORBX,LSPN)
C      CALL GTTIME(TIME2)
! Recovering mesh
       IF(MESH_SAV) THEN
         OPEN(99,FILE='VMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
         REWIND(99)
         READ(99) NMSH,JCALC
         READ(99)((RMSH(J,I),J=1,3),I=1,NMSH)
         READ(99)(WMSH(I),I=1,NMSH)
         CLOSE(99)
       END IF

! deallocate scaling factor
       IF(scaledsic) call deallocatersicfactor

C LB
C      CALL FLUSH(INITOUT)
C      CLOSE(INITOUT)
C LB
       INQUIRE(FILE='SICSTOP',EXIST=EXIST)
       IF(EXIST)CALL STOPIT
!      WRITE(6+IRANK,*)'SIC COUL:',SIC_COU,' EXCH:',SIC_EXO,
!     &                'ORBITAL INDEX:',IORBX
C      CALL FRMIORB(2,2)
       IRBSIC=0
       RETURN
       END
