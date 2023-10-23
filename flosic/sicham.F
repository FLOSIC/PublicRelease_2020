C ***************************************************************
C
       SUBROUTINE SICHAM(IFLO,MEQV,IFND)
C ORIGINAL VERSION BY MARK R PEDERSON (1988)
#ifndef NOCOMM
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
#else
       use FORRx, only : NFLO,KSPX,TMAT
       use SIC_ENERGY, only : SICENERGY
!       use NEWTIMES, only : TMWRD,TMNRM,CLWRD,CLNRM
!       use PTIME, only : TSETUP,TPOLYS,TFMTTM,TRUNIT,FLDNR,ASYMP,ACHRG
!     &  ,ATIME,TWASTED,COSTPROJ
       use TSTCHG, only : TOTCHG,COUSIC
!       use COUPDATA
       !use TMP1, only : CMB,RHOG
       use pot_dens, only : CMB=>COULOMB,RHOG
       !use TMP4, only : POT 
       use mixpot1,only : POTIN,POT=>POTOUT
       use mesh1,only : rmsh,nmsh,wmsh
       use debug1
       use coupdata1
       use common2,only : RIDT, BFCON, BFALP, N_BARE, N_CON,
     &   LSYMMAX, N_POS, NFNCT, IGGA, ISPN, NSPN
       use common3,only : RMAT,NGRP
       use common5,only : PSI, NWF, NWFS
       use common7,only : MODDEN, GAUSS_CUT
       use mpidat1,only : NPROC, NCALLED 
       include 'PARAMA2'
       INTEGER,INTENT(INOUT) :: IFLO,MEQV,IFND
       INTEGER :: I,I_LOC,I_LOC_BEG,I_POS,IALP,IATOM,IBASE,IBEG,
     &            ICELL,ICNT,
     &            ICNTR,IEND,IFNCT,IGP,IGRAD,IGRP,II,ILNG,ILOOP,IM,IMIN,
     &            IMSH,IPTS,IRHOG,ISHDUM,
     &            ISHELLA,ISITE,ISPD,ISPFAC,IST,ITOT,IWF,IWFS,IX,IZM,
     &            J,J_LOC,J_LOC_BEG,J_POS,JALP,JALP_BEG,JATOM,JBASE,
     &            JCELL,JFLO,JFNCT,JGP,JJCALL,JSHDUM,JSHELLA,JSITE,
     &            JSITE_BEG,JST,K_SITEI,K_SITEJ,KNTNAN,L,LC_LNG,LC_SHT,
     &            LI,LJ,LMSH,LNG,LSPN,LST,MAP,MCALL,MCNT,MD,MGRP,MINC,
     &            MM,MMSH,MSITES,MST,MTOT,MUI,MUIMAX,MUJ,MUJMAX,MXLG,
     &            MXPR,MXSIZ,NBLOCK,NCALL,NCELL,ND,NEEDED,NEEDED_PTS,
     &            NEEDED_SAV,NF1,NF2,NGRAD,NLOOP,
     &            NMAST,NMSH_SAV,NRWRD,NSET,NTID,NWRK
       REAL(8) :: TMWRD,TMNRM,CLWRD,CLNRM,
     &            TSETUP,TPOLYS,TFMTTM,TRUNIT,FLDNR,ASYMP,ACHRG,ARG,
     &            ATIME,TWASTED,COSTPROJ,
     &            ADDON,AI,AJ,ALCLNG,ALCSHT,ALPHAI,ALPHAJ,
     &            C,
     &            CHARGE,CHARGE_LONG,CHARGE_LR,CHARGE_NUM,COU,
     &            CU,DESIC,DVSIC,
     &            DIST,DMAT,DMIN,DN,ECT,ERRNAN,ERROR,EXC,FF,
     &            GAUSS_CUTI,GAUSS_CUTJ,
     &            OG,ORB_SCI,ORB_SIC,PI,PROD,PROEQV,PROIDT,PSIPSI,
     &            PSIPSI_NUM,RCALP,RLEN,RR_MAX,
     &            SAVINGS,SIC_COU,SIC_EXC,SS,TDD,
     &            TIME,TIME1,TIME2,TIMEP,TIMES_ALPHA,TIMES_SITES,
     &            TOT_CHG,TOT_SIC,TP1,TP2,VCDN,VCUP,VEX,VSIC,VXC,
     &            ZED,ZERO
#endif
       DIMENSION IFND(MX_GRP)
c  KAJ  get FLO info
c
       REAL*8 ::EC,DKF,ZET,DKS,RS,T,UU,VV,WW,ECRS,ECZET,ALFC
#ifndef NOCOMM
       COMMON/FORRx/NFLO,KSPX,TMAT(NDH,NDH,2)
       COMMON/SIC_ENERGY/SICENERGY
       COMMON/NEWTIMES/TMWRD,TMNRM,CLWRD,CLNRM
       COMMON/PTIME/TSETUP,TPOLYS,TFMTTM,TRUNIT,FLDNR,ASYMP,ACHRG
     &  ,ATIME,TWASTED,COSTPROJ
#endif
       PARAMETER (MXPR=MXPOISS)
       PARAMETER (MXLG=3)
#ifndef NOCOMM
       COMMON/TSTCHG/TOTCHG,COUSIC
       COMMON/COUPDATA/
     &  AIV(3,MXPR,MX_CNT+2,MXLG),AJV(3,MXPR,MX_CNT+2,MXLG)
     &  ,DMTV(10,10,MXPR,MX_CNT+2,MXLG),ALPIV(MXPR,MX_CNT+2,MXLG)
     &  ,ALPJV(MXPR,MX_CNT+2,MXLG),CENTER(3,MX_CNT+2)
     &  ,ADD(MAXUNSYM,MAXUNSYM,2),RVECI(3,MX_GRP),RVECJ(3,MX_GRP)
     &  ,NPAIRS(MX_CNT+2),IP(MX_CNT+2,MXLG)
#endif
       LOGICAL FIRST_ALPJ,CALLEDJ,NWRD,MWRD
#ifndef NOCOMM
       COMMON/TMP1/CMB    (MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
       COMMON/TMP4/POT    (MAX_PTS)!,RHOG(MAX_PTS,10,MXSPN)
#endif
       DIMENSION VSIC(MAX_PTS)
       DIMENSION IBEG(3),IEND(3)
       DIMENSION AI(3),AJ(3)
       DIMENSION SS(10,10),DMAT(10,10,2),C(3)
       DIMENSION OG(MAX_PTS,100),MAP(MX_GRP,1000)
       DIMENSION PROIDT(MAX_PTS),PROEQV(MAX_PTS)
       DIMENSION NEEDED(MX_GRP)
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA ZED/1.0D-30/
       GAUSS_CUT=1.0D30
       LMSH=NMSH ! NUMBER OF POINTS IN THE FUNDAMENTAL WEDGE
       NF1=1          
       NF2=MEQV        
       MGRP=0
       NMSH_SAV=NMSH
C EVENTUALLY GENERATE 100 RANDOM POINTS IN RMSH, SET NMSH=100 AND RUN
C THE GROUP ASSIGNER:
C 
C DETERMINE NEEDED REPLICAS OF RMSH:
       NMSH_SAV=NMSH
       MMSH=NMSH
       DO IMSH=1,NMSH
         WMSH(IMSH)=WMSH(IMSH)/NGRP
       END DO
       DO IGP=2,NGRP 
         DO IMSH=1,NMSH
           MMSH=MMSH+1
           WMSH(MMSH)=WMSH(IMSH)
           DO J=1,3
             RMSH(J,MMSH)=0.0D0
             DO L=1,3
               RMSH(J,MMSH)=RMSH(J,MMSH)+RMAT(J,L,IGP)*RMSH(L,IMSH)
             END DO
           END DO
         END DO
       END DO
       NMSH=MMSH !TOTAL NUMBER OF POINTS
       PRINT*,'POINTS BEFORE FLONASE:',NMSH
#ifdef NOCOMM
       ALLOCATE(AIV(3,MXPR,MX_CNT+2,MXLG))
       ALLOCATE(AJV(3,MXPR,MX_CNT+2,MXLG))
       ALLOCATE(DMTV(10,10,MXPR,MX_CNT+2,MXLG))
       ALLOCATE(ALPIV(MXPR,MX_CNT+2,MXLG))
       ALLOCATE(ALPJV(MXPR,MX_CNT+2,MXLG))
       ALLOCATE(CENTER(3,MX_CNT+2))
       ALLOCATE(ADD(MAXUNSYM,MAXUNSYM,2))
       ALLOCATE(RVECI(3,MX_GRP))
       ALLOCATE(RVECJ(3,MX_GRP))
       ALLOCATE(NPAIRS(MX_CNT+2))
       ALLOCATE(IP(MX_CNT+2,MXLG))
  
       call allocate_sic(.false.) 
       !reallocate to adjust the array size
       call allocate_sic(.true.)
#endif
C CHECK ORTHONORMALITY OF FLOS
C      MFLO=NFLO
C      NFLO=1
C      PRINT*,'ORTHOGONALITY CHECK'
C      CALL FLONASE(TIME1)
C      NFLO=MFLO
C      STOP
       DO JFLO=NF1,NF2 !THIS WILL BE DRIVE BY INFORMATION IN FRMORB
         IFLO=IFND(JFLO)
         NFLO=-IFLO
         POT=0.0D0
         PRINT*,'NFLO:',NFLO,KSPX,IFLO,NF1,NF2
         CALL FLONASE(TIME)
         TDD=0.0D0
         DO IMSH=1,MMSH
           TDD=TDD+ABS(POT(IMSH))
           OG(IMSH,IFLO)=POT(IMSH)
         END DO
         PRINT*,'ADD:',MMSH,JFLO,IFLO,NFLO,TDD
         IF(ABS(TDD).LT.1.0D-10)THEN
           STOP  !POT became zero
         END IF
       END DO
C      NMSH=NMSH_SAV
       DO JFLO=NF1,NF2
         IFLO=IFND(JFLO)
         PRINT*,'IFLO:',IFLO,MMSH
       END DO
       DO JFLO=NF1,NF2
         IFLO=IFND(JFLO)
         DO IGP=1,NGRP
           MAP(IGP,IFLO)=NGRP+1
         END DO
         DO JGP=1,NGRP              
           DO IGP=1,NGRP
             ERROR=0.0D0
             DO IMSH=1,LMSH
             ERROR=ERROR+ABS(OG(IMSH+(IGP-1)*LMSH,IFND(NF1))  !ALSO HERE
     &                      -OG(IMSH+(JGP-1)*LMSH,IFLO))
             END DO 
             IF(ERROR.LT.1.0D-5)THEN
               MAP(JGP,IFLO)=MIN(IGP,MAP(JGP,IFLO))
               MGRP=MAX(MAP(JGP,IFLO),MGRP)
               PRINT 610,LMSH,IFLO,IGP,JGP,ERROR
             END IF
           END DO
         END DO
         PRINT 615,IFLO,(MAP(JGP,IFLO),JGP=1,NGRP),0,0,IFND(NF1)
 610     FORMAT(4I5,G15.6)
 615     FORMAT('MAP IFLO:',I4,48I3)
       END DO
       PRINT*,'BACK FROM FLONASE'
C FIND NEEDED CELLS
       NCELL=0
       NEEDED=0
       DO JFLO=NF1,NF2
         IFLO=IFND(JFLO)
         PRINT 615,IFLO,(MAP(JGP,IFLO),JGP=1,NGRP),IFND(NF1)
         DO JGP=1,NGRP
           IMIN=1000
           DO ICELL=1,NCELL
             IMIN=MIN(ABS(MAP(JGP,IFLO)-NEEDED(ICELL)),IMIN)
           END DO
           IF(IMIN.NE.0)THEN
             NCELL=NCELL+1
             NEEDED(NCELL)=MAP(JGP,IFLO)
           END IF
         END DO
       END DO
       DO ICELL=1        ,NCELL
         DO JCELL=ICELL+1,NCELL
           IF(NEEDED(JCELL).LT.NEEDED(ICELL))THEN
             NEEDED_SAV=NEEDED(JCELL)
             NEEDED(JCELL)=NEEDED(ICELL)
             NEEDED(ICELL)=NEEDED_SAV
           END IF
         END DO
         PRINT*,ICELL,NEEDED(ICELL),' CELLS FOR THIS ADVENTURE'
       END DO
       call flush(6)
C FIGURE OUT THE FUNDAMENTAL WEDGE FOR FLO IFLO:
       PRINT*,'NEED MGRP OPERATIONS TO COMPLETE TASKS:',NCELL
       PRINT '(A,I10)','NUMBER OF INEQUIVALENT POINTS IN VMESH: ',NMSH
       OPEN(32,FILE='MESH_UNSYM',FORM='UNFORMATTED')
       WRITE(32)NMSH
       WRITE(32)RMSH,WMSH
       CLOSE(32)
       NMSH=0           
       DO ICELL=1,NCELL
         DO IMSH=1,LMSH
           NMSH=NMSH+1
           DO J=1,3
             RMSH(J,NMSH)=RMSH(J,IMSH+(NEEDED(ICELL)-1)*LMSH )
           END DO
           WMSH(  NMSH)=WMSH(  IMSH+(NEEDED(ICELL)-1)*LMSH )
         END DO
       END DO
       PRINT*,'POINTS REQUIRED 4 THIS SET OF FLOS:',NMSH
C PLACE POINTS FOR EACH CELL IN RMSH
       IF(NSPN*(MXPR/NSPN).NE.MXPR)THEN
         PRINT *,'MXPR MUST BE EVEN IN COUPOT1'
         CALL STOPIT
       END IF
       NFLO=-IFND(NF1) ! THIS WILL BE INFORMED FROM FRMORB
       TMWRD=0.0D0
       TMNRM=0.0D0
       IF (DEBUG) THEN
         PRINT *, 'GAUSS_CUT IN COUPOT'
         PRINT 9000,(GAUSS_CUT(IX),IX=1,MAX_IDENT)
 9000    FORMAT(4(1X,G15.8))
       END IF
       IZM=1                 
       IF(IZM.EQ.0)THEN      
C
C SETUP ATOMIC LOCATIONS AND BASIS SET INFORMATION
C
       CALL GTTIME(TIME1)
       LC_LNG=0
       LC_SHT=0
       TSETUP=0.0D0
       TPOLYS=0.0D0
       TFMTTM=0.0D0
       TRUNIT=0.0D0
       FLDNR=0.0D0
       ASYMP=0.0D0
       ACHRG=0.0D0
       ATIME=0.0D0
       TWASTED=0.0D0
       TIMEP=0.0D0
       COSTPROJ=0.0D0
       NCALL=0
       MCALL=0
       MCNT=0
       ISHELLA=0
       DMIN=1.0D30
       DO IFNCT=1,NFNCT
        DO I_POS=1,N_POS(IFNCT)
         ISHELLA=ISHELLA+1
         CALL GASITES(1,RIDT(1,ISHELLA),MTOT,RVECI,MSITES)
         DO ITOT=1,MTOT
          MCNT=MCNT+1
          CENTER(1,MCNT)=RVECI(1,ITOT)
          CENTER(2,MCNT)=RVECI(2,ITOT)
          CENTER(3,MCNT)=RVECI(3,ITOT)
          DIST=CENTER(1,MCNT)**2+CENTER(2,MCNT)**2+CENTER(3,MCNT)**2
          IF(DIST.LT.DMIN)DMIN=DIST
          NPAIRS(MCNT)=0
          IF (DEBUG) THEN
           PRINT 990,MCNT,(CENTER(J,MCNT),J=1,3)
  990      FORMAT(' PCENT:',I5,3G15.6)
          END IF
         END DO
        END DO
       END DO
C
       IF(SQRT(DMIN) .GT. 0.7D0)THEN
        MCNT=MCNT+1 
        CENTER(1,MCNT)=0.0D0
        CENTER(2,MCNT)=0.0D0
        CENTER(3,MCNT)=0.0D0
        NPAIRS(MCNT)=0
        IF (DEBUG) PRINT 990,MCNT,(CENTER(J,MCNT),J=1,3)
       END IF
       IF((MCNT-1).GT.MX_CNT)THEN
        PRINT*,'COUPOT: MX_CNT MUST BE AT LEAST:',MCNT-1
        CALL STOPIT
       END IF
c
#ifdef MPI
c
c MPI: send general data and initialize pot and rhog for all slaves
c
       CALL SENDDATA(103)
c
#endif
C
C ZERO POT AND RHOG
C
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
       DO IPTS=1,NMSH
        POT(IPTS)=0.0D0
        CMB(IPTS)=0.0D0
       END DO  
       IF (MODDEN .EQ. 1) THEN
        DO ISPN=1,NSPN
         DO IRHOG=1,NGRAD
          DO IPTS=1,NMSH
           RHOG(IPTS,IRHOG,ISPN)=0.0D0
          END DO
         END DO
        END DO
       END IF
C
       CHARGE=0.0D0
       CHARGE_LONG=0.0D0
       DO LNG=1,MXLG
        DO ICNTR=1,MCNT+1
         IP(ICNTR,LNG)=0
        END DO
       END DO
       ISITE=0
       ISHELLA=0
       DO 120 IFNCT=1,NFNCT
       DO 120 I_POS=1,N_POS(IFNCT)
        JSITE_BEG=ISITE
        ISHELLA=ISHELLA+1
        CALL OBINFO(1,RIDT(1,ISHELLA),RVECI,IST,ISHDUM)
C
C  TRANSLATE SHELLS TO ATOMS; GET GAUSS CUTOFF FOR ATOM I
C
        IATOM=ISHELLA
        GAUSS_CUTI=GAUSS_CUT(IATOM)
        DO 118 K_SITEI=1,IST
         ISITE=ISITE+1
c
c  KAJ  use RxRAVEL 
c
c        CALL UNRAVEL(IFNCT,ISHELLA,K_SITEI,RIDT(1,ISHELLA),
c    &                RVECI,LST,1)
C
           CALL RxRAVEL(IFNCT,ISHELLA,K_SITEI,RIDT(1,ISHELLA),
     &                  RVECI,LST,1)

         IF(LST.NE.IST)THEN
          PRINT *,'COUPOT: PROBLEM IN UNRAVEL'
          CALL STOPIT
         END IF
         I_LOC_BEG=0
         AI(1)=RVECI(1,K_SITEI)
         AI(2)=RVECI(2,K_SITEI)
         AI(3)=RVECI(3,K_SITEI)
         JSITE=JSITE_BEG
         JSHELLA=0
         DO 116 JFNCT=1,NFNCT
         DO 116 J_POS=1,N_POS(JFNCT)
          JSHELLA=JSHELLA+1
          IF(JSHELLA.LT.ISHELLA) GO TO 116
          CALL OBINFO(1,RIDT(1,JSHELLA),RVECJ,JST,JSHDUM)
C
C  TRANSLATE SHELLS TO ATOMS;  GET GAUSS CUT OFF  FOR ATOM J
C
          JATOM=JSHELLA
          GAUSS_CUTJ=GAUSS_CUT(JATOM)
          JJCALL=0 
          DO 114 K_SITEJ=1,JST
           CALLEDJ=.FALSE.
           JSITE=JSITE+1
           FIRST_ALPJ=.TRUE.
           IF(JSITE.LT.ISITE) GO TO 113
           J_LOC_BEG=0
           IF(JSITE.EQ.ISITE)J_LOC_BEG=I_LOC_BEG
           TIMES_SITES=1.0D0
           IF(JSITE.GT.ISITE)TIMES_SITES=2.0D0
           AJ(1)=RVECJ(1,K_SITEJ)
           AJ(2)=RVECJ(2,K_SITEJ)
           AJ(3)=RVECJ(3,K_SITEJ)
           DO 110 IALP=1,N_BARE(IFNCT)
            ALPHAI=BFALP(IALP,IFNCT)
            JALP_BEG=1
            IF(ISITE.EQ.JSITE)JALP_BEG=IALP
            DO 105 JALP=JALP_BEG,N_BARE(JFNCT)
             TIMES_ALPHA=1.0D0
             IF(ISITE.EQ.JSITE.AND.JALP.GT.IALP)TIMES_ALPHA=2.0D0
             ALPHAJ=BFALP(JALP,JFNCT)
             RCALP=1.0D0/(ALPHAI+ALPHAJ)
             ARG=ALPHAI*ALPHAJ*RCALP
     &        *((AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2)
             IF(ARG .GT.CUTEXP) GO TO 105
             LC_SHT=LC_SHT+1
             IF(ALPHAI.GT.GAUSS_CUTI.OR.ALPHAJ.GT.GAUSS_CUTJ)
     &         GO TO 105
             C(1)=(ALPHAI*AI(1)+ALPHAJ*AJ(1))*RCALP
             C(2)=(ALPHAI*AI(2)+ALPHAJ*AJ(2))*RCALP
             C(3)=(ALPHAI*AI(3)+ALPHAJ*AJ(3))*RCALP
             DMIN=1.0D30
             DO ICNT=1,MCNT
              DIST=SQRT((C(1)-CENTER(1,ICNT))**2+
     &                  (C(2)-CENTER(2,ICNT))**2+
     &                  (C(3)-CENTER(3,ICNT))**2)
              IF(DIST.LT.DMIN)THEN
               DMIN=DIST
               MINC=ICNT
              END IF
             END DO
C
C BLOCK 1:
C
             RLEN=(ALPHAI+ALPHAJ)*DMIN**2
             IF(RLEN .GT. 0.5D0)MINC=MCNT+1
             IF(ALPHAI+ALPHAJ .LT. 3.0D0)MINC=MCNT+1
             LNG=1
             IF(ALPHAI+ALPHAJ.GT. 10.0D0)LNG=2
             IF(ALPHAI+ALPHAJ.GT. 40.0D0)LNG=3
C
C END OF BLOCK 2:
C
             IF(.NOT.CALLEDJ)THEN
              JJCALL=JJCALL+1
              IF(JJCALL.EQ.1.AND.K_SITEJ.NE.1)THEN
c              CALL UNRAVEL(JFNCT,JSHELLA,1,RIDT(1,JSHELLA),
c    &                      RVECJ,MST,2)
c KAJ  Rxravel
               CALL RxRAVEL(JFNCT,JSHELLA,1,RIDT(1,JSHELLA),
     &                      RVECJ,MST,2)
               IF(MST.NE.JST)THEN
                PRINT *,'COUPOT: PROBLEM IN UNRAVEL'
                CALL STOPIT
               END IF
              END IF
c KAJ  Rxravel
c             CALL UNRAVEL(JFNCT,JSHELLA,K_SITEJ,RIDT(1,JSHELLA),
c    &                     RVECJ,MST,2)
              CALL RxRAVEL(JFNCT,JSHELLA,K_SITEJ,RIDT(1,JSHELLA),
     &                     RVECJ,MST,2)
              IF(MST.NE.JST)THEN
               PRINT *,'COUPOT: PROBLEM IN UNRAVEL'
               CALL STOPIT
              END IF
              CALLEDJ=.TRUE.
             END IF
C
C END OF NEW SETUP...
C
C
C CALCULATE DENSITY MATRIX:
C

             DO ISPD=1,NSPN
              DO MUI=1,10
               DO MUJ=1,10
                DMAT(MUJ,MUI,ISPD)=0.0D0
               END DO
              END DO
             END DO
             MUIMAX=0
             DO LI =0,LSYMMAX(IFNCT)
              DO MUI=IBEG(LI+1),IEND(LI+1)
               DO IBASE=1,N_CON(LI+1,IFNCT)
                IF(ABS(BFCON(IALP,IBASE,LI+1,IFNCT)).GT.ZED
     &            .AND.MUI.GT.MUIMAX) MUIMAX=MUI
               END DO
              END DO
             END DO
             MUJMAX=0
             DO LJ =0,LSYMMAX(JFNCT)
              DO MUJ=IBEG(LJ+1),IEND(LJ+1)
               DO JBASE=1,N_CON(LJ+1,JFNCT)
                IF(ABS(BFCON(JALP,JBASE,LJ+1,JFNCT)).GT.ZED
     &            .AND.MUJ.GT.MUJMAX) MUJMAX=MUJ
               END DO
              END DO
             END DO
C
             IF(FIRST_ALPJ)THEN
              I_LOC    =I_LOC_BEG
              DO LI =0,LSYMMAX(IFNCT)
               DO MUI=IBEG(LI+1),IEND(LI+1)
                DO IBASE=1,N_CON(LI+1,IFNCT)
                 I_LOC=I_LOC+1
                 J_LOC=J_LOC_BEG
                 DO LJ =0,LSYMMAX(JFNCT)
                  DO MUJ=IBEG(LJ+1),IEND(LJ+1)
                   DO JBASE=1,N_CON(LJ+1,JFNCT)
                    J_LOC=J_LOC+1
                    IWF=0
                    DO ISPN=1,NSPN
                     ADD(I_LOC,J_LOC,ISPN)=0.0D0
                     DO IWFS=1,NWFS(ISPN)
                      IWF=IWF+1
                       ADD(I_LOC,J_LOC,ISPN)=ADD(I_LOC,J_LOC,ISPN)
     &                +PSI(I_LOC,IWF,1)*PSI(J_LOC,IWF,2)
                     END DO
                    END DO
                   END DO
                  END DO
                 END DO
                END DO
               END DO
              END DO
              FIRST_ALPJ=.FALSE.
              IF(I_LOC.GT.MAXUNSYM)THEN
               PRINT*,'COUPOT: MAXUNSYM MUST BE AT LEAST:',I_LOC
               CALL STOPIT
              END IF
              IF(J_LOC.GT.MAXUNSYM)THEN
               PRINT*,'COUPOT: MAXUNSYM MUST BE AT LEAST:',J_LOC
               CALL STOPIT
              END IF
             END IF
             I_LOC    =I_LOC_BEG
             DO LI =0,LSYMMAX(IFNCT)
              DO MUI=IBEG(LI+1),IEND(LI+1)
               DO IBASE=1,N_CON(LI+1,IFNCT)
                I_LOC=I_LOC+1
                J_LOC=J_LOC_BEG
                DO LJ =0,LSYMMAX(JFNCT)
                 DO MUJ=IBEG(LJ+1),IEND(LJ+1)
                  DO JBASE=1,N_CON(LJ+1,JFNCT)
                   J_LOC=J_LOC+1
                   PROD=BFCON(IALP,IBASE,LI+1,IFNCT)
     &                 *BFCON(JALP,JBASE,LJ+1,JFNCT)
                   DO ISPN=1,NSPN
                    DMAT(MUI,MUJ,ISPN)=DMAT(MUI,MUJ,ISPN)
     &                                +PROD*ADD(I_LOC,J_LOC,ISPN)
                   END DO
                  END DO
                 END DO
                END DO
               END DO
              END DO
             END DO
C
             CALL OVMXSF(ALPHAI,ALPHAJ,AI,AJ,SS)
             ADDON=0.0D0
             DO ISPN=1,NSPN
              DO MUJ=1,MUJMAX
               DO MUI=1,MUIMAX
                ADDON=ADDON+DMAT(MUI,MUJ,ISPN)*SS(MUI,MUJ)
     &               *TIMES_SITES*TIMES_ALPHA
               END DO
              END DO
             END DO
             CHARGE=CHARGE+ADDON
             DO ISPN=1,NSPN
              DO MUJ=1,MUJMAX
               DO MUI=1,MUIMAX
                DMAT(MUI,MUJ,ISPN)=DMAT(MUI,MUJ,ISPN)
     &                            *TIMES_SITES*TIMES_ALPHA
               END DO
              END DO
             END DO
             CHARGE_LONG=CHARGE_LONG+ADDON
C
             DO 10 ISPN=1,NSPN
              IP(MINC,LNG)=IP(MINC,LNG)+1
              NPAIRS(MINC)=NPAIRS(MINC)+1
              LC_LNG=LC_LNG+1
              ND=MUIMAX
              MD=MUJMAX
              ALPIV(IP(MINC,LNG),MINC,LNG)=ALPHAI
              ALPJV(IP(MINC,LNG),MINC,LNG)=ALPHAJ
              AIV(1,IP(MINC,LNG),MINC,LNG)=AI(1)
              AIV(2,IP(MINC,LNG),MINC,LNG)=AI(2)
              AIV(3,IP(MINC,LNG),MINC,LNG)=AI(3)
              AJV(1,IP(MINC,LNG),MINC,LNG)=AJ(1)
              AJV(2,IP(MINC,LNG),MINC,LNG)=AJ(2)
              AJV(3,IP(MINC,LNG),MINC,LNG)=AJ(3)
              DO MUJ=1,10
               DO MUI=1,10
                DMTV(MUI,MUJ,IP(MINC,LNG),MINC,LNG)=DMAT(MUI,MUJ,ISPN)
               END DO
              END DO
   10        CONTINUE
             MXSIZ=MXPR
c
#ifdef MPI
c
             IF (MINC.EQ.MCNT+1) MXSIZ=MAX(1,MXPR/10)
c
#endif
c
             IF (IP(MINC,LNG).EQ.MXSIZ) THEN
              IF (IP(MINC,LNG).NE.NSPN) THEN
               ND=10
               MD=10
              END IF
              NCALL=NCALL+1
              MCALL=MCALL+1
              CALL GTTIME(TP1)
   
              IF(MINC.EQ.MCNT+1)THEN
               NWRD=.TRUE.
              ELSE
               NWRD=.FALSE.
              END IF
c
#ifndef MPI
c
              CALL POISSON2(NWRD,IP(MINC,LNG),ND,MD,
     &                      ALPIV(1,MINC,LNG),AIV(1,1,MINC,LNG),
     &                      ALPJV(1,MINC,LNG),AJV(1,1,MINC,LNG),
     &                      DMTV(1,1,1,MINC,LNG))
c
#else
c
c dvp+mrp 10/98: let master do some work when all children are busy
c                first deal with single weirdos, then in larger blocks
c                with non-weirdos. The block size for non-weirdos is
c                reduced as the number of processors increases
c
              IF (NPROC .EQ. 0) THEN
               CALL POISSON2(NWRD,IP(MINC,LNG),ND,MD,
     &                       ALPIV(1,MINC,LNG),AIV(1,1,MINC,LNG),
     &                       ALPJV(1,MINC,LNG),AJV(1,1,MINC,LNG),
     &                       DMTV(1,1,1,MINC,LNG))
              ELSE
               NMAST= MXPR/(4*NPROC*NSPN)
               IF(NMAST.GT.4)NMAST= 0              
               IF (NMAST .LE. 0) THEN
#ifdef NOCOMM
                CALL CKWORKER(2,NTID)
#else
                CALL CKCHILD(2,NTID)
#endif
                GOTO 20
               END IF
               DO ICNT= MCNT+1, 1, -1
                MWRD= (ICNT .EQ. MCNT+1)
                IF (MWRD) THEN
                 NBLOCK= NSPN
                ELSE
                 NBLOCK= NSPN*NMAST
                END IF
                DO ILNG= 1,3
                 NLOOP= (IP(ICNT,ILNG)+NBLOCK-1)/NBLOCK
                 DO ILOOP= 1,NLOOP
#ifdef NOCOMM
                  CALL CKWORKER(1,NTID)
#else
                  CALL CKCHILD(1,NTID)
#endif
                  IF (NCALLED .NE. NPROC) GOTO 20
                  NSET= MIN(IP(ICNT,ILNG),NBLOCK)
                  I= IP(ICNT,ILNG)-NSET+1
                  CALL POISSON2(MWRD,NSET,10,10,
     &                          ALPIV(I,ICNT,ILNG),AIV(1,I,ICNT,ILNG),
     &                          ALPJV(I,ICNT,ILNG),AJV(1,I,ICNT,ILNG),
     &                          DMTV(1,1,I,ICNT,ILNG))
                  IP(ICNT,ILNG)= IP(ICNT,ILNG)-NSET
                 END DO
                END DO
               END DO
   20          IF (IP(MINC,LNG) .GT. 0) THEN
                CALL PAMPOISS(1,ND,MD,MINC,LNG,NWRD,POT,RHOG)
               END IF
              END IF
c
#endif
c
              CALL GTTIME(TP2)
              TIMEP=TIMEP+TP2-TP1
              IP(MINC,LNG)=0
             END IF
  105       CONTINUE
  110      CONTINUE
  113     CONTINUE
  114     CONTINUE
  116    CONTINUE
  118   CONTINUE
  120  CONTINUE
C
C EMPTY TABLE OF PAIRS
C
       DO 800 LNG=1,MXLG
        DO 790 MINC=1,MCNT+1
         NCALL=NCALL+IP(MINC,LNG)
         IF(IP(MINC,LNG).NE.0)THEN
          MCALL=MCALL+1
          ND=10
          MD=10
          CALL GTTIME(TP1)
          IF(MINC.EQ.MCNT+1)THEN
           NWRD=.TRUE.
          ELSE
           NWRD=.FALSE.
          END IF
c
#ifndef MPI
c
          CALL POISSON2(NWRD,IP(MINC,LNG),ND,MD,
     &                  ALPIV(1,MINC,LNG),AIV(1,1,MINC,LNG),
     &                  ALPJV(1,MINC,LNG),AJV(1,1,MINC,LNG),
     &                  DMTV(1,1,1,MINC,LNG))
c
#else
c
c dvp+mrp 10/98: let master do some work when all children are busy
c                first deal with single weirdos, then in larger blocks
c                with non-weirdos. The block size for non-weirdos is
c                reduced as the number of processors increases
c
          IF (NPROC .EQ. 0) THEN
           CALL POISSON2(NWRD,IP(MINC,LNG),ND,MD,
     &                   ALPIV(1,MINC,LNG),AIV(1,1,MINC,LNG),
     &                   ALPJV(1,MINC,LNG),AJV(1,1,MINC,LNG),
     &                   DMTV(1,1,1,MINC,LNG))
          ELSE
           NMAST= MXPR/(4*NPROC*NSPN)
           IF(NMAST.GT.4)NMAST=0
C FOR MANY PROCESSOR CALCULATIONS ELIMINATE POSSIBILITY
C OF ALLOWING MASTER TO DO WORK...
C
C
C DISCARD CKCHILD(2,NTID) MODE.....FOR NOW.  12/30/2003
C PUT CKCHILD BACK IN                        01/07/2004
C          NMAST=0
           IF (NMAST .LE. 0) THEN
#ifdef NOCOMM
            CALL CKWORKER(2,NTID)
#else
            CALL CKCHILD(2,NTID)
#endif
            GOTO 30
           END IF
           DO ICNT= MCNT+1, 1, -1
            MWRD= (ICNT .EQ. MCNT+1)
            IF (MWRD) THEN
             NBLOCK= NSPN
            ELSE
             NBLOCK= NSPN*NMAST
            END IF
            DO ILNG= 1,3
             NLOOP= (IP(ICNT,ILNG)+NBLOCK-1)/NBLOCK
             DO ILOOP= 1,NLOOP
#ifdef NOCOMM
              CALL CKWORKER(1,NTID)
#else
              CALL CKCHILD(1,NTID)
#endif
              IF (NCALLED .NE. NPROC) GOTO 30
              NSET= MIN(IP(ICNT,ILNG),NBLOCK)
              I= IP(ICNT,ILNG)-NSET+1
              CALL POISSON2(MWRD,NSET,10,10,
     &                      ALPIV(I,ICNT,ILNG),AIV(1,I,ICNT,ILNG),
     &                      ALPJV(I,ICNT,ILNG),AJV(1,I,ICNT,ILNG),
     &                      DMTV(1,1,I,ICNT,ILNG))
              IP(ICNT,ILNG)= IP(ICNT,ILNG)-NSET
             END DO
            END DO
           END DO
   30      IF (IP(MINC,LNG) .GT. 0) THEN
            CALL PAMPOISS(1,ND,MD,MINC,LNG,NWRD,POT,RHOG)
            IP(MINC,LNG)= 0
           END IF
          END IF
c
#endif
c
          CALL GTTIME(TP2)
          TIMEP=TIMEP+TP2-TP1
         END IF
  790   CONTINUE
  800  CONTINUE
c
#ifdef MPI
c
       CALL PAMPOISS(2,ND,MD,MINC,LNG,NWRD,POT,RHOG)   !NRWRD -> NWRD
c
#endif
c
       ELSE    
        CONTINUE
        RHOG=0.0D0
        CMB=0.0D0
        CALL COUPOT1
        CMB=CMB*(NSPN/2.0D0)
        CHARGE=0.0D0
        DO IMSH=1,NMSH
          CHARGE=CHARGE+RHOG(IMSH,1,1)*WMSH(IMSH)
        END DO
        PRINT*,'CHARGE_T ON FUNDAMENTAL WEDGE:',CHARGE,NMSH
        OPEN(32,FILE='MESH_UNSYM',FORM='UNFORMATTED')
        READ (32)NMSH
        READ (32)RMSH,WMSH
        CLOSE(32)
        DO ICELL=NCELL,1,-1
         PRINT*,'ICELL, NEEDED:',ICELL,NEEDED(ICELL)
         DO IMSH=1,LMSH
           CMB (IMSH+(NEEDED(ICELL)-1)*LMSH)    
     &            =CMB (IMSH+(ICELL-1)*LMSH)    
           IF(ICELL.LT.NEEDED(ICELL))THEN
             CMB (IMSH+(ICELL-1)*LMSH)=0.0D0
           END IF
           DO LSPN=1,NSPN
             DO IGRAD=1,10
               RHOG(IMSH+(NEEDED(ICELL)-1)*LMSH,IGRAD,LSPN)
     &        =RHOG(IMSH+(ICELL-1        )*LMSH,IGRAD,LSPN)
               IF(ICELL.LT.NEEDED(ICELL))THEN
               RHOG(IMSH+(ICELL-1        )*LMSH,IGRAD,LSPN)=0.0D0
               END IF
             END DO
           END DO
         END DO
        END DO
        DO IGP=1,NGRP     
          DO IMSH=1,LMSH
            CMB (IMSH+(IGP-1)          *LMSH           )=
     &      CMB (IMSH+(MAP(IGP,IFND(1))-1)*LMSH           )
            DO LSPN=1,NSPN
              DO IGRAD=1,10
                RHOG(IMSH+(IGP       -1)*LMSH,IGRAD,LSPN)=
     &          RHOG(IMSH+(MAP(IGP,IFND(1))-1)*LMSH,IGRAD,LSPN) 
              END DO
            END DO
          END DO
        END DO
        JFLO=1
        IFLO=IFND(JFLO)
        IF(IFLO.EQ.1.AND.KSPX.EQ.1)TOT_CHG=0.0D0 !KW
        IF(IFLO.EQ.1.AND.KSPX.EQ.1)TOT_SIC=0.0D0 !KW          
        IF(JFLO.EQ.1)SIC_COU=0.0D0
        IF(JFLO.EQ.1)SIC_EXC=0.0D0
        IF(JFLO.EQ.1)ORB_SIC=0.0D0
        IF(JFLO.EQ.1)ORB_SCI=0.0D0
        CHARGE=0.0D0
        DESIC=0.0D0
        DVSIC=0.0D0 
        PI=4.0D0*ATAN(1.0D0)
        DO IMSH=1,LMSH*NGRP
         DN=RHOG(IMSH,1,KSPX)
         CHARGE=CHARGE+RHOG(IMSH,1,KSPX)*WMSH(IMSH)
         ZET=1.0D0
         DKF=(DN*3.0D0*PI**2)**(1./3.)
         DKS=SQRT(1.2732395447351627D0*DKF)
         RS=1.9191582926775128D0/DKF
         T=0.0D0
         UU=0.0D0
         VV=0.0D0
         WW=0.0D0
         VCUP=0.0D0
         VCDN=0.0D0
         EC=0.0D0
         VEX=0.0D0
         IF(DN.GT.1.0D-20)THEN        
           CALL PW91LC(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
           VEX=-(6.0D0*DN/PI)**(1./3.)
         END IF
         VSIC(IMSH)=-(CMB(IMSH)+VEX+VCUP)
         DVSIC=DVSIC+VSIC(IMSH)*DN*WMSH(IMSH)
         DESIC=DESIC-(0.5*CMB(IMSH)+0.75*VEX+EC)*DN*WMSH(IMSH)           
         SIC_COU=SIC_COU-0.5*CMB(IMSH)*DN*WMSH(IMSH)*(NF2-NF1+1)
         SIC_EXC=SIC_EXC-(0.75*VEX+EC)*DN*WMSH(IMSH)*(NF2-NF1+1)
        END DO
        ORB_SCI=DVSIC
        TOT_SIC=TOT_SIC+SIC_COU+SIC_EXC
        TOT_CHG=TOT_CHG+CHARGE*(NF2-NF1+1)
        PRINT 317,CHARGE,DVSIC,DESIC,TOT_CHG,TOT_SIC                                     
 317    FORMAT('THE TALLY:',10F15.6)
       END IF
       SICENERGY=TOT_SIC*(2/NSPN)
       CALL FLUSH(6)
       PRINT*,'TOT_SIC AND SICENERGY:',SICENERGY
 50    FORMAT(' WEDGE NEW COUSIC:',2I8,8F15.6)
 52    FORMAT(' WEDGE TOT_SIC, CHRGE, IC_COU, SIC_EXC, ORB_SIC, ORBSCI')
C      END DO
C      PRINT 50,LMSH,IFLO,TOT_SIC,TOT_CHG,
C    &     SIC_COU,SIC_EXC,ORB_SIC,ORB_SCI
C      FF=27.2116
C      PRINT 51,LMSH,IFLO,ORB_SIC,ORB_SCI,ORB_SIC*FF,ORB_SCI*FF
C51    FORMAT(' ORBITAL ENERGY AND SHIFT:',2I8,8F15.6)
CCC END OF BUG
       RR_MAX=0.0D0
       PSIPSI_NUM=0.0D0
       ZERO=0.0D0
! 20    FORMAT(10F12.5)
       PRINT*,'NMSH BEFORE FLONASE:',NMSH
       NFLO=-IFND(1)
       CALL FLONASE(TIME) ! RETRIEVE INEQUIVALENT FLO
       PRINT*,'NMSH AFTER  FLONASE:',NMSH
       DO 150 IPTS=1,NMSH
       PSIPSI_NUM=PSIPSI_NUM+POT(IPTS)*POT(IPTS)*WMSH(IPTS)
 150   CONTINUE
C            PSIPSI_NUM=(PSIPSI_NUM*NGRP)/MGRP
C            CHARGE_NUM=(CHARGE_NUM*NGRP)/MGRP
       PRINT*,"TST: PSI/CHG (NUM):",(NMSH/LMSH),PSIPSI_NUM,CHARGE_NUM
       PRINT*,"TST: LONG  CHARGE IN COUPOT: ",CHARGE_LONG
       PRINT*,"TOTAL CHARGE IN COUPOT:      ",CHARGE
       CALL FLUSH(6)
C      MM=1 !SIC Scisster
       MM=2 !PRP JCP and Support of Multiplicative Hamiltonian
       IF(MM.EQ.1)ORB_SCI=SQRT(ABS(ORB_SCI))
       IF(MM.EQ.2)THEN
         ORB_SCI=1.0D0            
C        IF(NF2-NF1.NE.0)THEN
C          PRINT*,'COMPLETE VSIC MAPPING'
C          CALL STOPIT
C        END IF
       END IF
C NOTE THAT POT=FLO
       NWRK=NMSH
       DO JFLO=NF1,NF2
        DO IM=1,MM
         PROIDT=0.0D0
         DO IMSH=1,NWRK
          PROIDT(IMSH)=POT(IMSH)*WMSH(IMSH)*ORB_SCI !W(I)|FLO_I><FLO_I|VI|FLO_I>
         END DO
         IF(IM.EQ.2)THEN
          DO IMSH=1,NWRK
           PROIDT(IMSH)=-PROIDT(IMSH)*VSIC(IMSH) !To match sign convention in 
          END DO 
         END IF
         IFLO=IFND(JFLO)
         DO JGP=1,NGRP   !MOST RECENT CHANGE
          DO IMSH=1,LMSH         
           PROEQV(IMSH+(JGP-1)          *LMSH)
     &    =PROIDT(IMSH+(MAP(JGP,IFLO)-1)*LMSH)
          END DO
         END DO
         NMSH=NGRP*LMSH
         II=IM
         CALL APATCH(II,MM,PROEQV)
        END DO
       END DO
       OPEN(60,FILE='VMOLD',FORM='UNFORMATTED')
       READ(60)NMSH
       READ(60)((RMSH(J,IMSH),J=1,3),IMSH=1,NMSH)
       READ(60) (WMSH(IMSH),IMSH=1,NMSH)
       CLOSE(60) 
C CALCULATE EXCHANGE CORRELATION POTENTIAL FOR ORBITAL
       CHARGE=CHARGE!*ISPFAC
       CHARGE_LONG=CHARGE_LONG!*ISPFAC
       IF(ABS(CHARGE_LR-CHARGE) .GT. 0.1D0)THEN
        PRINT *,'************** WARNING ****************'
        PRINT *,'CHARGE VIA ASYMPTOTIC FORMULA: ',CHARGE_LR
       END IF
       CALL GTTIME(TIME2)
       IF (DEBUG) THEN
        PRINT*,'TIME FOR COULOMB POTENTIAL,NMSH:',TIME2-TIME1,NMSH
        PRINT*,'POISSON DIAG: PT   :',TIMEP
        PRINT*,'POISSON DIAG: CT   :',TIME2-TIME1-TIMEP
        PRINT*,'POISSON DIAG, SETUP:',TSETUP
        PRINT*,'POISSON DIAG, POLYS:',TPOLYS
        PRINT*,'POISSON DIAG, FMTTM:',TFMTTM
        PRINT*,'POISSON DIAG, RUNIT:',TRUNIT
        PRINT*,'POISSON DIAG, ATIME:',ATIME
        PRINT*,'POISSON DIAG, WASTED:',TWASTED
        PRINT*,'POISSON DIAG, FLDNR:',FLDNR,FLDNR/(FLDNR+ASYMP+0.001D0)
        PRINT*,'POISSON DIAG, ASYMP:',ASYMP,ASYMP/(FLDNR+ASYMP+0.001D0)
        PRINT*,'POISSON DIAG, ACHRG:',ACHRG,CHARGE_LONG
        PRINT*,'POISSON DIAG, PRTME:',COSTPROJ
        PRINT 997,TMWRD,TMNRM,CLWRD,CLNRM,
     &            TMWRD/(CLWRD+0.001D0),TMNRM/(CLNRM+0.001D0)
        DO ICNT=1,MCNT+1
         PRINT 998,(CENTER(J,ICNT),J=1,3),NPAIRS(ICNT),NCALL,MCALL
        END DO
       END IF
 998   FORMAT(' POISSON:',3(1X,G18.6)/9X,3(1X,I18))
 997   FORMAT(' NEW TIMES:',6G12.4)
#ifdef NOCOMM
       print *,'deallocate aiv'
       IF(ALLOCATED(AIV))    DEALLOCATE(AIV)
       IF(ALLOCATED(AJV))    DEALLOCATE(AJV)
       IF(ALLOCATED(DMTV))   DEALLOCATE(DMTV)
       IF(ALLOCATED(ALPIV))  DEALLOCATE(ALPIV)
       IF(ALLOCATED(ALPJV))  DEALLOCATE(ALPJV)
       IF(ALLOCATED(CENTER)) DEALLOCATE(CENTER)
       IF(ALLOCATED(ADD))    DEALLOCATE(ADD)
       IF(ALLOCATED(RVECI))  DEALLOCATE(RVECI)
       IF(ALLOCATED(RVECJ))  DEALLOCATE(RVECJ)
       IF(ALLOCATED(NPAIRS)) DEALLOCATE(NPAIRS)
       IF(ALLOCATED(IP))     DEALLOCATE(IP)
       call allocate_sic(.false.) 
       !reallocate to adjust the array size
       call allocate_sic(.true.)
#endif
       PRINT *,'MADE IT END' 
       RETURN   
       END
C
