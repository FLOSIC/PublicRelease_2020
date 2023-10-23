C UTEP Electronic Structure Lab (2020)
C
C ***************************************************************
C
       SUBROUTINE COUPOT_DMAT
C ORIGINAL VERSION BY MARK R PEDERSON (1988)
! Version that uses density matrix LBasurto (2017)
       use global_inputs,only : dosjnt1,INBAS
       use pot_dens,only : COULOMB,RHOG
       use coupdata1
       use mesh1,only : rmsh,nmsh
       use debug1
       use common2,only : RIDT, BFCON, BFALP, N_BARE, N_CON,
     &   LSYMMAX, N_POS, NFNCT, IGGA, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI, NWF, NWFS
       use common7,only : MODDEN, GAUSS_CUT
       use mpidat1,only : NPROC, NCALLED
       use den_mat,only : DMAT_LOCAL
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:47 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I_LOC, I_LOC_BEG, I_POS, IALP, IATOM, IBASE, IBEG,
     & ICNT, IEND, IERR, IFNCT, ILNG, ILOOP, IPTS, ISHDUM, ISHELLA,
     & ISITE, ISPFAC, IST, ITOT, IWF, IWFS, IX, J, J_LOC, J_LOC_BEG,
     & JALP, JALP_BEG, JATOM, JBASE, JFNCT, JJCALL, JSHDUM, JSHELLA,
     & JSHELLATMP, JSITE, JSITE_BEG, JST, K_SITEI, K_SITEJ, LC_LNG,
     & LC_SHT, LI, LJ, LNG, LST, MCALL, MCNT, MD, MINC, MSITES, MST,
     & MTOT, MUI, MUIMAX, MUJ, MUJMAX, MXLG, MXPR, MXSIZ, NBLOCK,
     & NCALL, ND, NGRAD, NLOOP, NMAST, NRWRD, NSET, NTID
       REAL*8 :: SYMBOL , ACHRG, ADDON, AI, AJ, ALCLNG, ALCSHT, ALPHAI,
     & ALPHAJ, ARG, ASYMP, ATIME, C, CHARGE, CHARGE_LONG, CHARGE_LR,
     & CLNRM, CLWRD, COSTPROJ, DIST, DMAT2, DMIN, FLDNR,GAUSS_CUTI,
     & GAUSS_CUTJ, PROD, RARG, RCALP, RLEN, RR_MAX, RRR, SAVINGS, SS,
     & TFMTTM, TIME1, TIME2, TIMEP, TIMES_ALPHA, TIMES_SITES, TMNRM,
     & TMWRD, TP1, TP2, TPOLYS, TRUNIT, TSETUP, TWASTED, ZED
       SAVE
!       COMMON/NEWTIMES/TMWRD,TMNRM,CLWRD,CLNRM
!       COMMON/PTIME/TSETUP,TPOLYS,TFMTTM,TRUNIT,FLDNR,ASYMP,ACHRG
!     &  ,ATIME,TWASTED,COSTPROJ
       PARAMETER (MXPR=MXPOISS)
       PARAMETER (MXLG=3)
!       COMMON/COUPDATA/
!     &  AIV(3,MXPR,MX_CNT+2,MXLG),AJV(3,MXPR,MX_CNT+2,MXLG)
!     &  ,DMTV(10,10,MXPR,MX_CNT+2,MXLG),ALPIV(MXPR,MX_CNT+2,MXLG)
!     &  ,ALPJV(MXPR,MX_CNT+2,MXLG),CENTER(3,MX_CNT+2)
!     &  ,ADD(MAXUNSYM,MAXUNSYM,2),RVECI(3,MX_GRP),RVECJ(3,MX_GRP)
!     &  ,NPAIRS(MX_CNT+2),IP(MX_CNT+2,MXLG)
       LOGICAL FIRST_ALPJ,CALLEDJ,NWRD,MWRD,EXIST
       CHARACTER LETTER
!       COMMON/TMP1/ACOULOMB(MAX_PTS),ARHOG(MAX_PTS,KRHOG,MXSPN)
       DIMENSION IBEG(3),IEND(3)
       DIMENSION AI(3),AJ(3)
       DIMENSION SS(10,10),DMAT2(50,50,2),C(3)
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA ZED/1.0D-30/
C
       PRINT '(A,I10)','NUMBER OF POINTS IN COUPOT: ',NMSH
       IF(NSPN*(MXPR/NSPN).NE.MXPR)THEN
        write(6,*)'MXPR MUST BE EVEN IN COUPOT1'
        CALL STOPIT

       END IF
       TMWRD=0.0D0
       TMNRM=0.0D0
       IF (DEBUG) THEN
        write(6,*) 'GAUSS_CUT IN COUPOT'
        PRINT 9000,(GAUSS_CUT(IX),IX=1,MAX_IDENT)
 9000   FORMAT(4(1X,G15.8))
       END IF
!
! ALLOCATE ARRAYS FROM MODULE COUPDATA1
!
       IF(.NOT.ALLOCATED(AIV))THEN
         ALLOCATE(AIV(3,MXPR,MX_CNT+2,MXLG),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING AIV'
       IF(.NOT.ALLOCATED(AJV))THEN
         ALLOCATE(AJV(3,MXPR,MX_CNT+2,MXLG),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING AJV'
       IF(.NOT.ALLOCATED(DMTV))THEN
         ALLOCATE(DMTV(10,10,MXPR,MX_CNT+2,MXLG),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING DMTV'
       IF(.NOT.ALLOCATED(ALPIV))THEN
         ALLOCATE(ALPIV(MXPR,MX_CNT+2,MXLG),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING ALPIV'
       IF(.NOT.ALLOCATED(ALPJV))THEN
         ALLOCATE(ALPJV(MXPR,MX_CNT+2,MXLG),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING ALPJV'
       IF(.NOT.ALLOCATED(CENTER))THEN
         ALLOCATE(CENTER(3,MX_CNT+2),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING CENTER'
       IF(.NOT.ALLOCATED(ADD))THEN
         ALLOCATE(ADD(MAXUNSYM,MAXUNSYM,2),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING ADD'
       IF(.NOT.ALLOCATED(RVECI))THEN
         ALLOCATE(RVECI(3,MX_GRP),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING RVECI'
       IF(.NOT.ALLOCATED(RVECJ))THEN
         ALLOCATE(RVECJ(3,MX_GRP),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING RVECJ'
       IF(.NOT.ALLOCATED(NPAIRS))THEN
         ALLOCATE(NPAIRS(MX_CNT+2),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING NPAIRS'
       IF(.NOT.ALLOCATED(IP))THEN
         ALLOCATE(IP(MX_CNT+2,MXLG),STAT=IERR)
       ENDIF
       IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR ALLOCATING IP'
!
! SETUP ATOMIC LOCATIONS AND BASIS SET INFORMATION
!
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
          CENTER(:,MCNT)=RVECI(:,ITOT)
          !CENTER(1,MCNT)=RVECI(1,ITOT)
          !CENTER(2,MCNT)=RVECI(2,ITOT)
          !CENTER(3,MCNT)=RVECI(3,ITOT)
          DIST=CENTER(1,MCNT)**2+CENTER(2,MCNT)**2+CENTER(3,MCNT)**2
          IF(DIST.LT.DMIN)DMIN=DIST
          NPAIRS(MCNT)=0
!          IF (DEBUG) THEN    !NoVec flow assumed
!           PRINT 990,MCNT,(CENTER(J,MCNT),J=1,3)
  990      FORMAT(' PCENT:',I5,3G15.6)
!          END IF
         END DO
        END DO
       END DO
C
       IF(SQRT(DMIN) .GT. 0.7D0)THEN
        MCNT=MCNT+1
        CENTER(:,MCNT)=0.0d0 
!        CENTER(1,MCNT)=0.0D0
!        CENTER(2,MCNT)=0.0D0
!        CENTER(3,MCNT)=0.0D0
        NPAIRS(MCNT)=0
        IF (DEBUG) PRINT 990,MCNT,(CENTER(J,MCNT),J=1,3)
       END IF
       IF((MCNT-1).GT.MX_CNT)THEN
        write(6,*)'COUPOT: MX_CNT MUST BE AT LEAST:',MCNT-1
        CALL STOPIT
       END IF
c
#ifdef MPI
c
c MPI: send general data and initialize pot and rhog for all workers
c
!YY senddata 101 and 102 are in sic code
       CALL SENDDATA(103)
c
#endif
C
C ZERO POT AND RHOG
C
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=KRHOG
       COULOMB(:)=0.0d0
!       DO IPTS=1,NMSH
!        COULOMB(IPTS)=0.0D0
!       END DO  
       IF (MODDEN .EQ. 1) THEN
        RHOG(:,:,:)=0.0d0
        !DO ISPN=1,NSPN
        ! DO IRHOG=1,NGRAD
        !  DO IPTS=1,NMSH
        !   RHOG(IPTS,IRHOG,ISPN)=0.0D0
        !  END DO
        ! END DO
        !END DO
       END IF
C
       CHARGE=0.0D0
       CHARGE_LONG=0.0D0
       IP(:,:)=0
!       DO LNG=1,MXLG
!        DO ICNTR=1,MCNT+1
!         IP(ICNTR,LNG)=0
!        END DO
!       END DO
       ISITE=0
       ISHELLA=0
       loop120: DO IFNCT=1,NFNCT 
       loop119: DO I_POS=1,N_POS(IFNCT) 
        JSITE_BEG=ISITE
!        ISHELLA=ISHELLA+1
        ISHELLA=sum(N_POS(1:IFNCT)) - N_POS(IFNCT)+ I_POS
        CALL OBINFO(1,RIDT(1,ISHELLA),RVECI,IST,ISHDUM)
C
C  TRANSLATE SHELLS TO ATOMS; GET GAUSS CUTOFF FOR ATOM I
C
        IATOM=ISHELLA
        GAUSS_CUTI=GAUSS_CUT(IATOM)
        loop118: DO K_SITEI=1,IST
         ISITE=ISITE+1
!         CALL UNRAVEL(IFNCT,ISHELLA,K_SITEI,RIDT(1,ISHELLA),
!     &                RVECI,LST,1)
!         IF(LST.NE.IST)THEN
!          write(6,*)'COUPOT: PROBLEM IN UNRAVEL'
!          CALL STOPIT
!         END IF
         I_LOC_BEG=0
         AI(1)=RVECI(1,K_SITEI)
         AI(2)=RVECI(2,K_SITEI)
         AI(3)=RVECI(3,K_SITEI)
         JSITE=JSITE_BEG
         JSHELLA=0

!YY. Replacing double loops to a single JSHELLA loop
! This double loops were causing -O2 and -O3 problem under intel compiler.         
!         loop116: DO JFNCT=1,NFNCT
!         loop115: DO J_POS=1,N_POS(JFNCT)
!          JSHELLA=JSHELLA+1
         loop116: DO JSHELLA=ISHELLA,sum(N_POS(1:NFNCT))
          JSHELLAtmp=JSHELLA
          do JFNCT=1,NFNCT
           JSHELLAtmp=JSHELLAtmp-N_POS(JFNCT)
           if(JSHELLAtmp <= 0) exit
          end do

          !JSHELLA=sum(N_POS(1:JFNCT)) - N_POS(JFNCT)+ J_POS
!          IF(JSHELLA.LT.ISHELLA) cycle !not exit loop115 !GOTO 115
!          write(*,*) 'opt', ISHELLA, JSHELLA, JFNCT
          CALL OBINFO(1,RIDT(1,JSHELLA),RVECJ,JST,JSHDUM)
C
C  TRANSLATE SHELLS TO ATOMS;  GET GAUSS CUT OFF  FOR ATOM J
C
          JATOM=JSHELLA
          GAUSS_CUTJ=GAUSS_CUT(JATOM)
          JJCALL=0 
          loop114: DO K_SITEJ=1,JST
           CALLEDJ=.FALSE.
           JSITE=JSITE+1
           FIRST_ALPJ=.TRUE.
           IF(JSITE.LT.ISITE) cycle !GOTO 114
           J_LOC_BEG=0
           IF(JSITE.EQ.ISITE)J_LOC_BEG=I_LOC_BEG
           TIMES_SITES=1.0D0
           IF(JSITE.GT.ISITE)TIMES_SITES=2.0D0
           AJ(1)=RVECJ(1,K_SITEJ)
           AJ(2)=RVECJ(2,K_SITEJ)
           AJ(3)=RVECJ(3,K_SITEJ)
           rarg=(AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2

           loop110: DO IALP=1,N_BARE(IFNCT) 
            ALPHAI=BFALP(IALP,IFNCT)
            JALP_BEG=1
            IF(ISITE.EQ.JSITE)JALP_BEG=IALP
            loop105: DO JALP=JALP_BEG,N_BARE(JFNCT) 
             TIMES_ALPHA=1.0D0
             IF(ISITE.EQ.JSITE.AND.JALP.GT.IALP)TIMES_ALPHA=2.0D0
             ALPHAJ=BFALP(JALP,JFNCT)
             RCALP=1.0D0/(ALPHAI+ALPHAJ)
             ARG=ALPHAI*ALPHAJ*RCALP
     &        * rarg
!     &        *((AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2)
             IF(ARG .GT.CUTEXP) cycle !GOTO 105
             LC_SHT=LC_SHT+1
             IF(ALPHAI.GT.GAUSS_CUTI.OR.ALPHAJ.GT.GAUSS_CUTJ)
     &         cycle !GOTO 105
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
!              IF(JJCALL.EQ.1.AND.K_SITEJ.NE.1)THEN
!!               CALL UNRAVEL(JFNCT,JSHELLA,1,RIDT(1,JSHELLA),
!     &                      RVECJ,MST,2)
!               IF(MST.NE.JST)THEN
!                write(6,*)'COUPOT: PROBLEM IN UNRAVEL'
!                CALL STOPIT
!               END IF
!              END IF
!              CALL UNRAVEL(JFNCT,JSHELLA,K_SITEJ,RIDT(1,JSHELLA),
!     &                     RVECJ,MST,2)
!              IF(MST.NE.JST)THEN
!               write(6,*)'COUPOT: PROBLEM IN UNRAVEL'
!               CALL STOPIT
!              END IF
              CALLEDJ=.TRUE.
             END IF
C
C END OF NEW SETUP...
C
C
C CALCULATE DENSITY MATRIX:
C
             DMAT2(:,:,:)=0.0d0
!             DO ISPN=1,NSPN
!              DO MUI=1,10
!               DO MUJ=1,10
!                DMAT(MUJ,MUI,ISPN)=0.0D0
!               END DO
!              END DO
!             END DO
!             CALL TRACER('GETTING ATOM I',ISHELLA)
!             CALL TRACER('GETTING ATOM J',JSHELLA)
             CALL GET_DENMAT_ATOMS(ISHELLA,JSHELLA)
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
!              DO LI =0,LSYMMAX(IFNCT)
!               DO MUI=IBEG(LI+1),IEND(LI+1)
!                DO IBASE=1,N_CON(LI+1,IFNCT)
!                 I_LOC=I_LOC+1
!                 J_LOC=J_LOC_BEG
!                 DO LJ =0,LSYMMAX(JFNCT)
!                  DO MUJ=IBEG(LJ+1),IEND(LJ+1)
!                   DO JBASE=1,N_CON(LJ+1,JFNCT)
!                    J_LOC=J_LOC+1
!                    IWF=0
!                    DO ISPN=1,NSPN
!                     ADD(I_LOC,J_LOC,ISPN)=0.0D0
!                     DO IWFS=1,NWFS(ISPN)
!                      IWF=IWF+1
!                       ADD(I_LOC,J_LOC,ISPN)=ADD(I_LOC,J_LOC,ISPN)
!     &                +PSI(I_LOC,IWF,1)*PSI(J_LOC,IWF,2)
!                     END DO
!                    END DO
!                   END DO
!                  END DO
!                 END DO
!                END DO
!               END DO
!              END DO
              FIRST_ALPJ=.FALSE.
 !             IF(I_LOC.GT.MAXUNSYM)THEN
 !              write(6,*)'COUPOT: MAXUNSYM MUST BE AT LEAST:',I_LOC
 !              CALL STOPIT
 !             END IF
 !             IF(J_LOC.GT.MAXUNSYM)THEN
 !              write(6,*)'COUPOT: MAXUNSYM MUST BE AT LEAST:',J_LOC
 !              CALL STOPIT
 !             END IF
             END IF
             I_LOC    =I_LOC_BEG
             DO LI =0,LSYMMAX(IFNCT)
              DO IBASE=1,N_CON(LI+1,IFNCT)
               DO MUI=IBEG(LI+1),IEND(LI+1)
                I_LOC=I_LOC+1
                J_LOC=J_LOC_BEG
                DO LJ =0,LSYMMAX(JFNCT)
                 DO JBASE=1,N_CON(LJ+1,JFNCT)
                   PROD=BFCON(IALP,IBASE,LI+1,IFNCT)
     &                 *BFCON(JALP,JBASE,LJ+1,JFNCT)
                  DO MUJ=IBEG(LJ+1),IEND(LJ+1)
                   J_LOC=J_LOC+1
                   DO ISPN=1,NSPN
                    DMAT2(MUI,MUJ,ISPN)=DMAT2(MUI,MUJ,ISPN)
     &                                +PROD*DMAT_LOCAL(I_LOC,J_LOC,ISPN)
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
              !DO MUJ=1,MUJMAX  !loop reordering suggested by opt report
               DO MUI=1,MUIMAX
                DO MUJ=1,MUJMAX
                ADDON=ADDON+DMAT2(MUI,MUJ,ISPN)*SS(MUI,MUJ)
     &               *TIMES_SITES*TIMES_ALPHA
               END DO
              END DO
             END DO
             CHARGE=CHARGE+ADDON
             DO ISPN=1,NSPN
              DO MUJ=1,MUJMAX
               DO MUI=1,MUIMAX
                DMAT2(MUI,MUJ,ISPN)=DMAT2(MUI,MUJ,ISPN)
     &                            *TIMES_SITES*TIMES_ALPHA
               END DO
              END DO
             END DO
             CHARGE_LONG=CHARGE_LONG+ADDON
C
             ND=MUIMAX !Moved to here
             MD=MUJMAX
             DO ISPN=1,NSPN
              IP(MINC,LNG)=IP(MINC,LNG)+1
              NPAIRS(MINC)=NPAIRS(MINC)+1
              LC_LNG=LC_LNG+1
              !ND=MUIMAX
              !MD=MUJMAX
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
                 DMTV(MUI,MUJ,IP(MINC,LNG),MINC,LNG)=
     &                      DMAT2(MUI,MUJ,ISPN)
               END DO
              END DO
             end do
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
c dvp+mrp 10/98: let manager do some work when all workers are busy
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
               IF (NMAST .LE. 0) THEN
                CALL CKWORKER(2,NTID)
                !GOTO 20
               !END IF
               else
                loop20: DO ICNT= MCNT+1, 1, -1
                 MWRD= (ICNT .EQ. MCNT+1)
                 IF (MWRD) THEN
                  NBLOCK= NSPN
                 ELSE
                  NBLOCK= NSPN*NMAST
                 END IF
                 DO ILNG= 1,3
                  NLOOP= (IP(ICNT,ILNG)+NBLOCK-1)/NBLOCK
                  DO ILOOP= 1,NLOOP
                   CALL CKWORKER(1,NTID)
                   IF (NCALLED .NE. NPROC) exit loop20 
                   NSET= MIN(IP(ICNT,ILNG),NBLOCK)
                   I= IP(ICNT,ILNG)-NSET+1
                   CALL POISSON2(MWRD,NSET,10,10,
     &                          ALPIV(I,ICNT,ILNG),AIV(1,I,ICNT,ILNG),
     &                          ALPJV(I,ICNT,ILNG),AJV(1,I,ICNT,ILNG),
     &                          DMTV(1,1,I,ICNT,ILNG))
                   IP(ICNT,ILNG)= IP(ICNT,ILNG)-NSET
                  END DO
                 END DO
                END DO loop20
               end if
               IF (IP(MINC,LNG) .GT. 0) THEN
                CALL PAMPOISS(1,ND,MD,MINC,LNG,NWRD,COULOMB,RHOG)
               END IF
              END IF
c
#endif
c
              CALL GTTIME(TP2)
              TIMEP=TIMEP+TP2-TP1
              IP(MINC,LNG)=0
             END IF
            end do loop105  
           end do loop110  
          end do loop114
!         end do loop115
         end do loop116
        end do loop118
       end do loop119
       end do loop120
!
! EMPTY TABLE OF PAIRS
!
       loop800: DO LNG=1,MXLG     
        loop790: DO MINC=1,MCNT+1 
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
c dvp+mrp 10/98: let manager do some work when all workers are busy
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
           IF (NMAST .LE. 0) THEN
            CALL CKWORKER(2,NTID)
           ! GOTO 30
           !END IF
           else
            loop30: DO ICNT= MCNT+1, 1, -1
             MWRD= (ICNT .EQ. MCNT+1)
             IF (MWRD) THEN
              NBLOCK= NSPN
             ELSE
              NBLOCK= NSPN*NMAST
             END IF
             DO ILNG= 1,3
              NLOOP= (IP(ICNT,ILNG)+NBLOCK-1)/NBLOCK
              DO ILOOP= 1,NLOOP
               CALL CKWORKER(1,NTID)
               IF (NCALLED .NE. NPROC) exit loop30 
               NSET= MIN(IP(ICNT,ILNG),NBLOCK)
               I= IP(ICNT,ILNG)-NSET+1
               CALL POISSON2(MWRD,NSET,10,10,
     &                       ALPIV(I,ICNT,ILNG),AIV(1,I,ICNT,ILNG),
     &                       ALPJV(I,ICNT,ILNG),AJV(1,I,ICNT,ILNG),
     &                       DMTV(1,1,I,ICNT,ILNG))
               IP(ICNT,ILNG)= IP(ICNT,ILNG)-NSET
              END DO
             END DO
            END DO loop30
           end if
!   30      continue
           IF (IP(MINC,LNG) .GT. 0) THEN
            CALL PAMPOISS(1,ND,MD,MINC,LNG,NWRD,COULOMB,RHOG)
            IP(MINC,LNG)= 0
           END IF
          END IF
c
#endif
c
          CALL GTTIME(TP2)
          TIMEP=TIMEP+TP2-TP1
         END IF
        end do loop790  
       end do loop800
c
#ifdef MPI
c
       CALL PAMPOISS(2,ND,MD,MINC,LNG,NWRD,COULOMB,RHOG)
c
#endif
c
       ALCLNG=LC_LNG
       ALCSHT=LC_SHT
       ISPFAC=2/NSPN
       SAVINGS=ALCLNG/ALCSHT
       IF (DEBUG) THEN
        write(6,*)'!!!!!!!!!!! SAVINGS:',SAVINGS,' !!!!!!!!!!!!'
        write(6,*)'NSPN:',NSPN
       END IF
       CHARGE=CHARGE*ISPFAC
       CHARGE_LONG=CHARGE_LONG*ISPFAC
       PRINT '(A,F16.6)','LONG  CHARGE IN COUPOT: ',CHARGE_LONG
       PRINT '(A,F16.6)','TOTAL CHARGE IN COUPOT: ',CHARGE
C
       RR_MAX=0.0D0
       DO IPTS=1,NMSH
        RRR=SQRT(RMSH(1,IPTS)**2+RMSH(2,IPTS)**2+RMSH(3,IPTS)**2)
        IF(RRR.GT.RR_MAX)THEN
         RR_MAX=RRR
         CHARGE_LR=RRR*COULOMB(IPTS)
        END IF
       END DO
       IF(ABS(CHARGE_LR-CHARGE) .GT. 0.1D0)THEN
        write(6,*)'************** WARNING ****************'
        write(6,*)'CHARGE VIA ASYMPTOTIC FORMULA: ',CHARGE_LR
       END IF
       CALL GTTIME(TIME2)
        write(6,'(A,F12.3,2X,I8)')'TIME FOR COULOMB POTENTIAL,NMSH:',
     %                      TIME2-TIME1,NMSH
       IF (DEBUG) THEN
        write(6,*)'TIME FOR COULOMB POTENTIAL,NMSH:',TIME2-TIME1,NMSH
        write(6,*)'POISSON DIAG: PT   :',TIMEP
        write(6,*)'POISSON DIAG: CT   :',TIME2-TIME1-TIMEP
        write(6,*)'POISSON DIAG, SETUP:',TSETUP
        write(6,*)'POISSON DIAG, POLYS:',TPOLYS
        write(6,*)'POISSON DIAG, FMTTM:',TFMTTM
        write(6,*)'POISSON DIAG, RUNIT:',TRUNIT
        write(6,*)'POISSON DIAG, ATIME:',ATIME
        write(6,*)'POISSON DIAG, WASTED:',TWASTED
        write(6,*)'POISSON DIAG, FLDNR:',
     &             FLDNR,FLDNR/(FLDNR+ASYMP+0.001D0)
        write(6,*)'POISSON DIAG, ASYMP:',
     &             ASYMP,ASYMP/(FLDNR+ASYMP+0.001D0)
        write(6,*)'POISSON DIAG, ACHRG:',ACHRG,CHARGE_LONG
        write(6,*)'POISSON DIAG, PRTME:',COSTPROJ
        write(6,997)TMWRD,TMNRM,CLWRD,CLNRM,
     &            TMWRD/(CLWRD+0.001D0),TMNRM/(CLNRM+0.001D0)
 997   FORMAT(' NEW TIMES:',6G12.4)
        DO ICNT=1,MCNT+1
         write(6,998)(CENTER(J,ICNT),J=1,3),NPAIRS(ICNT),NCALL,MCALL
 998   FORMAT(' POISSON:',3(1X,G18.6)/9X,3(1X,I18))
        END DO
       END IF

!
! DEALLOCATE ARRAYS FROM MODULE COUPDATA1 IF NO LONGER NEEDED
!
       CALL CHECK_INPUTS
       IF(.NOT.DOSJNT1)THEN
         DEALLOCATE(AIV,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING AIV'
         DEALLOCATE(AJV,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING AJV'
         DEALLOCATE(DMTV,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING DMTV'
         DEALLOCATE(ALPIV,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING ALPIV'
         DEALLOCATE(ALPJV,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING ALPJV'
         DEALLOCATE(CENTER,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING CENTER'
         DEALLOCATE(ADD,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING ADD'
         DEALLOCATE(RVECI,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING RVECI'
         DEALLOCATE(RVECJ,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING RVECJ'
         DEALLOCATE(NPAIRS,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING NPAIRS'
         DEALLOCATE(IP,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'COUPOT1:ERROR DEALLOCATING IP'
       ENDIF
       RETURN   
       END
