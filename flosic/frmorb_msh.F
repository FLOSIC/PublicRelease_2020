C UTEP Electronic Structure Lab (2020)
       SUBROUTINE SICLAGM_MSH(LSPN,LFM)
       use debug1
       use mesh1,only : NMSH,FO_MESH
       use common2,only : RIDT,N_CON,LSYMMAX,N_POS,NFNCT,NSPN
       use common5,only : NWF,NWFS,PSI
       use mixpot1,only : POTIN,POT=>POTOUT
       use pot_dens,only : COULOMB,RHOG
       use for_diag1
!SIC modules
       use SICMAT,only : SIC
       use FRM,only    : BFRM,RESULTS,LFRM
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LSPN, LFM, I, IFRM, IPTS, ISIZE, ISPN, ITASKB,
     & ITASKE, J, JWF, LPTS, MPTS, NGRAD, NMAX, NTID, I_POS, ICON,
     & IERR, IFNCT, IGR, ILOC, ISHDUM, ISHELLA, IWF, J_POS, JFM, JFN,
     & JPTS, KPTS, L_NUC, LFN, LI, LMAX1, LPV, M_NUC, MU, NDERV, NPV
       REAL*8 :: SYMBOL , AFRM, APT1, APT2, CHGDN, CHGE, CHGUP, FMAT,
     & RHI, RPTS, SICP, SMCHG, TMKIN, TMP, VER, WMSA, ADD1, ADD2, ADD3,
     & ADD4, AGRAD, FACTOR
!       INCLUDE 'commons.inc'
       PARAMETER (NMAX=MPBLOCK)
!       LOGICAL ICOUNT,EXIST,FIRST
!       LOGICAL LGGA,IUPDAT
!       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C       COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),JJJJJJ,ZTZL
!       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
!       COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
!     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
!     &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
!       COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)
!       COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!     &   LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)
C      COMMON/HMATSIC/ OVTM(MAX_OCC,MAX_OCC,2),HMTM(MAX_OCC,MAX_OCC,2)
!YY Additional common block
!       COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
!     &  EVAL(NDH),SC1(NDH),SC2(NDH)
!       COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),ZTZL,JJJJJJ

       LOGICAL EXIST,FIRST
!       LOGICAL LGGA,IUPDAT

       DIMENSION WMSA(NMAX),RPTS(3,NMAX)
       DIMENSION SICP(NMAX)
       DIMENSION AFRM(3,MAX_OCC)
       DIMENSION FMAT(MAX_OCC,MAX_OCC,4,2),RHI(4)
       DIMENSION TMP(MAX_OCC,MAX_OCC),VER(MAX_OCC,MAX_OCC,4)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
       DATA FIRST/.TRUE./
       INTEGER MFRM
       character*10 :: fotxt
       character*20 :: filename
       logical      :: exist1
       SAVE

       write(fotxt,'(I0.3)') lfm
       filename='FOMSH'//trim(fotxt)
       inquire(file=filename,exist=exist1)
       if(exist1) then
         call tracer('skipping creating submesh',lfm)
         return
       else
#ifdef MPI
!         call tracer('initializing fo_mesh to false')
         call global_call(62)
#else
         fo_mesh(:)=.false.
#endif
       endif
!       IF(FIRST)THEN
!         OPEN(53,FILE='SICD.DAT')
!       ELSE
!         IF(LSPN.EQ.1.AND.LFM.EQ.1)REWIND(53)
!         IF(LSPN.EQ.2.AND.LFM.EQ.1)WRITE(53,*)' '
!       ENDIF
!       FIRST=.FALSE.
!       MFRM=0
!       DO ISPN=1,NSPN
!         NFRM(ISPN)=LFRM(ISPN)
!         DO IFRM=1,LFRM(ISPN)
!           MFRM=MFRM+1
!           DO J=1,3
!             AFRM(J,MFRM)=BFRM(J,IFRM,ISPN)
!           END DO
!         END DO
!       END DO
C      DO ISPN=1,NSPN
C        PRINT *, " >>> SIC: FRMORB2, NFRM(ISPN)", ISPN, NFRM(ISPN)
C      END DO
!       CHGE=0.0D0
!       SMCHG=0.0D0
!       TMKIN=0.0D0
!       CHGUP=0.0D0
!       CHGDN=0.0D0
!       ITASKB=2
!       ITASKE=2
!       DO I=1,NDH
!         DO J=1,NDH
!           OVER(J,I)=0.0D0
!         END DO
!       END DO
!       LGGA= .TRUE.
       NGRAD=4
C
C LOOP OVER ALL POINTS
C#ifdef MPI
C       CALL SENDDATA(101)
C       CALL SENDDATA(104)
C       CALL SENDDATA(105)
C       CALL SENDDATA(208)
C#endif
!       CALL TRACER('STARTING MESH LOOP')
       CALL GTTIME(APT1)
       LPTS=0
 10    CONTINUE
       IF(LPTS+NMAX.LT.NMSH)THEN
         MPTS=NMAX
       ELSE
         MPTS=NMSH-LPTS
       END IF
!      FORALL (IWF=1:NWF, IGR=1:NGRAD, IPTS=1:MPTS)
!        PSIG(IPTS,IGR,IWF) = 0.0D0
!      END FORALL
!       DO IPTS=1,MPTS
!         WMSA(IPTS)=WMSH(LPTS+IPTS)
!         SICP(IPTS)=POT(LPTS+IPTS)
!         DO J=1,3
!           RPTS(J,IPTS)=RMSH(J,LPTS+IPTS)
!           END DO
!       END DO
#ifdef MPI
       CALL PAMLM_MSH(1,LSPN,LFM,MPTS,LPTS)
!       CALL TRACER('SENT JOB OUT',LPTS)
#else
       CALL SICLM_MSH(LSPN,LFM,MPTS,LPTS)
#endif
       LPTS=LPTS+MPTS
       IF (LPTS .LT. NMSH) GOTO 10
#ifdef MPI
       CALL CKWORKER(3,NTID)
       CALL GTTIME(APT2)
       PRINT*,'TIME FOR SIC MATRIX:',LFM,APT2-APT1
       CALL PAMLM_MSH(2,LSPN,LFM,MPTS,LPTS)
#endif
       CALL WRITE_FOMESH(LFM)
!      PRINT*,'SMCHG:',CHGUP,CHGDN,SMCHG,TMKIN,NWFS(1)
!       LGGA=.FALSE.
!       PRINT*,'SIC MATRIX FOR ORBITAL LFM:',LFM
!       WRITE(53, 1000)(SIC(JWF,LFM,LSPN),JWF=1,NFRM(LSPN))
!       PRINT     1000,(SIC(JWF,LFM,LSPN),JWF=1,NFRM(LSPN))
! 1000  FORMAT(' ',15F12.5)
       RETURN
       END

!####################################################################

!       SUBROUTINE SICLM_MSH(LSPN,LFM,MPTS,WMSA,SICP,RPTS)
       SUBROUTINE SICLM_MSH(LSPN,LFM,MPTS,LPTS)
       use debug1
       use mesh1,only : RMSH,FO_MESH,FOMESH_CUTOFF
       use common2,only : RIDT,N_CON,LSYMMAX,N_POS,NFNCT,NSPN
       use common5,only : NWF,NWFS,PSI
       use mixpot1,only : POTIN,POT=>POTOUT
       use pot_dens,only : COULOMB,RHOG
       use for_diag1
!SIC modules
       use SICMAT,only : SIC
       use FRM,only    : BFRM,RESULTS,LFRM
       use mpidat1,only : IRANK
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LSPN, LFM, I, IFRM, IPTS, ISIZE, ISPN, ITASKB,
     & ITASKE, J, JWF, LPTS, MPTS, NGRAD, NMAX, NTID, I_POS, ICON,
     & IERR, IFNCT, IGR, ILOC, ISHDUM, ISHELLA, IWF, J_POS, JFM, JFN,
     & JPTS, KPTS, L_NUC, LFN, LI, LMAX1, LPV, M_NUC, MU, NDERV, NPV
       REAL*8 :: SYMBOL , AFRM, APT1, APT2, CHGDN, CHGE, CHGUP, FMAT,
     & RHI, RPTS, SICP, SMCHG, TMKIN, TMP, VER, WMSA, ADD1, ADD2, ADD3,
     & ADD4, AGRAD, FACTOR
!      INCLUDE 'commons.inc'
       PARAMETER (NMAX=MPBLOCK)
!      LOGICAL ICOUNT,EXIST,FIRST
       LOGICAL EXIST,FIRST
       LOGICAL LGGA,IUPDAT
!      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
!C     COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),JJJJJJ,ZTZL
!      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
!      COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
!    &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
!    &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
!      COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)
!      COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!    &   LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)
!C     COMMON/HMATSIC/ OVTM(MAX_OCC,MAX_OCC,2),HMTM(MAX_OCC,MAX_OCC,2)
       DIMENSION WMSA(NMAX),RPTS(3,NMAX)
       DIMENSION SICP(NMAX)
       DIMENSION AFRM(3,MAX_OCC)
       DIMENSION FMAT(MAX_OCC,MAX_OCC,4,2),RHI(4)
       DIMENSION TMP(MAX_OCC,MAX_OCC),VER(MAX_OCC,MAX_OCC,4)

!       REAL*8,allocatable :: PSIG(:,:,:), PTS(:,:)
!     &                      ,GRAD(:,:,:,:,:),RVECA(:,:)
       REAL*8 :: PTS(NSPEED,3),GRAD(NSPEED,KRHOG,6,MAX_CON,3),
     &           RVECA(3,MX_GRP)
!       LOGICAL,allocatable :: ICOUNT(:,:)
       LOGICAL :: ICOUNT(MAX_CON,3)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
       DATA FIRST/.TRUE./

       INTEGER MXXD,KXXS,KXXO
       REAL*8 :: PSIG_VALUE(MPTS)

       SAVE  !YY to save or not save?

!       CALL TRACER('FOMESH_CUTOFF',1,FOMESH_CUTOFF)
!       PRINT*,'FIRST SICLMSLV:',IRANK,NFNCT,NWFS,N_POS(1),MPTS,LFM
!       IF(LFM.EQ.0)THEN
!         PRINT*,IRANK, LFM,' FIRST LFM IS ZERO'
!         CALL STOPIT
!       END IF
       AGRAD=0.0D0
       LGGA= .TRUE.
       NGRAD=4
!       FORALL (IWF=1:NWF, IGR=1:NGRAD, IPTS=1:MPTS)
!         PSIG(IPTS,IGR,IWF) = 0.0D0
!       END FORALL
       ISHELLA=0
!       CALL TRACER('NFNCT',NFNCT)
       DO 86 IFNCT=1,NFNCT
         LMAX1=LSYMMAX(IFNCT)+1
C
C FOR ALL POSITIONS OF THIS CENTER
C
!         CALL TRACER('N_POS(IFNCT)',N_POS(IFNCT))
         DO 84 I_POS=1,N_POS(IFNCT)
          ISHELLA=ISHELLA+1
C
C GET SYMMETRY INFO
C
          CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
          IF(NWF.GT.MAX_OCC)THEN
           PRINT *,'APTSLV: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
          MXXD=2
          KXXS=LSPN     !This is not currently needed
          KXXO=0        !This should not be needed
!          CALL TRACER('M_NUC',M_NUC)
          DO 82 J_POS=1,M_NUC
           CALL UNRAVEL2(MXXD,KXXS,KXXO,IFNCT,ISHELLA,J_POS,
     &                   RIDT(1,ISHELLA),RVECA,L_NUC,1)
!           WRITE(6+IRANK,*)'AFTER UNRAVEL'
           IF(L_NUC.NE.M_NUC)THEN
            PRINT *,'APTSLV: PROBLEM IN UNRAVEL'
           CALL STOPIT
           END IF
C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
           KPTS=0
           DO 80 JPTS=1,MPTS,NSPEED
            NPV=MIN(NSPEED,MPTS-JPTS+1)
            DO LPV=1,NPV
             KPTS=KPTS+1
             PTS(LPV,1)=RMSH(1,LPTS+KPTS)-RVECA(1,J_POS)
             PTS(LPV,2)=RMSH(2,LPTS+KPTS)-RVECA(2,J_POS)
             PTS(LPV,3)=RMSH(3,LPTS+KPTS)-RVECA(3,J_POS)
            END DO
            PSIG_VALUE(:)=0.0D0
C
C GET ORBITS AND DERIVATIVES
C
            NDERV=1
            CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            IF (IUPDAT) THEN
             IPTS=JPTS-1
             ILOC=0
             DO 78 LI=1,LMAX1
              DO MU=1,ISIZE(LI)
               DO ICON=1,N_CON(LI,IFNCT)
                ILOC=ILOC+1
                IF (ICOUNT(ICON,LI)) THEN
!                 DO IWF=LFM
                  FACTOR=PSI(ILOC,LFM,1)
                  if(abs(FACTOR) .GT. 1.0d-10) then
                  DO IGR=1,NGRAD
                   DO LPV=1,NPV
!                    PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF)
!     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
!                    AGRAD=AGRAD+ABS(GRAD(LPV,IGR,MU,ICON,LI))
!                     AGRAD=AGRAD+ABS(FACTOR)
                     PSIG_VALUE(IPTS+LPV)=PSIG_VALUE(IPTS+LPV)
     &               +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
                   END DO
                  END DO
                  end if
!                 END DO
                END IF
               END DO
              END DO
   78        CONTINUE
              DO LPV=1,NPV
                IF(ABS(PSIG_VALUE(IPTS+LPV))>=FOMESH_CUTOFF) THEN
                  FO_MESH(IPTS+LPV+LPTS)=.true.
                ENDIF
              ENDDO
            END IF
   80      CONTINUE
   82     CONTINUE
   84    CONTINUE
   86   CONTINUE

       RETURN
      END
