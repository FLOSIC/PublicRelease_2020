C UTEP Electronic Structure Lab (2020)
       SUBROUTINE SICLAGM(LSPN)
       use debug1
       use global_inputs,only : CALCTYPE1
       use mpidat1,only : IRANK,NPROC,SHM_SIZE
       use mesh1,only : wmsh,rmsh,nmsh
       use common2,only : NSPN,IGGA   !,RIDT,N_CON,LSYMMAX,N_POS,NFNCT
       use common5,only : NWFS,CONVERGENCE!,NWF,PSI
       use mixpot1,only : POT=>POTOUT   !POTIN
!      use pot_dens,only : COULOMB,RHOG
       use for_diag1
!SIC modules
       use SICMAT,only : SIC,ZPOT,ZMGGAS !,SIC_COL
       use FRM,only    : BFRM,LFRM   !,RESULTS,DEBDAX
       use MOCORB,only : NFRM        !,SLAT,ZTZL,JJJJJJ
!MGGA module
       use global_inputs,only : libxc1
       USE XTMP2A,ONLY : ISMGGA
       use scaledpzsic, only : SDSICSAV,SDSICON
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LSPN, I, IFRM, IRBOFS, IPTS, ISPN, ITASKB,
     & ITASKE, J, JWF, LFM, LPTS, MPTS, NGRAD, NMAX, NTID
       REAL(8) :: AFRM, APT1, APT2, CHGDN, CHGE, CHGUP,
     & RPTS, SMCHG, TMKIN, WMSA

       REAL(8), ALLOCATABLE :: POTFILE(:,:),POTSMALL(:,:)
       CHARACTER*12 :: ZPOTSTR,MGGASTR

       PARAMETER (NMAX=MPBLOCK)

       LOGICAL FIRST  !,EXIST
       LOGICAL LGGA   !,IUPDAT

       DIMENSION WMSA(NMAX),RPTS(3,NMAX)
       !DIMENSION SICP(NMAX)
       REAL(8) :: SICP(NMAX,MAX_OCC)
       DATA FIRST/.TRUE./
       INTEGER MFRM

       SAVE

       FIRST=.FALSE.
       MFRM=0
      
!ISMGGA statement is needed for group calculation at iter=1
       ISMGGA= ((IGGA(1) >= 3) .OR. (IGGA(2) >= 3))
       CHGE=0.0D0
       SMCHG=0.0D0
       TMKIN=0.0D0
       CHGUP=0.0D0
       CHGDN=0.0D0
       ITASKB=2
       ITASKE=2
!
!
       ALLOCATE(POTSMALL(NMAX,MAX_OCC))
       ALLOCATE(ZMGGAS(4,NMAX,MAX_OCC))
       LGGA= .TRUE.
       NGRAD=4
C
C LOOP OVER ALL POINTS
       CALL GTTIME(APT1)
       IRBOFS=(LSPN-1)*NWFS(1)
       DO LFM=1,NFRM(LSPN)
        WRITE(ZPOTSTR,'(A,I4.4)')'ZPOT',LFM+IRBOFS
        print *,206+LFM+IRBOFS,"OPENING ",ZPOTSTR
        OPEN(206+LFM+IRBOFS,FILE=ZPOTSTR,FORM='UNFORMATTED',
     &       STATUS='OLD')
       END DO
       IF(LIBXC1.OR.ISMGGA)THEN
        DO LFM=1,NFRM(LSPN)
         WRITE(MGGASTR,'(A,I4.4)')'MGGA',LFM+IRBOFS
         print *,106+LFM+IRBOFS,"OPENING ",MGGASTR
         OPEN(106+LFM+IRBOFS,FILE=MGGASTR,FORM='UNFORMATTED',
     &        STATUS='OLD')
        END DO
       ENDIF
       LPTS=0
 10    CONTINUE
       IF(LPTS+NMAX.LT.NMSH)THEN
         MPTS=NMAX
       ELSE
         MPTS=NMSH-LPTS
       END IF

! v3 Read potential from files
       DO LFM=1,NFRM(LSPN)
        DO IPTS=1,MPTS
         READ(206+LFM+IRBOFS) POTSMALL(IPTS,LFM)
        END DO
       END DO

       DO IPTS=1,MPTS
         WMSA(IPTS)=WMSH(LPTS+IPTS)
         DO LFM=1,NFRM(LSPN)
          SICP(IPTS,LFM+IRBOFS)=POTSMALL(IPTS,LFM)*WMSA(IPTS)
         END DO
!MGGA - We assume spin down is 0 (need to check this).  RHOG down is 0.
         IF(LIBXC1.OR.ISMGGA)THEN
          DO LFM=1,NFRM(LSPN)
           READ(106+LFM+IRBOFS) ZMGGAS(1:4,IPTS,LFM+IRBOFS)
           IF(SDSICON) ZMGGAS(:,IPTS,LFM+IRBOFS)=
     &                 ZMGGAS(:,IPTS,LFM+IRBOFS)*SDSICSAV(LFM,LSPN)
          END DO
         ENDIF
         DO J=1,3
           RPTS(J,IPTS)=RMSH(J,LPTS+IPTS)
         END DO
       END DO
#ifndef MPI
       CALL SICLMSLV(LSPN,MPTS,WMSA,SICP,RPTS)
#else
       IF(NPROC > 0) THEN !Prevent dead lock when NPROC=0
        CALL PAMLMSIC(1,LSPN,NFRM(LSPN),MPTS,WMSA,SICP,RPTS)
       ELSE
        CALL SICLMSLV(LSPN,MPTS,WMSA,SICP,RPTS)
       END IF
#endif


       LPTS=LPTS+MPTS
       IF (LPTS .LT. NMSH) GOTO 10
#ifdef MPI
       CALL CKWORKER(3,NTID)
       CALL GTTIME(APT2)
       WRITE(6+IRANK,*)'TIME FOR SIC MATRIX:',LFM,APT2-APT1
#ifndef MPI_3
!      CALL TRACER('CALL PAMLMSIC MODE 2',LPTS)
       IF(NPROC > 0) THEN
        CALL PAMLMSIC(2,LSPN,NFRM(LSPN),MPTS,WMSA,SICP,RPTS)
       END IF
#endif
#endif
C      PRINT*,'SMCHG:',CHGUP,CHGDN,SMCHG,TMKIN,NWFS(1)
       LGGA=.FALSE.
       WRITE(6,*)'SIC MATRIX FOR SPIN LSPN:',LSPN
       DO LFM=1,NFRM(LSPN)
        WRITE(6,1000)(SIC(JWF,LFM,LSPN),JWF=1,NFRM(LSPN))
       END DO
       PRINT *,"SICLAGM DONE",NFRM(1:2), LSPN
 1000  FORMAT(' ',15F12.5)
       IF(CONVERGENCE.AND. CALCTYPE1.NE.2) THEN
        !Keep ZPOT FILES for post conv. proc.
        DO LFM=1,NFRM(LSPN)
         CLOSE(206+LFM+IRBOFS)
         print *,206+LFM+IRBOFS," CLOSED"
        END DO
        IF(LIBXC1.OR.ISMGGA)THEN
         DO LFM=1,NFRM(LSPN)
          CLOSE(106+LFM+IRBOFS)
          print *,106+LFM+IRBOFS," CLOSED"
         END DO
        ENDIF
       ELSE !.not. convergence
        !Clean up ZPOT files for the next iter. step.
        DO LFM=1,NFRM(LSPN)
         CLOSE(206+LFM+IRBOFS,STATUS='DELETE')
         print *,206+LFM+IRBOFS," CLOSED AND DELETED"
        END DO
        IF(LIBXC1.OR.ISMGGA)THEN
         DO LFM=1,NFRM(LSPN)
          CLOSE(106+LFM+IRBOFS,STATUS='DELETE')
          print *,106+LFM+IRBOFS," CLOSED AND DELETED"
         END DO
        ENDIF
       ENDIF
       DEALLOCATE(ZMGGAS)
       DEALLOCATE(POTSMALL)
       RETURN
       END

!####################################################################
! commenting unused variables
! optimized with blas routines 11/2017: cmd
 
       SUBROUTINE SICLMSLV(LSPN,MPTS,WMSA,SICP,RPTS)
!      use mesh1,only : wmsh,rmsh,nmsh
       use common2,only : RIDT,N_CON,LSYMMAX,N_POS,NFNCT !,NSPN
       use common5,only : NWF,NWFS,PSI
!      use mixpot1,only : POTIN,POT=>POTOUT
!      use pot_dens,only : COULOMB,RHOG
       use for_diag1
!SIC modules
       use SICMAT,only : SIC,ZMGGAS  !SIC_COL
       use MOCORB,only : NFRM
!      use FRM,only    : BFRM,RESULTS,LFRM,DEBDAX
       use mpidat1,only : IRANK
!MGGA module
       use global_inputs,only : libxc1
       USE XTMP2A,ONLY : ISMGGA
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
!      PARAMETER (NMAX=MPBLOCK)
       INTEGER,PARAMETER :: NMAX=MPBLOCK
       INTEGER,INTENT(IN) :: LSPN,MPTS
       REAL(8),INTENT(IN) :: WMSA(NMAX),SICP(NMAX,MAX_OCC),
     &                       RPTS(3,NMAX)
!      REAL(8) :: WMSA, RPTS, SICP
!      DIMENSION WMSA(NMAX),RPTS(3,NMAX)
!      DIMENSION SICP(NMAX)
       INTEGER :: IPTS, ISIZE, NGRAD, I_POS, ICON, IGR,
     & IERR, IFNCT, ILOC, ISHDUM, ISHELLA, IWF, J_POS, JFM, JFN,
     & JPTS, KPTS, L_NUC, LFN, LI, LMAX1, LPV, M_NUC, MU, NDERV, NPV
     & ,NOFS
! unused variables
!      INTEGER :: I, IFRM, ISPN, ITASKB, ITASKE, J, LPTS, JWF, NMAX,
!    & NTID, 
!      REAL(8) :: SYMBOL , AFRM, APT1, APT2, CHGDN, CHGE, CHGUP, FMAT,
!    & RHI, SMCHG, TMKIN, TMP, VER, ADD1, ADD2, ADD3,
!    & ADD4, AGRAD, FACTOR
!      DIMENSION AFRM(3,MAX_OCC)
!      DIMENSION FMAT(MAX_OCC,MAX_OCC,4,2),RHI(4)
!      DIMENSION TMP(MAX_OCC,MAX_OCC),VER(MAX_OCC,MAX_OCC,4)
       INTEGER :: N_CONS,NLOC
! blas variables
       CHARACTER(1) :: transa
       REAL(8) :: psiadd,alpha,beta
       REAL(8),external :: ddot
!      INCLUDE 'commons.inc'
!      LOGICAL ICOUNT,EXIST,FIRST
!      LOGICAL EXIST,FIRST
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

       REAL*8,allocatable :: PSIG(:,:,:), PSIG1(:,:),PTS(:,:)
!      REAL(8),allocatable :: PSIG(:,:), PTS(:,:)
     &                      ,GRAD(:,:,:,:,:),RVECA(:,:)
     &                      ,PSIGLFN(:,:),G1(:,:),G4(:,:,:)
       LOGICAL,allocatable :: ICOUNT(:,:) !,IC1(:)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
!      DATA FIRST/.TRUE./

       INTEGER MXXD,KXXS,KXXO


       SAVE  !YY to save or not save?
!YY allocate tmp2 arrays
       !For LDA and GGA you only need 1 dim
       ALLOCATE(PSIG(NMAX,4,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING PSIG'
       ALLOCATE(PSIG1(NMAX,MAX_OCC),STAT=IERR) !never needs derivatives
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING PSIG1'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING NSPEED'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING GRAD'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING RVECA'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING ICOUNT'

!YY I need to add libxc case
       if(ISMGGA.OR.LIBXC1) then
        ALLOCATE(PSIGLFN(NMAX,8),STAT=IERR)
       else
        ALLOCATE(PSIGLFN(NMAX,1),STAT=IERR)
       end if
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING PSIGLFN'
       ALLOCATE(G4(NSPEED,MAXUNSYM,4),STAT=IERR) 
       ALLOCATE(G1(NSPEED,MAXUNSYM),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR ALLOCATING G1'
!      ALLOCATE(IC1(MAXUNSYM),STAT=IERR)

!      write(6+IRANK,*)'SIC_COL ALLOC',ALLOCATED(SIC_COL)
C      PRINT*,'FIRST SICLMSLV:',IRANK,NFNCT,NWFS,N_POS(1),MPTS
!      IF(LFM.EQ.0)THEN
!        PRINT*,IRANK, LFM,' FIRST LFM IS ZERO'
!        CALL STOPIT
!      END IF
!      AGRAD=0.0D0
       LGGA= .TRUE.
       CALL CHECK_INPUTS !Added for safety
       if(ISMGGA.OR.LIBXC1) then
        NGRAD=4  !meta-GGA
       else
        NGRAD=1  !LDA
       end if
       FORALL (IWF=1:NWF, IGR=1:NGRAD, IPTS=1:MPTS)
         PSIG(IPTS,IGR,IWF) = 0.0D0
       END FORALL
       FORALL (IWF=1:NWF, IPTS=1:MPTS)
         PSIG1(IPTS,IWF) = 0.0D0
       END FORALL
       ISHELLA=0
       DO 86 IFNCT=1,NFNCT
         LMAX1=LSYMMAX(IFNCT)+1
C
C FOR ALL POSITIONS OF THIS CENTER
C
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
          DO 82 J_POS=1,M_NUC
           CALL UNRAVEL2(MXXD,KXXS,KXXO,IFNCT,ISHELLA,J_POS,
     &                   RIDT(1,ISHELLA),RVECA,L_NUC,1)
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
             PTS(LPV,1)=RPTS(1,KPTS)-RVECA(1,J_POS)
             PTS(LPV,2)=RPTS(2,KPTS)-RVECA(2,J_POS)
             PTS(LPV,3)=RPTS(3,KPTS)-RVECA(3,J_POS)
            END DO
            IPTS=JPTS-1
C
C GET ORBITS AND DERIVATIVES
C
            IF(ISMGGA.OR.LIBXC1) THEN
             NDERV=1   !GGA/meta-GGA
            ELSE
             NDERV=0   !LDA
            ENDIF
            CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            IF (IUPDAT) THEN
             IPTS=JPTS-1
!         ---------------------------------------------------------
! PUT GRAD INTO MATRIX FORM (DIM:NPVxNLOC)
!
             ILOC=0
             DO78:DO LI=1,LMAX1
              N_CONS=N_CON(LI,IFNCT)
              DO MU=1,ISIZE(LI)
               DO ICON=1,N_CONS
                ILOC=ILOC+1
!                IC1(ILOC)=ICOUNT(ICON,LI)
                DO LPV=1,NPV
                 DO IGR=1,NGRAD
                  !G1(LPV,ILOC)=GRAD(LPV,1,MU,ICON,LI)
                  G4(LPV,ILOC,IGR)=GRAD(LPV,IGR,MU,ICON,LI)
                  !need deriv for meta-GGA
                 END DO 
                END DO
               END DO
              END DO
             END DO DO78
             NLOC=ILOC
!
! MULTIPLY GRADxPSI AND ADD CONTRIBUTION TO BLOCK OF PSIG (DIM:NPVxNWF)
! PSIG block starts at PSIG(IPTS+1,1)
! PSIG(MPTS,NWF)=PSIG(MPTS,NWF)+ G1(NPV,NLOC)xPSI(NLOC,NWF)
             transa='N'
             alpha=1.0d0
             beta=1.0d0 !accumulate PSIG
             DO IGR=1,NGRAD
              G1(:,:)=G4(:,:,IGR)
              PSIG1(:,:)=PSIG(:,IGR,:)
              call dgemm(transa,transa,NPV,NWF,NLOC,alpha,
     &                     G1,NSPEED, PSI(1,1,1),MAXUNSYM, beta,
     &                  PSIG1(IPTS+1,1),NMAX)
              PSIG(:,IGR,:)=PSIG1(:,:)
             END DO

!         ---------------------------------------------------------
! pre-dgemm simplification
!             DO IWF=1,NWF !changing loop order
!              DO LPV=1,NPV
!                tmp=0.0d0
!                 DO ILOC=1,NLOC
!                  tmp=tmp +G1(LPV,ILOC)*PSI(ILOC,IWF,1)
!                 END DO
!                 PSIG(IPTS+LPV,IWF)=PSIG(IPTS+LPV,IWF)+tmp
!!                 PSIG(IPTS+LPV,1,IWF)=PSIG(IPTS+LPV,1,IWF)+tmp
!              END DO
!             END DO
!         ---------------------------------------------------------
! original code
!            DO IWF=1,NWF !changing loop order
!             IPTS=JPTS-1
!             ILOC=0
!             DO 78 LI=1,LMAX1
!              DO MU=1,ISIZE(LI)
!               DO ICON=1,N_CON(LI,IFNCT)
!                ILOC=ILOC+1
!                IF (ICOUNT(ICON,LI)) THEN
!                 DO IWF=1,NWF
!                  FACTOR=PSI(ILOC,IWF,1)
!                  if(abs(FACTOR) .GT. 1.0d-10) then
!                  DO IGR=1,NGRAD
!                   DO LPV=1,NPV
!                    PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF)
!     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
!C                    AGRAD=AGRAD+ABS(GRAD(LPV,IGR,MU,ICON,LI))
!                     AGRAD=AGRAD+ABS(FACTOR)
!                   END DO
!                  END DO
!                  end if
!!                 END DO !DO IWF=1,NWF
!                END IF
!               END DO
!              END DO
!   78        CONTINUE
!            END IF
!            END DO !DO IWF=1,NWF
!         ---------------------------------------------------------
            END IF !IF (IUPDAT)
   80      CONTINUE
   82     CONTINUE
   84    CONTINUE
   86   CONTINUE
! 50     FORMAT(' ',21F10.4)
        !YY debug message to check MGGA MPI
        !print *,"ISMGGA and IRANK", ISMGGA,IRANK
C CONSTRUCT OVERLAP MATRIX OF FERMI ORBITAL
C
! ADDS ONLY USED FOR DEBUGGING. SLOWS DOWN CODE 
!        ADD1=0.0d0
!        ADD2=0.0d0
!        ADD3=0.0d0
!        ADD4=0.0d0
        IF(LSPN.NE.0)THEN
!        LFN=LFM+NWFS(1)*(LSPN-1)
!<<<<<<<<< ADD LFN LOOP HERE <<<<<<<<<<<
         NOFS=NWFS(1)*(LSPN-1) 
         DO LFN=1+NOFS,NFRM(LSPN)+NOFS
          PSIGLFN=0.0d0
          DO IPTS=1,MPTS
!          PSIGLFN(IPTS)=PSIG(IPTS,1,LFN)*WMSA(IPTS)*SICP(IPTS)
           PSIGLFN(IPTS,1)=PSIG(IPTS,1,LFN)*SICP(IPTS,LFN) !WMSA multiplication done in SICLAGM
           if(ISMGGA.OR.LIBXC1) then
            PSIGLFN(IPTS,2)=PSIG(IPTS,1,LFN)*WMSA(IPTS)
     &                     *ZMGGAS(1,IPTS,LFN)
            PSIGLFN(IPTS,3)=PSIG(IPTS,1,LFN)*WMSA(IPTS)
     &                     *ZMGGAS(2,IPTS,LFN)
            PSIGLFN(IPTS,4)=PSIG(IPTS,1,LFN)*WMSA(IPTS)
     &                     *ZMGGAS(3,IPTS,LFN)
            PSIGLFN(IPTS,5)=WMSA(IPTS)*(
     &                      PSIG(IPTS,2,LFN)*ZMGGAS(1,IPTS,LFN)+
     &                      PSIG(IPTS,3,LFN)*ZMGGAS(2,IPTS,LFN)+
     &                      PSIG(IPTS,4,LFN)*ZMGGAS(3,IPTS,LFN))
            PSIGLFN(IPTS,6)=PSIG(IPTS,2,LFN)*WMSA(IPTS)
     &                     *ZMGGAS(4,IPTS,LFN)
            PSIGLFN(IPTS,7)=PSIG(IPTS,3,LFN)*WMSA(IPTS)
     &                     *ZMGGAS(4,IPTS,LFN)
            PSIGLFN(IPTS,8)=PSIG(IPTS,4,LFN)*WMSA(IPTS)
     &                     *ZMGGAS(4,IPTS,LFN)
           end if
          END DO

          DO JFM=1,NWFS(LSPN) ! CHANGES THIS TO NWFS(LSPN) ZENGUI LOOK HERE
           JFN=JFM+NWFS(1)*(LSPN-1)

           psiadd=ddot(MPTS,PSIG(:,1,JFN),1,PSIGLFN(:,1),1) !dot product over MPTS
           if(ISMGGA .or. libxc1) THEN
            psiadd=psiadd
!GGA term 1
!>   /f$  equation1  /f$
     &      +ddot(MPTS,PSIG(:,2,JFN),1,PSIGLFN(:,2),1)
     &      +ddot(MPTS,PSIG(:,3,JFN),1,PSIGLFN(:,3),1)
     &      +ddot(MPTS,PSIG(:,4,JFN),1,PSIGLFN(:,4),1)
!GGA term 2
!>   /f$  equation2  /f$
     &      +ddot(MPTS,PSIG(:,1,JFN),1,PSIGLFN(:,5),1)
!Meta-GGA term
!>   /f$  equation3  /f$
            if(ISMGGA) then
             psiadd=psiadd
     &        +ddot(MPTS,PSIG(:,2,JFN),1,PSIGLFN(:,6),1)
     &        +ddot(MPTS,PSIG(:,3,JFN),1,PSIGLFN(:,7),1)
     &        +ddot(MPTS,PSIG(:,4,JFN),1,PSIGLFN(:,8),1)
            end if
           end if
! LB ONLY A SINGLE ROW OF SIC IS SHARED NOW
!          SIC(JFM,LFN,LSPN)=SIC(JFM,LFN,LSPN)-psiadd
! YY Why there is NOFS when SIC array has a SPIN dimension? 
           SIC(JFM,LFN-NOFS,LSPN)=SIC(JFM,LFN-NOFS,LSPN)-psiadd
!          SIC_COL(JFM)=SIC_COL(JFM)-psiadd
!LB
          END DO
!         ---------------------------------------------------------
! original code
!          DO JFM=1,NWFS(LSPN) ! CHANGES THIS TO NWFS(LSPN) ZENGUI LOOK HERE
!           LFN=LFM+NWFS(1)*(LSPN-1)
!           JFN=JFM+NWFS(1)*(LSPN-1)
!           DO IPTS=1,MPTS
!!            ADD1=ADD1+ABS(PSIG(IPTS,1,LFN))
!!            ADD2=ADD2+ABS(PSIG(IPTS,1,JFN))
!!            ADD3=ADD3+1.
!!            ADD4=ADD4+WMSA(IPTS)*SICP(IPTS)
!            SIC(JFM,LFM,LSPN)=SIC(JFM,LFM,LSPN)-
!     &      PSIG(IPTS,1,JFN)*PSIG(IPTS,1,LFN)*WMSA(IPTS)*SICP(IPTS)
!           END DO
!          END DO
!          ---------------------------------------------------------
C         DO ISPN=1,NSPN
C          JBEG= (ISPN-1)*NWFS(1)
C          DO JWF=1,NWFS(ISPN)
C           JLOC=JWF+JBEG
C           DO IPTS=1,MPTS
C            SMCHG=SMCHG+WMSA(IPTS)*PSIG(IPTS,1,JLOC)**2
C            TMKIN=TMKIN+WMSA(IPTS)*PSIG(IPTS,2,JLOC)**2
C            TMKIN=TMKIN+WMSA(IPTS)*PSIG(IPTS,3,JLOC)**2
C            TMKIN=TMKIN+WMSA(IPTS)*PSIG(IPTS,4,JLOC)**2
C           END DO
C          END DO
C         END DO
!>>>>>>>>>> LFM LOOP ENDS HERE >>>>>>>>>>>>>>
         END DO
        END IF
! 100    CONTINUE
! 1000  FORMAT(' ',15F12.5)
C      PRINT *,LFN,LFM,'SIC MAT:',ADD1,ADD2,ADD3,ADD4!(SIC(JWF,LFM,LSPN),JWF=1,5)
!      IF(LFM.EQ.0)THEN
!        PRINT*,'SIC MAT LFM IS ZERO',LFM
!        CALL STOPIT
!      END IF

!YY Deallocate tmp2
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING PSIG'
       DEALLOCATE(PSIG1,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING PSIG1'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING NSPEED'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING GRAD'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING RVECA'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING ICOUNT'
       DEALLOCATE(PSIGLFN, G1, G4, STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB3:ERROR DEALLOCATING PSIGLFN, G1'
       RETURN
       END
