C UTEP Electronic Structure Lab (2020)
C ----- FRMORB2() ------------
!> YY. FRMORB2 subroutine is called right before SIC apotnl.
!> with ITASKB=1,ITASKE=1 for frmorb2(0,0).
!> Also in apotnl when convergence is met. -> ITASKB=3 ITASKE=2 and -LSPN.
!> Computes PSIG similarly to DFT DENSOLD.
!> Construct FMAT
!> @param[in] LSPN spin index.
!> @param[in] LFM orbital index. Should be same as IORBX.
       SUBROUTINE FRMORB2(LSXX,LFM)
       use debug1
       use mesh1,only : wmsh,rmsh,nmsh
       use common2,only : RIDT,N_CON,LSYMMAX,N_POS,NFNCT,NSPN
       use common3,only : RMAT,NGRP
       use common5,only : NWF,NWFS,PSI
       use mixpot1,only : POT=>POTOUT
       use for_diag1,only : HAM,OVER,EVAL,SC1,FILO
       use mpidat1,only : IRANK
!SIC module
       use SICMAT,only : SIC,FMAT
       use FRM,only    : BFRM,RESULTS,LFRM,DEBDAX
       use HMATSIC,only : OVTM,HMTM
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
!MGGA moduel
       USE XTMP2A,ONLY : ISMGGA,MIXINS
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:46 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LSXX, LFM, I, I_POS, IA, IB, ICON, IERR, IF, IFM,
     & IFNCT, IFRM, IGR, ILOC, IPTS, ISHDUM, ISHELLA, ISIZE, ISP, ISPN,
     & ITASK, ITASKB, ITASKE, IWF, IX, J, J_POS, JBEG, JF, JFM, JFN,
     & JLOC, JPTS, JSPN, JWF, K, KF, KPTS, KWF, KX, L_NUC, LF, LFN, LI,
     & LMAX1, LPTS, LPV, LSPN, M_NUC, MF, MFRM, MGRP, MMSH, MPTS, MU,
     & NDERV, NF, NGRAD, NMAX, NPV
       REAL*8 :: SYMBOL , APT1, CHGDN, CHGE, CHGUP, DET, FACTOR,
     & RHI, RPTS, SICP, SMCHG, TIME1, TIME2, TMKIN, 
     & WMSA, MGGASUM1,MGGASUM2,MGGASUM3
!      INCLUDE 'commons.inc'
       PARAMETER (NMAX=MPBLOCK)
!      LOGICAL ICOUNT,EXIST
!      LOGICAL LGGA,IUPDAT
!      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C      COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),JJJJJJ,ZTZL
!      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)

       REAL*8,allocatable :: PSIG(:,:,:),PTS(:,:)
     &                      ,GRAD(:,:,:,:,:),RVECA(:,:)
       LOGICAL,allocatable :: ICOUNT(:,:)
!      COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
!    &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
!    &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
!Common for SIC
!      COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)
!      COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!    &   LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)
!      COMMON/HMATSIC/ OVTM(MAX_OCC,MAX_OCC,2),HMTM(MAX_OCC,MAX_OCC,2)
!      COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
!    &  EVAL(NDH),SC1(NDH),SC2(NDH)
!      COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),ZTZL,JJJJJJ

       LOGICAL EXIST
       LOGICAL LGGA,IUPDAT

       DIMENSION WMSA(NMAX),RPTS(3,NMAX)
       DIMENSION SICP(NMAX)
       !DIMENSION AFRM(3,MAX_OCC)
       DIMENSION RHI(4)
       !DIMENSION TMP(MAX_OCC,MAX_OCC),VER(MAX_OCC,MAX_OCC,4)
       REAL*8,ALLOCATABLE :: AFRM(:,:),TMP(:,:),VER(:,:,:)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/

       SAVE
       LSPN=ABS(LSXX)

!Allocate local arrays 
       ALLOCATE(TMP(MAX_OCC,MAX_OCC)) !Do we really need this array?
       TMP(:,:)=0.0d0
       ALLOCATE(VER(MAX_OCC,MAX_OCC,4)) !Do we really need this?
       VER(:,:,:)=0.0d0
       ALLOCATE(AFRM(3,MAX_OCC))
       AFRM(:,:)=0.0d0

!YY allocate for_diag1
!LB used mostly for MPI runs
       if(.not.allocated(HAM)) ALLOCATE(HAM(NDH,NDH))
       if(.not.allocated(OVER)) ALLOCATE(OVER(NDH,NDH))
       if(.not.allocated(EVAL)) ALLOCATE(EVAL(NDH))
       if(.not.allocated(SC1)) ALLOCATE(SC1(NDH))
       if(.not.allocated(FILO))ALLOCATE(FILO(NDH,NDH))
!YY allocate tmp2 arrays
       ALLOCATE(PSIG(NMAX,10,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR ALLOCATING PSIG'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR ALLOCATING NSPEED'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR ALLOCATING GRAD'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR ALLOCATING RVECA'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR ALLOCATING ICOUNT'

C      CALL FIXTMAT(2)
C      WRITE(INITOUT,*)'IRANK:EXECUTING FRMORB2',IRANK,LSPN,LFM
       do isp=1,mxspn
CJK     do i=1,mx_CNT
        do i=1,max_occ
         do j=1,3
          debdax(j,i,isp)=0.0d0
         end do
        end do
       end do
       MFRM=0
       DO ISPN=1,NSPN
        NFRM(ISPN)=LFRM(ISPN)
        DO IFRM=1,LFRM(ISPN)
         MFRM=MFRM+1
         DO J=1,3
          AFRM(J,MFRM)=BFRM(J,IFRM,ISPN)
         END DO
        END DO
       END DO
!      DO ISPN=1,NSPN
C       PRINT *, " >>> SIC: FRMORB2, NFRM(ISPN)", ISPN, NFRM(ISPN)
!      END DO
C      CLOSE(15)
C      call system('rm ZZZ')
C
       CHGE=0.0D0
       SMCHG=0.0D0
       TMKIN=0.0D0
       CHGUP=0.0D0
       CHGDN=0.0D0
       CALL GTTIME(APT1)
       ITASKB=1
       ITASKE=2
       IF(LSPN.NE.0)ITASKB=ITASKE
       IF(LFM.EQ.0 )ITASKE=1
       IF(LSXX.LT.0)ITASKB=ITASKE+1
       CALL GTTIME(TIME1)
       DO 1000 ITASK=ITASKB,ITASKE
        DO I=1,NDH
         DO J=1,NDH
          OVER(J,I)=0.0D0
         END DO
        END DO
        IF(ITASK.EQ.1)THEN
         MGRP=1
         MMSH=MFRM
        ELSE
         MGRP=NGRP
         MMSH=NMSH
        END IF
        LGGA= .TRUE.
        NGRAD=4
C
C LOOP OVER ALL POINTS
C
        LPTS=0
 10     CONTINUE
        IF(LPTS+NMAX.LT.MMSH)THEN
         MPTS=NMAX
        ELSE
         MPTS=MMSH-LPTS
        END IF
C
C INITIALIZE PSIG AND RHOB
C
        FORALL (IWF=1:NWF, IGR=1:NGRAD, IPTS=1:MPTS)
          PSIG(IPTS,IGR,IWF) = 0.0D0
        END FORALL

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF(ITASK.EQ.2)THEN
         DO IPTS=1,MPTS
          WMSA(IPTS)=WMSH(LPTS+IPTS)/MGRP
          IF(LSPN.EQ.0)SICP(IPTS)=1.0D0
          IF(LSPN.NE.0)SICP(IPTS)=POT(LPTS+IPTS)
          DO J=1,3
           RPTS(J,IPTS)=0.0D0
           DO K=1,3
C THa fix: change IGRP tp MGRP, because IGRP is 0 here !!
             RPTS(J,IPTS)=RPTS(J,IPTS)
     &                +RMAT(K,J,MGRP)*RMSH(K,LPTS+IPTS)
           END DO
          END DO
         END DO
        ELSE
         DO IPTS=1,MPTS
          WMSA(IPTS)=0.0D0
          DO J=1,3
           RPTS(J,IPTS)=AFRM(J,IPTS+LPTS)
          END DO
         END DO
        END IF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C FOR ALL CENTER TYPES
C
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
           WRITE(6+IRANK,*)'APTSLV: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
          DO 82 J_POS=1,M_NUC
C
C UNSYMMETRIZE
C
           IF(LSPN.EQ.0)THEN
            CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                   RVECA,L_NUC,1)
           ELSE
C           CALL LORAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
C     &                  RVECA,L_NUC,1)
            CALL UNRAVEL3(2,0,0,IFNCT,ISHELLA,J_POS,
     &                    RIDT(1,ISHELLA),
     &                   RVECA,L_NUC,1)
           END IF
           IF(L_NUC.NE.M_NUC)THEN
            WRITE(6+IRANK,*)'APTSLV: PROBLEM IN UNRAVEL'
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
                 DO IWF=1,NWF
                  FACTOR=PSI(ILOC,IWF,1)
                  if(abs(FACTOR) .GT. 1.0D-10) then
                  !print *,"FACTOR",ILOC,IWF,LPV,FACTOR
                  DO IGR=1,NGRAD
                   DO LPV=1,NPV
C                    WRITE(INITOUT,*)'PSIG',PSIG(IPTS+LPV,IGR,IWF)
                    PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF)
     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
C                    WRITE(INITOUT,*)'FACTOR',FACTOR
C                    WRITE(INITOUT,*)'GRAD',GRAD(LPV,IGR,MU,ICON,LI)
                   END DO
                  END DO
                  end if
                 END DO
                END IF
               END DO
              END DO
   78        CONTINUE
            END IF
   80      CONTINUE
   82     CONTINUE
   84    CONTINUE
   86   CONTINUE
C       CALL STOPIT
!############################################
!#               ITASK == 1                 #
!############################################
C CONSTRUCT FMAT:
        IF(ITASK.EQ.1)THEN
C        PRINT*,'PSI:'
C        WRITE(INITOUT,*)'PSI:IRANK',IRANK
!         DO IPTS=1,MPTS
C         PRINT 50,(RPTS(J,IPTS),J=1,3)
C         WRITE(INITOUT,50)(RPTS(J,IPTS),J=1,3)
!         END DO
!         DO IPTS=1,MPTS
C         PRINT 50,(PSIG(IPTS,1,JWF),JWF=1,5)
C         WRITE(INITOUT,50)(PSIG(IPTS,1,JWF),JWF=1,5)
!         END DO
         DO IPTS=1,MPTS
          JFM=IPTS+LPTS
                            JSPN=1
          IF(JFM.GT.NFRM(1))JSPN=2
          JFM=JFM-(JSPN-1)*NFRM(1)
C         JFMX=MAX(JFM,JFMX)
          JBEG= (JSPN-1)*NWFS(1)
          DO K=1,4
           RHI(K)=0.0D0
          END DO
          DO JWF=1,NWFS(JSPN)
           JLOC=JWF+JBEG
           DO K=1,4
            RHI(K)=RHI(K)+PSIG(IPTS,1,JLOC)
     &                   *PSIG(IPTS,K,JLOC)
           END DO
          END DO
C         WRITE(INITOUT,*)'IPTS,RHI:',IPTS,RHI(1),JSPN
          DO K=1,4
           DO JWF=1,NWFS(JSPN)
            JLOC=JWF+JBEG
            FMAT(JWF,JFM,K,JSPN)=PSIG(IPTS,K,JLOC)/SQRT(RHI(1))
            IF(K.GE.2)THEN
!> \f$ \frac {d} {dx} \frac {\Psi_i} {\sqrt{\rho}} \f$
             FMAT(JWF,JFM,K,JSPN)=FMAT(JWF,JFM,K,JSPN)-
     &       FMAT(JWF,JFM,1,JSPN)*RHI(K)/RHI(1) !NOTE THAT RHI(K)=0.5 (Grad RHO)
            END IF
           END DO
          END DO

         END DO
C        DO JSPN=1,NSPN
C         !WRITE(INITOUT,*)'JFMX:',NWFS(JSPN),JSPN,LSPN
C         print *,'JFMX:',NWFS(JSPN),JSPN,LSPN
C         DO IFM=1,NWFS(JSPN)
C          !WRITE(INITOUT,50)(FMAT(JWF,IFM,1,JSPN),JWF=1,NWFS(JSPN))
C          print 50,(FMAT(JWF,IFM,1,JSPN),JWF=1,NWFS(JSPN))
C         END DO
C        END DO
        END IF
 50     FORMAT(' ',21F10.4)
!#############################################
!#              ITASK == 2                   #
!#############################################
C CONSTRUCT OVERLAP MATRIX OF FERMI ORBITAL
        IF(ITASK.EQ.2)THEN
         IF(LSPN.NE.0)THEN





          DO JFM=1,NWFS(LSPN) ! CHANGES THIS TO NWFS(LSPN) ZENGUI LOOK HERE
           LFN=LFM+NWFS(1)*(LSPN-1)
           JFN=JFM+NWFS(1)*(LSPN-1)
C          DO IPTS=1,MPTS
C           SIC(JFM,LFM,LSPN)=SIC(JFM,LFM,LSPN)-
C    &      PSIG(IPTS,1,JFM)*PSIG(IPTS,1,LFM)*WMSA(IPTS)*SICP(IPTS)
C          END DO
           DO IPTS=1,MPTS
            SIC(JFM,LFM,LSPN)=SIC(JFM,LFM,LSPN)-
     &      PSIG(IPTS,1,JFN)*PSIG(IPTS,1,LFN)*WMSA(IPTS)*SICP(IPTS)
           END DO
!YY If MGGA add the extra terms here
           print *,"frmorb2, ismgga", ISMGGA
           IF(ISMGGA) THEN
            MGGASUM1=0.0d0
            MGGASUM2=0.0d0
            MGGASUM3=0.0d0
            DO IPTS=1,MPTS
             MGGASUM1=MGGASUM1
     &      +PSIG(IPTS,2,JFN)*PSIG(IPTS,1,LFN)*WMSA(IPTS)*MIXINS(1,IPTS)
     &      +PSIG(IPTS,3,JFN)*PSIG(IPTS,1,LFN)*WMSA(IPTS)*MIXINS(2,IPTS)
     &      +PSIG(IPTS,4,JFN)*PSIG(IPTS,1,LFN)*WMSA(IPTS)*MIXINS(3,IPTS)
!
             MGGASUM2=MGGASUM2
     &      +PSIG(IPTS,1,JFN)*PSIG(IPTS,2,LFN)*WMSA(IPTS)*MIXINS(1,IPTS)
     &      +PSIG(IPTS,1,JFN)*PSIG(IPTS,3,LFN)*WMSA(IPTS)*MIXINS(2,IPTS)
     &      +PSIG(IPTS,1,JFN)*PSIG(IPTS,4,LFN)*WMSA(IPTS)*MIXINS(3,IPTS)
!
             MGGASUM3=MGGASUM3
     &      +PSIG(IPTS,2,JFN)*PSIG(IPTS,2,LFN)*WMSA(IPTS)*MIXINS(4,IPTS)
     &      +PSIG(IPTS,3,JFN)*PSIG(IPTS,3,LFN)*WMSA(IPTS)*MIXINS(4,IPTS)
     &      +PSIG(IPTS,4,JFN)*PSIG(IPTS,4,LFN)*WMSA(IPTS)*MIXINS(4,IPTS)
            END DO
          SIC(JFM,LFM,LSPN)=SIC(JFM,LFM,LSPN)-MGGASUM1-MGGASUM2-MGGASUM3
           END IF
!ISMGGA ends here
          END DO

         END IF
         DO IPTS=1,MPTS
          DO ISPN=1,NSPN
           RHI=0.0D0
           JBEG= (ISPN-1)*NWFS(1)
           DO JWF=1,NWFS(ISPN)
            JLOC=JWF+JBEG
            RHI=RHI+PSIG(IPTS,1,JLOC)**2
           END DO
           IF(ISPN.EQ.1   )CHGUP=CHGUP+RHI(1)*WMSA(IPTS)
           IF(ISPN.EQ.NSPN)CHGDN=CHGDN+RHI(1)*WMSA(IPTS)
          END DO
         END DO
C
C UPDATING RHOG, START WITH DENSITY
C
         DO ISPN=1,NSPN
          JBEG= (ISPN-1)*NWFS(1)
          DO JWF=1,NWFS(ISPN)
           JLOC=JWF+JBEG
           DO IPTS=1,MPTS
            SMCHG=SMCHG+WMSA(IPTS)*PSIG(IPTS,1,JLOC)**2
            TMKIN=TMKIN+WMSA(IPTS)*PSIG(IPTS,2,JLOC)**2
            TMKIN=TMKIN+WMSA(IPTS)*PSIG(IPTS,3,JLOC)**2
            TMKIN=TMKIN+WMSA(IPTS)*PSIG(IPTS,4,JLOC)**2
           END DO
          END DO
         END DO
        END IF
!###################################
!#          ITASK==2 ENDS HERE     #
!###################################
 100    CONTINUE
        LPTS=LPTS+MPTS
        IF (LPTS .LT. MMSH) GOTO 10
        IF(LSPN.NE.0.AND.LFM.GT.0)THEN
         DO JFM=1,NFRM(LSPN)
          HMTM(JFM,LFM,LSPN)=SIC(JFM,LFM,LSPN)
         END DO
C        PRINT*,'SMCHG:',CHGUP,CHGDN,SMCHG,TMKIN,NWFS(1)
C        WRITE(INITOUT,*)'HMTM:',LFM,LSPN
C        WRITE(INITOUT,*)'HMTM:',(HMTM(JFM,LFM,LSPN),JFM=1,NFRM(LSPN))
C        WRITE(INITOUT,*)'SMCHG:',CHGUP,CHGDN,SMCHG,TMKIN,NWFS(1)
        END IF
        LGGA=.FALSE.
        CALL GTTIME(TIME2)
        WRITE(6+IRANK,*)'TIME in 1000 LOOP:',TIME2-TIME1
 1000   CONTINUE
C       IF(TIME2.EQ.TIME2)CALL STOPIT
!       IF(LFM==1)THEN
!        WRITE(98,*)NSPN,NFRM(1),NWFS(1)
!        DO K=1,4
!         WRITE(98,*)'K:',K
!         DO ISPN=1,NSPN
!          DO IFM=1,NFRM(ISPN)
!           WRITE(98,50) (FMAT(JWF,IFM,K,ISPN),JWF=1,NWFS(ISPN))
!          END DO
!         END DO
!        END DO
!       ENDIF
C CALCULATE THE OVERLAP MATRIX OF THE FERMI ORBITALS
C       PRINT*,'(Xj,Yj,Zj) is the position of the jth FO'
C       WRITE(INITOUT,*)'(Xj,Yj,Zj) is the position of the jth FO'
        DO 1200 ISPN=1,NSPN
        DO 1200 K=1,4
C        PRINT*,'ISPN,K:',ISPN,K
C        IF(K.EQ.1)PRINT*,'<Fj|Fi>:'
C        IF(K.EQ.2)PRINT*,'<dFj/dXj|Fi>:'
C        IF(K.EQ.3)PRINT*,'<dFj/dYj|Fi>:'
C        IF(K.EQ.4)PRINT*,'<dFj/dZj|Fi>:'
C        IF(K.EQ.1)WRITE(INITOUT,*)'<Fj|Fi>:'
C        IF(K.EQ.2)WRITE(INITOUT,*)'<dFj/dXj|Fi>:'
C        IF(K.EQ.3)WRITE(INITOUT,*)'<dFj/dYj|Fi>:'
C        IF(K.EQ.4)WRITE(INITOUT,*)'<dFj/dZj|Fi>:'
         DO IFM=1,NFRM(ISPN)
          DO JFM=1,NFRM(ISPN)
           OVER(JFM,IFM)=0.0D0
           DO KWF=1,NWFS(ISPN)
            OVER(JFM,IFM)=OVER(JFM,IFM)+
     &      FMAT(KWF,IFM,1,ISPN)*FMAT(KWF,JFM,K,ISPN)
            !This is <phi|dphi/dam>? 
!> No, it is \f$ dSnm/dam = \langle F_n | dFm/dam  \rangle \f$ 
           END DO
          END DO
          IF(K.GT.1)THEN
           DO JFM=1,NFRM(ISPN)
            OVER(JFM,IFM)=2.0D0*OVER(JFM,IFM)
           END DO
          END IF
C         WRITE(INITOUT,50) (OVER(JFM,IFM),JFM=1,NFRM(ISPN))
         END DO
         IF(ISPN.EQ.LSPN)THEN
          DO IF=1,NFRM(ISPN)
           DO JF=1,NFRM(ISPN)
            VER(JF,IF,K)=OVER(JF,IF)
           END DO
          END DO
         END IF
         IF(K.EQ.1)THEN
          DO I=1,NFRM(ISPN)
           DO J=1,NFRM(ISPN)
            HAM(J,I)=0.0D0
           END DO
           HAM(I,I)=1.0D0
          END DO
          CALL LOWSIC(ISPN,NDH,NFRM(ISPN),OVER,HAM,FILO,EVAL,SC1)
          !YY You get FILO arrays back from LOWSIC
C THa:        no idea what this call does
          CALL SWAPFGI(ISPN,NFRM(ISPN),0)
          WRITE(6+IRANK,*)'BACK FROM LOWSIC'
          INQUIRE(FILE='DETONLY', EXIST=EXIST)
          IF(EXIST.AND.ISPN.EQ.NSPN)CALL STOPIT
C>  FIND EXPRESSION FOR LO IN TERMS OF CO:
          DO IWF=1,NFRM(ISPN)
           DO IF=1,NFRM(ISPN)
            SLAT(IF,IWF,ISPN)=0.0D0
           END DO
C          PRINT 50,(HAM(JF,IWF),JF=1,NFRM(ISPN))
C          PRINT 50,(FMAT(IWF,JF,1,ISPN),JF=1,NFRM(ISPN))
          END DO
          DO JF=1,NFRM(ISPN)
           DO IWF=1,NWFS(ISPN)
            DO IF=1,NFRM(ISPN)
             SLAT(IF,IWF,ISPN)=SLAT(IF,IWF,ISPN)+
     &        HAM(IF,JF)*FMAT(IWF,JF,1,ISPN)
            END DO
           END DO
          END DO
C CHECK IF SLAT ARE ORTHONORMAL
C         DO IWF=1,NWFS(ISPN)
C          DO JWF=1,NWFS(ISPN)
C           DOT=0.0D0
C           DO IF=1,NFRM(ISPN)
C            DOT=DOT+SLAT(IF,IWF,ISPN)*SLAT(IF,JWF,ISPN)
C           END DO
C           PRINT*, IWF,JWF,DOT
C          END DO
C         END DO
C         READ*,DUMMY
C
C
          DET=1.0d0
          DO I=1,NFRM(ISPN)
           DET=DET/EVAL(I)**2
C          PRINT*,1.0d0/EVAL(I)**2,DET
C          WRITE(INITOUT,*)1.0d0/EVAL(I)**2,DET
           EVAL(I)=1.0D0/EVAL(I)**2
           IF(EVAL(I).LT.0.01)THEN
            WRITE(6+IRANK,*)' BOTHERSOME STATE:'
            WRITE(6+IRANK,*)(HAM(JF,I),JF=1,NFRM(ISPN))
C           WRITE(INITOUT,*)' BOTHERSOME STATE:'
C           WRITE(INITOUT,*)(HAM(JF,I),JF=1,NFRM(ISPN))
           END IF
          END DO
C         PRINT 777,EVAL(1),EVAL(NFRM(ISPN)),DET
C         WRITE(INITOUT,777)EVAL(1),EVAL(NFRM(ISPN)),DET
 777      FORMAT(3G15.6,' =DETERMINANT')
         ELSE
C Calculate dQ_i/dAxj:  Derivative terms
          DO JF=1,NFRM(ISPN)
           DO KF=1,NFRM(ISPN)
            DO MF=1,NFRM(ISPN)
             SC1(MF)=OVER(MF,JF)+OVER(MF,KF) !or transpose?
            END DO
            DO IA=1,NFRM(ISPN)
             DO MF=1,NFRM(ISPN)
              !YY This is T_aj T_ak
              TMP(IA,MF)=TMP(IA,MF)+FILO(JF,IA)*FILO(KF,IA)*
     &                 SC1(MF)
             END DO
            END DO
           END DO
          END DO
C         PRINT*,'dQ_i/dA_mx',K
C         WRITE(INITOUT,*)'dQ_i/dA_mx',K
C         DO IA=1,MF
C          PRINT 50, (OVER(IA,MF),MF=1,NFRM(ISPN))
C          WRITE(INITOUT,50) (OVER(IA,MF),MF=1,NFRM(ISPN))
C         END DO
C>  YY I may need this.
C>  \f$dQ_i/dA_mx\f$
C>   OVER(IA,MF)
C>
         END IF
 1200   CONTINUE
C       PRINT *,'IRANK,LSPN',IRANK,LSPN
C       IF(LSPN.NE.0.AND.LFM.EQ.NFRM(LSPN))THEN
C        PRINT*,'SIC MATRIX:',LSPN,LFM
C        WRITE(INITOUT,*)'SIC MATRIX:',LSPN,LFM
C Return in cas FORMORB(0,X) is called
        IF(LSPN.EQ.0) then
         !Deallocate and return
         GOTO 2400
         RETURN
        ENDIF

C       DO IFM=1,NFRM(LSPN)
C        PRINT *,          (SIC(IFM,JFM,LSPN),JFM=1,NFRM(LSPN))
C       END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 2323  format('23SP:',23,12F10.3)
C THa:       again no idea what that call did ?
C       DO iii=1, NSPN
         CALL SWAPFGI(lspn,NFRM(lspn),1)
C       END DO

C       IF(LSPN.EQ.1)OPEN(47,FILE='FILO1',FORM='UNFORMATTED')
C       IF(LSPN.EQ.2)OPEN(47,FILE='FILO2',FORM='UNFORMATTED')
C       REWIND(47)
C       READ(47)((FILO(JFM,IFM),JFM=1,NFRM(LSPN)),IFM=1,NFRM(LSPN))
C       READ(47) (EVAL(JFM)    ,JFM=1,NFRM(LSPN))
C       CLOSE(47)

        do ia=1,nfrm(lspn)
         eval(ia)=1.0d0/eval(ia)**2
        end do
c       print 2323,lspn,nfrm(lspn),(eval(ia),ia=1,nfrm(lspn))
        do ia=1,nfrm(lspn)
         do ib=1,nfrm(lspn)
          HAM(ib,ia)=0.0d0
          do kf=1,nfrm(lspn)
           do lf=1,nfrm(lspn)
!>  Equation: 
!>  \f$ ( T_ak T_bl - Tal Tbk ) \f$
!>  \f$  \langle \phi_l | V_k | \phi_k \rangle \f$
            HAM(ib,ia)=HAM(ib,ia)+
     &      (FILO(kf,ia)*FILO(lf,ib)-FILO(lf,ia)*FILO(kf,ib))*
     &      sic(kf,lf,lspn)
           end do
          end do
!>  \f$ \sqrt{Q_a*Q_b} \f$
          HAM(ib,ia)=HAM(ib,ia)/(sqrt(EVAL(ib))*sqrt(EVAL(ia)))
c         HAM(ib,ia)=HAM(ib,ia)*(sqrt(EVAL(ib))-sqrt(EVAL(ia)))
c         HAM(ib,ia)=HAM(ib,ia)/(sqrt(EVAL(ib))+sqrt(EVAL(ia)))
         end do
c        print 50,(HAM(ib,ia),ib=1,nfrm(lspn))
        end do
        do ia=1,nfrm(lspn)
         do nf=1,nfrm(lspn)
          SC1(nf)=0.0d0
         end do
         do ib=1,nfrm(lspn)
          do nf=1,nfrm(lspn)
!>  \f$ HAM* T_{b,n}\f$
           SC1(nf)=HAM(ib,ia)*FILO(nf,ib)+SC1(nf)
          end do
         end do
         do kx=1,3
          do if=1,nfrm(lspn)
           do mf=1,nfrm(lspn)
            OVER(if,mf)=0.0D0
            do kwf=1,nwfs(lspn)
!>  \f$ dS_{n,m}/dam  \f$
             OVER(if,mf)=OVER(if,mf)+
     &       fmat(kwf,if,1,lspn)*fmat(kwf,mf,kx+1,lspn)
            end do
           end do
          end do
          do jf=1,nfrm(lspn)
           OVER(jf,if)=2.0d0*OVER(jf,if)
          end do
          do mf=1,nfrm(lspn)
           do nf=1,nfrm(lspn)
            debdax(kx,mf,lspn)=debdax(kx,mf,lspn)-FILO(mf,ia)*SC1(nf)*
     &        OVER(nf,mf) ! this has to be <Fn|dFm/dAmx> <- this is over.
           end do
          end do
         end do
        end do
C       open(39,file='fforce.dat')
C       rewind(39)
c       do if=1,nfrm(lspn)
c        write(39,*)(debdax(ix,if),ix=1,3)
c       end do
        write(6+irank,*)'right hand side of equation 51'
        do ia=1,nfrm(lspn)
         do ib=1,nfrm(lspn)
          HAM(ib,ia)=HAM(ib,ia)*(sqrt(EVAL(ib))-sqrt(EVAL(ia)))
          HAM(ib,ia)=HAM(ib,ia)/(sqrt(EVAL(ib))+sqrt(EVAL(ia)))
!> \f$ (\sqrt{Qb}-\sqrt{Qa})/(\sqrt{Qb}+\sqrt{Qa})  \f$
         end do
        end do
        do 600 kx=1,3
c construct OVER (left hand side of equation 51)
c        print*,'kx, OVER:',kx
         DO IF=1,NFRM(LSPN)
          DO MF=1,NFRM(LSPN)
           OVER(MF,IF)=0.0D0
           DO KWF=1,NWFS(LSPN)
            OVER(MF,IF)=OVER(MF,IF)+
     &      FMAT(KWF,IF,1,LSPN)*FMAT(KWF,MF,KX+1,LSPN)
           END DO
          END DO
C         IF(KX.GT.1)THEN  !KX=1?
           DO JF=1,NFRM(LSPN)
            OVER(JF,IF)=2.0D0*OVER(JF,IF)
           END DO
C         END IF
c         print 50,(OVER(jf,if),jf=1,nfrm(lspn))
         END DO
         do ia=1,nfrm(lspn)
          do ib=1,nfrm(lspn)
           do mf=1,nfrm(lspn)
            SC1(mf)=0.0d0
            do if=1,nfrm(lspn)
             SC1(mf)=SC1(mf)-
     &           (FILO(if,ib)*FILO(mf,ia)+FILO(mf,ib)*FILO(if,ia))*
     &           OVER(mf,if)*0.5d0
c            print*,if,mf,SC1(mf)
            end do
           end do
           do mf=1,nfrm(lspn)
            debdax(kx,mf,lspn)=debdax(kx,mf,lspn)
     &                  -0.5d0*SC1(mf)*HAM(ib,ia) !Last change
           end do
c          print*,mf,(debdax(ix,mf,lspn),ix=1,3)
          end do
         end do
 600    continue
#ifndef GROUP
        open(40,file='fforce.prt')
        rewind(40)
C       write(initout, *) 'FROM FRMORB2',lfm,lspn
        do if=1,nfrm(lspn)
C        write(initout,*)if,(debdax(ix,if,lspn),ix=1,3)
         write(40,39)(debdax(ix,if,lspn),ix=1,3)
        end do
 39     format(' ',3F20.12)
        close(40)
        call system('cat fforce.prt >>fforce.out')
#endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       END IF

2400   CONTINUE

!YY Deallocate tmp2
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR DEALLOCATING PSIG'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR DEALLOCATING PTS'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR DEALLOCATING GRAD'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR DEALLOCATING RVECA'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRMORB2:ERROR DEALLOCATING ICOUNT'

!      DEALLOCATE(HAM)
!      DEALLOCATE(OVER)
!      DEALLOCATE(EVAL)
!      DEALLOCATE(SC1)
!      DEALLOCATE(FILO)

       DEALLOCATE(AFRM)
       DEALLOCATE(VER)
       DEALLOCATE(TMP)

       RETURN
       END
