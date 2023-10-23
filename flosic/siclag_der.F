C UTEP Electronic Structure Lab (2020)
       SUBROUTINE SICLAG_DER
!
!  KAJ 1-8-2018 This computes <del(phi_i)/del(R_nu)|VSIC_i|phi_i>
!        for each atomic site, nu; compute for the orbital LFM
!
       use mesh1,only : wmsh,rmsh,nmsh
       use common2,only : RIDT,N_CON,LSYMMAX,N_POS,NFNCT,NSPN
       use common5,only : NWF,NWFS,PSI
       use mixpot1,only : POTIN,POT=>POTOUT
       use for_diag1
!SIC modules
       use SICMAT,only : SIC,DERSIC
       use FRM,only    : BFRM,RESULTS,LFRM,DEBDAX
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
! Converion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

       INCLUDE  'PARAMA2'  
       INTEGER :: LSPN, LFM, I, IID,IFRM, IPTS, ISIZE, ISPN, ITASKB,
     & ITASKE, J, JWF, LPTS, MPTS, NGRAD, NMAX, NTID, I_POS, ICON,IX,
     & IERR, IFNCT, IGR, ILOC, ISHDUM, ISHELLA, IWF, J_POS, JFM, JFN,
     & JPTS, KPTS, L_NUC, LFN, LI, LMAX1, LPV, M_NUC, MU, NDERV, NPV,
     & NTOT,NMSH1
       REAL*8 :: SYMBOL , APT1, APT2, CHGDN, CHGE, CHGUP, FMAT,
     & RHI, RPTS, SMCHG, TMKIN, WMSA, ADD1, ADD2, 
     & ADD3,ADD4, AGRAD, FACTOR,over_test,testsic,
     & TIMET1,TIMET2
       !Additional vars
       INTEGER :: IRBOFS
       REAL(8),ALLOCATABLE :: POTSMALL(:,:)
       CHARACTER*12 :: ZPOTSTR
       PARAMETER (NMAX=MPBLOCK)

       LOGICAL EXIST,FIRST,IUPDAT

       DIMENSION WMSA(NMAX),RPTS(3,NMAX)
       REAL(8) :: SICP(NMAX,MAX_OCC) !New size 
C
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
       DATA FIRST/.TRUE./
       INTEGER MFRM

       SAVE

       CALL GTTIME(TIMET1)
c
! allocate array to store V_SIC temporarily
       ALLOCATE(POTSMALL(NMAX,MAX_OCC))

       NGRAD=4
c      IF(LSPN.EQ.1) dersic = 0.0d0
       dersic = 0.0d0
       CALL GTTIME(APT1)
c
c  open files containing SIC potentials
c
       DO LSPN = 1,NSPN
       IRBOFS=(LSPN-1)*NWFS(1)
       DO LFM=1,NFRM(LSPN)
        WRITE(ZPOTSTR,'(A,I4.4)')'ZPOT',LFM+IRBOFS
!       print *,206+LFM+IRBOFS,"OPENING ",ZPOTSTR
        OPEN(206+LFM+IRBOFS,FILE=ZPOTSTR,FORM='UNFORMATTED',
     &       STATUS='OLD')
       END DO
       END DO
C
C LOOP OVER ALL POINTS
C
       LPTS=0
 10    CONTINUE
         IF(LPTS+NMAX.LT.NMSH)THEN
           MPTS=NMAX
         ELSE
           MPTS=NMSH-LPTS
         END IF
! Read potential from files
         DO LSPN = 1,NSPN
         IRBOFS=(LSPN-1)*NWFS(1)
           DO LFM=1,NFRM(LSPN)
            LFN = LFM + (LSPN-1)*NFRM(1)
            DO IPTS=1,MPTS
             READ(206+LFM+IRBOFS) POTSMALL(IPTS,LFN)
            END DO
           END DO
         end do

         DO IPTS=1,MPTS
           WMSA(IPTS)=WMSH(LPTS+IPTS)
           DO LSPN = 1,NSPN
           DO LFM=1,NFRM(LSPN)
            LFN = LFM + (LSPN-1)*NFRM(1)
            SICP(IPTS,LFN)=POTSMALL(IPTS,LFN)
           END DO
           END DO
           DO J=1,3
             RPTS(J,IPTS)=RMSH(J,LPTS+IPTS)
           END DO
         END DO
#ifndef MPI
c
         CALL SICLM_DERSLV(MPTS,WMSA,SICP,RPTS)
#else
c
c  NTOT and NSPN are now dummy indices in PAMLMSIC_DER
c  and PASLMSIC_DER
c
         NTOT = NFRM(1) + NFRM(2)
         CALL PAMLMSIC_DER(1,NSPN,NTOT,MPTS,WMSA,SICP,RPTS)
#endif
         LPTS=LPTS+MPTS
         IF (LPTS .LT. NMSH) GOTO 10
#ifdef MPI
         CALL CKWORKER(3,NTID)
         CALL GTTIME(APT2)
         CALL PAMLMSIC_DER(2,NSPN,NTOT,MPTS,WMSA,SICP,RPTS)
#endif
       print *, '--------------------------'
       print *, '        end of loop       '
       print *, '--------------------------'
       PRINT*,'DERSIC MATRIX FOR SPIN:',LSPN
       DO LSPN = 1,NSPN
       DO LFM=1,LFRM(LSPN)
         LFN=LFM+LFRM(1)*(LSPN-1)
         PRINT*,'DERSIC MATRIX FOR ORBITAL:',LFN
         do iid = 1,5
           PRINT 1000,(DERSIC(IX,LFN,IID),IX=1,3)
         end do
       END DO
       end do
       
       print *, 'NSPN in SICLAG_DER', nspn
       if(nspn.eq.1) then
        do iid = 1, MX_CNT
        do LFN = 1,LFRM(1) + LFRM(2)
           do ix = 1,3
             DERSIC(IX,LFN,IID) = 2.0d0*DERSIC(IX,LFN,IID)
           end do
        end do
        end do
       end if        
1000   FORMAT(' ',3F12.5)
       DO LSPN = 1,NSPN
       IRBOFS = (LSPN -1)*NWFS(1)
       DO LFM=1,LFRM(LSPN)
c        CLOSE(206+LFM+IRBOFS,STATUS='DELETE')
         CLOSE(206+LFM+IRBOFS)
       END DO
       END DO
       DEALLOCATE(POTSMALL)

       CALL GTTIME(TIMET2)
       WRITE(6,*)'TIME FOR SICLAG_DER/SPIN:',TIMET2-TIMET1,LSPN

       RETURN
       END

!####################################################################

       SUBROUTINE SICLM_DERSLV(MPTS,WMSA,SICP,RPTS)
       use common2,only : RIDT,N_CON,LSYMMAX,N_POS,NFNCT,NSPN
       use common3,only : RMAT, NGRP
       use common5,only : N_OCC,NWF,NWFS,PSI
       use common8,only : REP, N_REP, NDMREP
       use pot_dens,only : COULOMB,RHOG
       use for_diag1
!SIC modules
       use SICMAT,only : SIC,DERSIC
       use FRM,only : LFRM 
       use mpidat1,only : IRANK
! Converion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

       INCLUDE  'PARAMA2'  
       INTEGER :: LSPN, LFM, I, IFRM, IPT, ISIZE, ISPN, ITASKB,
     & ITASKE, J, JWF, LPTS, MPTS, NGRAD, NMAX, NTID, I_POS, ICON,
     & IERR, IFNCT, IGR, ILOC, ISHDUM, IID, IWF, J_POS, JFM, JFN,
     & MDUM,MDX,IGRP,IOCC,IROW,JROW,JWF_BEG,K_REP,KPTS,ix,LFN,LI,
     & LMAX1,LPV,m_nuc,mu,nderv,npv,jpts,l_nuc,jbeg,jend
       INTEGER :: ISPFAC
       REAL*8 :: SYMBOL , APT1, APT2, CHGDN, CHGE, CHGUP, FMAT,
     & RHI, RPTS, SICP, SMCHG, TMKIN, WMSA, ADD1, ADD2, ADD3,
     & ADD4, AGRAD, FACTOR, RR, X, Y, Z, SPFAC
       PARAMETER (NMAX=MPBLOCK)
       LOGICAL EXIST,FIRST,IUPDAT
       DIMENSION WMSA(NMAX),RPTS(3,NMAX)
       DIMENSION SICP(NMAX,MAX_OCC)
       DIMENSION FMAT(MAX_OCC,MAX_OCC,4,2),RHI(4)

      REAL*8,allocatable :: PSIG(:,:,:), PTS(:,:), PSID(:,:,:) 
     &    ,GRAD(:,:,:,:,:),GRAD1(:,:,:,:),RVECA(:,:),GSUB(:,:) 
      LOGICAL ,allocatable :: ICOUNT(:,:) 
C
       DIMENSION RR(3),JBEG(3),JEND(3)
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
       DATA FIRST/.TRUE./
       DATA JBEG,JEND/1,2,5,1,4,10/

       INTEGER MXXD,KXXS,KXXO


       SAVE  !YY to save or not save?

!YY allocate tmp2 arrays
       ALLOCATE(PSIG(NMAX,10,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING PSIG'
       ALLOCATE(PSID(3,MAX_OCC,NMAX),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING PSID'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING PTS'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING GRAD'
       ALLOCATE(GRAD1(3,MAX_CON,10,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING GRAD1'
       ALLOCATE(GSUB(3,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING GSUB'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING RVECA'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR ALLOCATING ICOUNT'

       LSPN = NSPN
!      ISPFAC = 2/NSPN
!      spfac = 2.0d0/float(nspn)
!      spfac = 1.0d0
c      if(lfrm(2).ne.0) then
c        print *, 'lfrm(2) ne 0', lfrm
c        stop
c      end if
       NGRAD=4
       FORALL (IWF=1:MAX_OCC,IPT=1:NMAX,IGR=1:10)
         PSIG(IPT,IGR,IWF) = 0.0D0
       END FORALL
c
c      Get correct index for local orbital wave function
c
       IID=0
       DO 86 IFNCT=1,NFNCT
         LMAX1=LSYMMAX(IFNCT)+1
C
C FOR ALL POSITIONS OF THIS CENTER
C
         DO 84 I_POS=1,N_POS(IFNCT)
          IID=IID+1
C
C GET SYMMETRY INFO
c
C
          CALL OBINFO(1,RIDT(1,IID),RVECA,M_NUC,ISHDUM)
          IF(NWF.GT.MAX_OCC)THEN
           PRINT *,'APTSLV: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C UNRAVEL2 gives coefficients of local orbitals in fixed basis
C
          MXXD=1        !Get coefs for occupied orbitals only -- both spins
          KXXS=LSPN     !
          KXXO=0        !This should not be needed
          DO 82 J_POS=1,M_NUC
           CALL UNRAVEL2(MXXD,KXXS,KXXO,IFNCT,IID,J_POS,
     &                   RIDT(1,IID),RVECA,L_NUC,1)
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
C
C GET ORBITS AND DERIVATIVES for local orbital lfn
C
            NDERV=1
            CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            IF (IUPDAT) THEN
             IPT=JPTS-1
             ILOC=0
             DO 78 LI=1,LMAX1
              DO MU=1,ISIZE(LI)
               DO ICON=1,N_CON(LI,IFNCT)
                ILOC=ILOC+1
                IF (ICOUNT(ICON,LI)) THEN
                 DO LFM=1,LFRM(1)+LFRM(2) !LFRM(LSPN) !now get all orbitals of the same spin at once
                  LFN = LFM !+ (LSPN-1)*LFRM(1)
                  FACTOR=PSI(ILOC,LFN,1)  !note LFN here, but LFM below in PSIG
                  if(abs(FACTOR) .GT. 1.0d-10) then
                  DO IGR=1,NGRAD
                   DO LPV=1,NPV
                    PSIG(IPT+LPV,IGR,LFM)=PSIG(IPT+LPV,IGR,LFM)
     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
                   END DO    !LPV:  fill in NSPEED points at a time
                  END DO     !IGR
                  end if
                 END DO      !IWF
                END IF
               END DO        !ICON
              END DO         !MU
   78        CONTINUE        !LI
            END IF
   80      CONTINUE          !JPTS   Loop over all NMAX points
   82     CONTINUE           !JPOS
   84    CONTINUE            !IPOS
   86   CONTINUE             !IFNCT
c
C  The following section computes dphi_i/dR = Sum_s {c^i_s * df_s/dR}
c      for each atomic coordinate
c  As each is built, it is used to compute matrix elements
c        <dphi_i/dR|VSIC_i|phi_i> 
C
c  psid stores dphi_i/dR
c  psig stores phi_i
c
C  IID is the atom index of the derivative;psid has an implicit IID index
c
      IID=0
      DO 485 IFNCT=1,NFNCT
       DO 483 I_POS=1,N_POS(IFNCT)
        IID=IID+1
        CALL OBINFO(1,RIDT(1,IID),RVECA,M_NUC,ISHDUM)
        PSID = 0.0D0
        J_POS=1
c
c  unravel2 returns coefficients for local orbitals
c
         MDX=1             !get coefs for all occupied orbitals
         MDUM=1            !dummy index when mdx = 1
         CALL UNRAVEL2(MDX,LSPN,MDUM,IFNCT,IID,J_POS,RIDT(1,IID),
     & RVECA,L_NUC,1)
        IF(L_NUC.NE.M_NUC)THEN
         write(6,*)'FRCSLV: PROBLEM IN UNRAVEL'
         CALL STOPIT
        END IF
        DO 495 IPT=1,MPTS
         DO 493 IGRP=1,NGRP
C
C ROTATE POINT WHEN SYMMETRY is used
c  ***NOTE*** point group symmetry is not implemented yet 1-28-2018 KAJ
c 
C
          DO I=1,3
           RR(I)=0.0D0
           DO J=1,3
            RR(I)=RR(I)+RMAT(I,J,IGRP)*RPTS(J,IPT)
           END DO
          END DO
          ILOC=0
C
C  INITIALIZE SITE DERIVS
C
          DO IWF=1,NWF
           DO IX=1,3
            GSUB(IX,IWF)=0.0D0
           END DO
          END DO
          X=RR(1)-RVECA(1,J_POS)
          Y=RR(2)-RVECA(2,J_POS)
          Z=RR(3)-RVECA(3,J_POS)
          CALL GTGRAD(X,Y,Z,IFNCT,GRAD1)
          DO 480 LI=1,3
           DO 478 MU=JBEG(LI),JEND(LI)
            DO 475 ICON=1,N_CON(LI,IFNCT)
             ILOC=ILOC+1
             DO 465 IWF=1,lfrm(1)+lfrm(2)   !NWF
              GSUB(1,IWF)=GSUB(1,IWF)
     &                   +PSI(ILOC,IWF,1)*GRAD1(1,ICON,MU,LI)
              GSUB(2,IWF)=GSUB(2,IWF)
     &                   +PSI(ILOC,IWF,1)*GRAD1(2,ICON,MU,LI)
              GSUB(3,IWF)=GSUB(3,IWF)
     &                   +PSI(ILOC,IWF,1)*GRAD1(3,ICON,MU,LI)
  465        CONTINUE
  475       CONTINUE
  478      CONTINUE
  480     CONTINUE
C
C UPDATE PSID:  implicit iid index;  lspn in call statement
C
c
c  The following lines are needed if there is point group symmetry
c  SKIPPING OVER THEM UNTIL THE CODE IS SYMMETRIZED
c
          if(n_rep.ne.1) then
          IWF=0
           DO K_REP=1,N_REP
            DO IOCC=1,N_OCC(K_REP,LSPN)  !n_occ = # of basis functions?
             JWF_BEG=IWF
             DO IROW=1,NDMREP(K_REP)
              IWF=IWF+1
              JWF=JWF_BEG
              DO JROW=1,NDMREP(K_REP)
               JWF=JWF+1
               PSID(1,IWF,IPT)=PSID(1,IWF,IPT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GSUB(1,JWF)
               PSID(2,IWF,IPT)=PSID(2,IWF,IPT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GSUB(2,JWF)
               PSID(3,IWF,IPT)=PSID(3,IWF,IPT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GSUB(3,JWF)
              END DO     !JROW
             END DO      !IROW
            END DO       !IOCC
           END DO        !K_REP
         else 
          do iwf = 1, lfrm(1) + lfrm(2)
           do ix = 1,3
             psid(ix,iwf,ipt) = gsub(ix,iwf)
           end do
          end do
         end if
C
C END OF UPDATE:
C
  493    CONTINUE
  495   CONTINUE
C
C  UPDATE MATRIX ELEMENTS
C
!
! Note PSIG index should be LFM; PSID index should be LFN
!  PSID is currently being computed for orbitals of both spins KAJ 2-14-23
!  
        DO LFM=1,LFRM(1)+LFRM(2) !LFRM(LSPN)
         LFN = LFM !+ (LSPN-1)*LFRM(1)
         DO 3100 IPT=1,MPTS
           DO 391 IX=1,3
            DERSIC(IX,LFN,IID)=DERSIC(IX,LFN,IID) + 
     &       PSID(IX,LFN,IPT)*SICP(IPT,LFM)*PSIG(IPT,1,LFM)*WMSA(IPT)
  391      CONTINUE
 3100    CONTINUE
        END DO  !LFM
  483  CONTINUE  !I_POS
  485 CONTINUE   !IFNCT
c
 1000  FORMAT(' ',i5,5F12.5)
c
!YY Deallocate tmp2
c
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING PSIG'
       DEALLOCATE(PSID,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING PSID'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING NSPEED'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING GRAD'
       DEALLOCATE(GRAD1,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING GRAD1'
       DEALLOCATE(GSUB,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING GSUB'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING RVECA'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'SICDERSLV:ERROR DEALLOCATING ICOUNT'
       RETURN
       END
