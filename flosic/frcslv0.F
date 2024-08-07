C UTEP Electronic Structure Lab (2020)
       SUBROUTINE FRCSLV0(MPTS,LPTS)
c
c     The following is how I call this subroutine in the current serial version 
c     SUBROUTINE FRCSLV(MPTS,LPTS,over1,ek,dftV,allH)  !serial version 
c     I want to change this so that we call it using the MPI approach.
C
C
C 02/19/97 David Clay Patton
C
      use global_inputs,only : LIBXC1
      use mixpot1,only : VXC=>POTIN,POT=>POTOUT
      use mesh1,only : wmsh,rmsh,nmsh
      use common2,only : RIDT, N_CON, LSYMMAX, N_POS,
     &  NFNCT, ISPN, NSPN, FRC1, FRC2, over1, ek, dftV,allH
      use common3,only : RMAT, NGRP
      use common5,only : N_OCC, PSI, NWF, NWFS, EVLOCC
      use common8,only : REP, N_REP, NDMREP
      use xtmp2a, only : MIXINS,ISMGGA
      use SICMAT,only  : SIC
      use FRM,only     : LFRM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MPTS, LPTS, I, I_CON, I_LOC, I_POS, I_PT, ICON, IERR,
     & IFNCT, IGRP, IID, ILOC, IOCC, IPT, IPTS, IROW, ISHDUM, ISIZE,
     & ISPFAC, IWF, IX, J, J_POS, JBEG, JEND, JPTS, JROW, JWF, JWF_BEG,
     & K_REP, L_NUC, LI, LMAX1, LPTS_BEG, LPV, M_NUC, MU, NMAX, NPV,
     & kspn, kwf, mdum, mdx, mwf,nwftot
       REAL*8 :: SYMBOL, FACTOR, GGRAD, GGSUB, GNAD, GRAD, GSUB, PTS, 
     & RR, RVECA, X, Y, Z, term_ibp
      SAVE
      PARAMETER (NMAX=MPBLOCK)
      INTEGER :: IOFS, IOFS3
      LOGICAL ICOUNT,IUPDAT
C     COMMON/TMP2/PSIN(MAX_OCC,NMAX),PSIG(3,MAX_OCC,NMAX)
C     & ,PSIL(MAX_OCC,NMAX)
      REAL*8,ALLOCATABLE :: PSIN(:,:),PSIG(:,:,:),PSIL(:,:)
      real*8,allocatable :: PPSIG(:,:,:),PSIR(:,:,:)
      real*8,allocatable :: PPSIGX(:,:),PPSIGY(:,:),PPSIGZ(:,:)
      real*8,allocatable :: FRC3(:,:),FRC4(:,:)
      REAL*8             :: vxcmrot(6)   !vxcmrot(1:3) [vxcmrot(4:6)] is rotated vxcm for spin up [down]
C     COMMON/MIXPOT/VXC(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
      DIMENSION GSUB(3,MAX_OCC)
      DIMENSION GRAD(3,MAX_CON,10,3)
      DIMENSION GGSUB(6,MAX_OCC)
! YY. GGRAD are derivative terms: d/dRv d/dx_i phi's
      DIMENSION GGRAD(6,MAX_CON,10,3)  
      DIMENSION RVECA(3,MX_GRP)
     & ,ISIZE(3),RR(3),JBEG(3),JEND(3)
      DIMENSION PTS(NSPEED,3),GNAD(NSPEED,10,6,MAX_CON,3) 
      DIMENSION ICOUNT(MAX_CON,3) 
      DATA JBEG,JEND/1,2,5,1,4,10/
      DATA ISIZE/1,3,6/
C
      IF (MPTS .GT. NMAX) THEN
       write(6,*)'FRCSLV: NMAX MUST BE AT LEAST: ',MPTS
       CALL STOPIT
      END IF
C
C ALLOCATE LOCAL ARRAYS
C
      ALLOCATE(PSIN(MAX_OCC,NMAX),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR ALLOCATING PSIN'
      ALLOCATE(PSIG(3,MAX_OCC,NMAX),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR ALLOCATING PSIG'
      ALLOCATE(PSIL(MAX_OCC,NMAX),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR ALLOCATING PSIL'
      IF(LIBXC1.OR.ISMGGA)THEN
       allocate(PPSIG(6,MAX_OCC,NMAX))
       allocate(PPSIGX(MAX_OCC,NMAX))
       allocate(PPSIGY(MAX_OCC,NMAX))
       allocate(PPSIGZ(MAX_OCC,NMAX))
       allocate(PSIR(MAX_OCC,NMAX,3))
      ENDIF 

      ISPFAC=2/NSPN
      LPTS_BEG=LPTS
C
c  nwftot = total number of occupied orbitals
c
      nwftot = lfrm(1) + lfrm(2)
c
c  Print lagrange multiplier matrix for test
c
c     if(lpts.eq.0) then
c      print *, 'frcslv'
c      print *, 'LM Matrix'
c      do ispn=1,nspn
c       do i =1,lfrm(ispn)
c        print 27, (sic(j,i,ispn),j=1,lfrm(ispn))
c       end do
c      end do
cc     print *, 'frcslv nwf nwfs', nwf, nwfs(1),nwfs(2)
c      print *, 'frcslv n_occ', (n_occ(k_rep,1),k_rep=1,5) 
c      print *, 'frcslv n_occ', (n_occ(k_rep,2),k_rep=1,5) 
c     end if
c27   format(1x, 6(f10.4,1x))

! Zero out FRC3 and FRC4 here
      IID=0
      !MAX_IDENT comes from PARAMS
      IF(LIBXC1.OR.ISMGGA)THEN
       allocate(FRC3(3,MAX_IDENT),FRC4(3,MAX_IDENT))
       FRC3(:,:)=0.0d0
       FRC4(:,:)=0.0d0
      END IF
C
C  CALCULATE WAVEFUNCTIONS, LAPLACIANS on current block of points
C  STORE WAVEFUNCTIONS IN PSIN AND NEGATIVE LAPLACIAN IN PSIL
C
c  PSIN = phi_i;  PSIL = -del^2/2 (phi_i)
c
      DO IPT=1,MPTS
       DO IWF=1,NWFTOT    
        PSIN(IWF,IPT)=0.0D0
        PSIL(IWF,IPT)=0.0D0
       END DO
      END DO
      IF(LIBXC1.OR.ISMGGA) PSIR(:,:,:)=0.0D0
!      PSIR(:,:,:)=0.0D0
      IID=0
      DO 386 IFNCT=1,NFNCT
       LMAX1=LSYMMAX(IFNCT)+1
       DO 384 I_POS=1,N_POS(IFNCT)
        IID=IID+1
        CALL OBINFO(1,RIDT(1,IID),RVECA,M_NUC,ISHDUM)
        DO 382 J_POS=1,M_NUC
c
c UNRAVEL2 returns LOCAL ORBITAL coefficients
c
         MDX=1             !=1 gets evecs only for occ orbitals
         KSPN=NSPN
         MDUM=1
         CALL UNRAVEL2(MDX,KSPN,MDUM,IFNCT,IID,J_POS,RIDT(1,IID),
     &                 RVECA,L_NUC,1)
         IF(NWFTOT.GT.MAX_OCC)THEN
          write(6,*)'FRCSLV: MAX_OCC MUST BE AT LEAST:',NWFTOT
          CALL STOPIT
         END IF
         IF(L_NUC.NE.M_NUC)THEN
          write(6,*)'FRCSLV: PROBLEM IN UNRAVEL'
          CALL STOPIT
         END IF
         LPTS=LPTS_BEG
         DO 380 JPTS=1,MPTS,NSPEED
          NPV=MIN(NSPEED,MPTS-JPTS+1)
          DO LPV=1,NPV  
           PTS(LPV,1)=RMSH(1,LPTS+LPV)-RVECA(1,J_POS)
           PTS(LPV,2)=RMSH(2,LPTS+LPV)-RVECA(2,J_POS)
           PTS(LPV,3)=RMSH(3,LPTS+LPV)-RVECA(3,J_POS)
          END DO
          CALL GORBDRV(2,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GNAD)
          IF(IUPDAT)THEN
           IPTS=JPTS-1 
           ILOC=0
           DO 378 LI=1,LMAX1
            DO MU=1,ISIZE(LI)
             DO ICON=1,N_CON(LI,IFNCT)
              ILOC=ILOC+1
              IF(ICOUNT(ICON,LI))THEN
               DO  IWF=1,NWFTOT         !sum over all occ wfs
                FACTOR=PSI(ILOC,IWF,1)
                DO LPV=1,NPV
                 PSIN(IWF,IPTS+LPV)=PSIN(IWF,IPTS+LPV) 
     &           +      FACTOR* GNAD(LPV,1,MU,ICON,LI)   
                 PSIL(IWF,IPTS+LPV)=PSIL(IWF,IPTS+LPV)
     &           -0.5D0*FACTOR*(GNAD(LPV,5,MU,ICON,LI)+  
     &                          GNAD(LPV,6,MU,ICON,LI)+  
     &                          GNAD(LPV,7,MU,ICON,LI))  
! YY. I need first derivatives: d/dx, d/dy, and d/dz
                 IF(LIBXC1.OR.ISMGGA)THEN
                  PSIR(IWF,IPTS+LPV,1)=PSIR(IWF,IPTS+LPV,1)
     &            +      FACTOR* GNAD(LPV,2,MU,ICON,LI)
                  PSIR(IWF,IPTS+LPV,2)=PSIR(IWF,IPTS+LPV,2)
     &            +      FACTOR* GNAD(LPV,3,MU,ICON,LI)
                  PSIR(IWF,IPTS+LPV,3)=PSIR(IWF,IPTS+LPV,3)
     &            +      FACTOR* GNAD(LPV,4,MU,ICON,LI)
                 END IF
                END DO
               END DO
              END IF
             END DO
            END DO
  378      CONTINUE  
          END IF 
          LPTS=LPTS+NPV
  380    CONTINUE
  382   CONTINUE
  384  CONTINUE
  386 CONTINUE
C
C  LOOP OVER SHELLS TO GET DERIVATIVES:  ANOTHER STAB AT GETTING 
C  THE SYM RIGHT
C
C  This section computes dphi_i/dR = Sum_s {c^i_s * df_s/dR}
c      for each atomic coordinate
c  As each is built, it is used to compute matrix elements
c        <dphi_i/dR|H_i|phi_i> and <dphi_i/dR|phi_j>
C
C  IID is the atom index of the derivative;
c
      IID=0
      DO 485 IFNCT=1,NFNCT
       DO 483 I_POS=1,N_POS(IFNCT)
        IID=IID+1
        CALL OBINFO(1,RIDT(1,IID),RVECA,M_NUC,ISHDUM)
C
C  IDENTIFY ATOM SHELL
C  STORE SYMMETRIZED GRADIENTS IN PSIG
C
        DO I_PT=1,MPTS
         DO IWF=1,NWFTOT    
          DO IX=1,3
           PSIG(IX,IWF,I_PT)=0.0D0
          END DO
         END DO
        END DO
        IF(LIBXC1.OR.ISMGGA)THEN
         PPSIG(:,:,:)=0.0d0
         PPSIGX(:,:)=0.0d0
         PPSIGY(:,:)=0.0d0
         PPSIGZ(:,:)=0.0d0
        ENDIF
        J_POS=1
c
c  unravel2 returns coefficients for local orbitals
c
        MDX=1
        KSPN=NSPN
        MDUM=1
        CALL UNRAVEL2(MDX,KSPN,MDUM,IFNCT,IID,J_POS,RIDT(1,IID),
     &                RVECA,L_NUC,1)
        IF(L_NUC.NE.M_NUC)THEN
         write(6,*)'FRCSLV: PROBLEM IN UNRAVEL'
         CALL STOPIT
        END IF
        LPTS=LPTS_BEG
        DO 495 I_PT=1,MPTS
         LPTS=LPTS+1
         DO 493 IGRP=1,NGRP
C
C ROTATE POINT
C MUST BE RMAT(I,J,IGRP) DUE TO SYMMETRIZING
C
          DO I=1,3
           RR(I)=0.0D0
           vxcmrot(I)=0.0d0
           if(NSPN == 2 ) vxcmrot(I+3)=0.0d0
           DO J=1,3
            RR(I)=RR(I)+RMAT(I,J,IGRP)*RMSH(J,LPTS)
            IF(LIBXC1.OR.ISMGGA)THEN
             !Rotate vxcm as well
             vxcmrot(I)=vxcmrot(I)+RMAT(I,J,IGRP)*MIXINS(J,I_PT) !Check this indeces
             if(NSPN == 2) vxcmrot(I+3)=vxcmrot(I+3)
     &                      +RMAT(I,J,IGRP)*MIXINS(J,I_PT+MPTS)
            ENDIF
           END DO
          END DO
          I_LOC=0
C
C  INITIALIZE SITE DERIVS
C
          DO IWF=1,NWFTOT    
           DO IX=1,3
            GSUB(IX,IWF)=0.0D0
           END DO
          END DO
          IF(LIBXC1.OR.ISMGGA) GGSUB(:,:)=0.0d0
          X=RR(1)-RVECA(1,J_POS)
          Y=RR(2)-RVECA(2,J_POS)
          Z=RR(3)-RVECA(3,J_POS)
          IF(LIBXC1.OR.ISMGGA) THEN
           CALL GTGRADDRV(X,Y,Z,IFNCT,GRAD,GGRAD)
          ELSE
           CALL GTGRAD(X,Y,Z,IFNCT,GRAD) 
          END IF
          DO 480 LI=1,3
           DO 478 MU=JBEG(LI),JEND(LI)
            DO 475 I_CON=1,N_CON(LI,IFNCT)
             I_LOC=I_LOC+1
             DO 465 IWF=1,NWFTOT
              GSUB(1,IWF)=GSUB(1,IWF)
     &                   +PSI(I_LOC,IWF,1)*GRAD(1,I_CON,MU,LI)
              GSUB(2,IWF)=GSUB(2,IWF)
     &                   +PSI(I_LOC,IWF,1)*GRAD(2,I_CON,MU,LI)
              GSUB(3,IWF)=GSUB(3,IWF)
     &                   +PSI(I_LOC,IWF,1)*GRAD(3,I_CON,MU,LI)
! sum of basis set functions
              IF(LIBXC1.OR.ISMGGA) THEN
               GGSUB(1,IWF)=GGSUB(1,IWF)+PSI(I_LOC,IWF,1)
     &                     *GGRAD(1,I_CON,MU,LI)
               GGSUB(2,IWF)=GGSUB(2,IWF)+PSI(I_LOC,IWF,1)
     &                     *GGRAD(2,I_CON,MU,LI)
               GGSUB(3,IWF)=GGSUB(3,IWF)+PSI(I_LOC,IWF,1)
     &                     *GGRAD(3,I_CON,MU,LI)
               GGSUB(4,IWF)=GGSUB(4,IWF)+PSI(I_LOC,IWF,1)
     &                     *GGRAD(4,I_CON,MU,LI)
               GGSUB(5,IWF)=GGSUB(5,IWF)+PSI(I_LOC,IWF,1)
     &                     *GGRAD(5,I_CON,MU,LI)
               GGSUB(6,IWF)=GGSUB(6,IWF)+PSI(I_LOC,IWF,1)
     &                     *GGRAD(6,I_CON,MU,LI)
              END IF
  465        CONTINUE
  475       CONTINUE
  478      CONTINUE
  480     CONTINUE
C
C UPDATE PSIG: note that psig has an implicit iid index
C
          IWF=0
          DO ISPN=1,NSPN
           IOFS3=3*(ISPN-1) ! Offset for rotated vxcm array
           DO K_REP=1,N_REP
            DO IOCC=1,N_OCC(K_REP,ISPN)
             JWF_BEG=IWF
             DO IROW=1,NDMREP(K_REP)
              IWF=IWF+1
              JWF=JWF_BEG
              DO JROW=1,NDMREP(K_REP)
               JWF=JWF+1
C
C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Look Here...
C The 3 Occurences of REP(IROW,JROW...) might be REP(JROW,IROW....)
C
               PSIG(1,IWF,I_PT)=PSIG(1,IWF,I_PT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GSUB(1,JWF)
               PSIG(2,IWF,I_PT)=PSIG(2,IWF,I_PT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GSUB(2,JWF)
               PSIG(3,IWF,I_PT)=PSIG(3,IWF,I_PT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GSUB(3,JWF)

               IF(LIBXC1.OR.ISMGGA)THEN
                PPSIG(1,IWF,I_PT)=PPSIG(1,IWF,I_PT)
     &           +REP(JROW,IROW,IGRP,K_REP)*GGSUB(1,JWF)
                PPSIG(2,IWF,I_PT)=PPSIG(2,IWF,I_PT)
     &           +REP(JROW,IROW,IGRP,K_REP)*GGSUB(2,JWF)
                PPSIG(3,IWF,I_PT)=PPSIG(3,IWF,I_PT)
     &           +REP(JROW,IROW,IGRP,K_REP)*GGSUB(3,JWF)
                PPSIG(4,IWF,I_PT)=PPSIG(4,IWF,I_PT)
     &           +REP(JROW,IROW,IGRP,K_REP)*GGSUB(4,JWF)
                PPSIG(5,IWF,I_PT)=PPSIG(5,IWF,I_PT)
     &           +REP(JROW,IROW,IGRP,K_REP)*GGSUB(5,JWF)
                PPSIG(6,IWF,I_PT)=PPSIG(6,IWF,I_PT)
     &           +REP(JROW,IROW,IGRP,K_REP)*GGSUB(6,JWF)

                PPSIGX(IWF,I_PT)=PPSIGX(IWF,I_PT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(1,JWF)*vxcmrot(1+IOFS3)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(4,JWF)*vxcmrot(2+IOFS3)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(5,JWF)*vxcmrot(3+IOFS3)

                PPSIGY(IWF,I_PT)=PPSIGY(IWF,I_PT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(4,JWF)*vxcmrot(1+IOFS3)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(2,JWF)*vxcmrot(2+IOFS3)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(6,JWF)*vxcmrot(3+IOFS3)

                PPSIGZ(IWF,I_PT)=PPSIGZ(IWF,I_PT)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(5,JWF)*vxcmrot(1+IOFS3)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(6,JWF)*vxcmrot(2+IOFS3)
     &          +REP(JROW,IROW,IGRP,K_REP)*GGSUB(3,JWF)*vxcmrot(3+IOFS3)
               END IF
              END DO
             END DO
            END DO
           END DO
          END DO
          IF ((IWF.NE.NWF).OR.(JWF.NE.NWF)) THEN
           write(6,*)'FRCSLV: OOPS: IWF, JWF AND NWF NOT EQUAL'
           CALL STOPIT
          END IF
C
C END OF UPDATE:
C
  493    CONTINUE
  495   CONTINUE
C
C  UPDATE MATRIX ELEMENTS
C  JPTS = 1, 2*nmsh -- needed for spin up and spin down potential 
c  LPTS = 1, nmsh -- used for wave functions
c
c  Here FRC1 contains sum_i<del(phi_i)/del(R_nu)|HDFT|phi_i>.  Note that
c   it is missing the SIC part that is computed in siclag_der.
c
c  Here FRC2 contains sum_ij Lambda_ij <del(phi_i)/del(R_nu)|phi_j>
C
        DO 30 ISPN=1,NSPN
         JPTS=LPTS_BEG + NMSH*(ISPN-1)
         LPTS=LPTS_BEG
         IOFS=(ISPN-1)*MPTS !Needed for MIXINS array spin offset
         DO 3100 I_PT=1,MPTS
          JPTS=JPTS+1
          LPTS=LPTS+1
          IWF=LFRM(1)*(ISPN-1)
          DO 390 JWF=1,LFRM(ISPN)
           IWF=IWF+1
c          write(6,*) 'E-VALUES',IWF, EVLOCC(IWF)
           DO 391 IX=1,3
            FRC1(IX,IID)=FRC1(IX,IID) + ISPFAC*PSIG(IX,IWF,I_PT)
     &       *(PSIL(IWF,I_PT) + POT(JPTS)*PSIN(IWF,I_PT))*WMSH(LPTS)
!           FRC2(IX,IID)=FRC2(IX,IID) + ISPFAC*PSIG(IX,IWF,I_PT)
!    &       *EVLOCC(IWF)*PSIN(IWF,I_PT)*WMSH(LPTS)
            mwf = LFRM(1)*(ISPN -1)
            do kwf = 1,LFRM(ISPN)
             mwf = mwf + 1
c
c  spin indexing of SIC arrary requires using
c  jwf and kwf for spin down
c
             FRC2(IX,IID)=FRC2(IX,IID) + ISPFAC*PSIG(IX,IWF,I_PT)
     &        *SIC(jwf,kwf,ispn)*PSIN(MWF,I_PT)*WMSH(LPTS)
            end do
  391      CONTINUE
           IF(LIBXC1.OR.ISMGGA)THEN

! GGA libxc/Ham. mixing first and second term combined
            FRC3(1,IID)=FRC3(1,IID)
     &       + ISPFAC* PPSIGX(IWF,I_PT)*PSIN(IWF,I_PT)*WMSH(LPTS)
     &       + ISPFAC*(PSIR(IWF,I_PT,1)*MIXINS(1,I_PT+IOFS) +
     &                 PSIR(IWF,I_PT,2)*MIXINS(2,I_PT+IOFS) +
     &                 PSIR(IWF,I_PT,3)*MIXINS(3,I_PT+IOFS))
     &                *PSIG(1,IWF,I_PT)*WMSH(LPTS)

            FRC3(2,IID)=FRC3(2,IID)
     &       + ISPFAC* PPSIGY(IWF,I_PT)*PSIN(IWF,I_PT)*WMSH(LPTS)
     &       + ISPFAC*(PSIR(IWF,I_PT,1)*MIXINS(1,I_PT+IOFS) +
     &                 PSIR(IWF,I_PT,2)*MIXINS(2,I_PT+IOFS) +
     &                 PSIR(IWF,I_PT,3)*MIXINS(3,I_PT+IOFS))
     &                *PSIG(2,IWF,I_PT)*WMSH(LPTS)

            FRC3(3,IID)=FRC3(3,IID)
     &        + ISPFAC* PPSIGZ(IWF,I_PT)*PSIN(IWF,I_PT)*WMSH(LPTS)
     &        + ISPFAC*(PSIR(IWF,I_PT,1)*MIXINS(1,I_PT+IOFS) +
     &                  PSIR(IWF,I_PT,2)*MIXINS(2,I_PT+IOFS) +
     &                  PSIR(IWF,I_PT,3)*MIXINS(3,I_PT+IOFS))
     &                 *PSIG(3,IWF,I_PT)*WMSH(LPTS)

! MGGA
            FRC4(1,IID)=FRC4(1,IID) + ISPFAC*
     &                (PPSIG(1,IWF,I_PT)*PSIR(IWF,I_PT,1)+ 
     &                 PPSIG(4,IWF,I_PT)*PSIR(IWF,I_PT,2)+ 
     &                 PPSIG(5,IWF,I_PT)*PSIR(IWF,I_PT,3))*
     &                 MIXINS(4,I_PT+IOFS)*WMSH(LPTS)

            FRC4(2,IID)=FRC4(2,IID) + ISPFAC*
     &                (PPSIG(4,IWF,I_PT)*PSIR(IWF,I_PT,1)+ 
     &                 PPSIG(2,IWF,I_PT)*PSIR(IWF,I_PT,2)+ 
     &                 PPSIG(6,IWF,I_PT)*PSIR(IWF,I_PT,3))*
     &                 MIXINS(4,I_PT+IOFS)*WMSH(LPTS)

            FRC4(3,IID)=FRC4(3,IID) + ISPFAC*
     &                (PPSIG(5,IWF,I_PT)*PSIR(IWF,I_PT,1)+ 
     &                 PPSIG(6,IWF,I_PT)*PSIR(IWF,I_PT,2)+ 
     &                 PPSIG(3,IWF,I_PT)*PSIR(IWF,I_PT,3))*
     &                 MIXINS(4,I_PT+IOFS)*WMSH(LPTS)
           ENDIF

          ! write(*,*) "FRC",IID,FRC3(1,IID)
  390     CONTINUE
 3100    CONTINUE
   30   CONTINUE
  483  CONTINUE
  485 CONTINUE

!YY. Adding FRC3 and FRC4 back
      IF(LIBXC1.OR.ISMGGA)THEN
       IID=0
       DO IFNCT=1,NFNCT
        DO I_POS=1,N_POS(IFNCT)
         IID=IID+1
         FRC1(1,IID)=FRC1(1,IID)+FRC3(1,IID)+FRC4(1,IID)
         FRC1(2,IID)=FRC1(2,IID)+FRC3(2,IID)+FRC4(2,IID)
         FRC1(3,IID)=FRC1(3,IID)+FRC3(3,IID)+FRC4(3,IID)
        END DO
       END DO
      END IF
C
C  TEST calculations: to test that orbitals are created correctly
c  Over1 = overlap matrix in local orbital basis:  nocc x nocc
c  ek = kinetic energy matrix in LO basis
c  dftV = DFT potential matrix in LO basis 
c  allH = DFT hamiltonian matrix in LO basis
c
C
      DO  ISPN=1,NSPN
       JPTS=LPTS_BEG + NMSH*(ISPN-1)
       LPTS=LPTS_BEG
       IOFS=(ISPN-1)*MPTS
       DO I_PT=1,MPTS
        JPTS=JPTS+1
        LPTS=LPTS+1
        IWF=LFRM(1)*(ISPN-1)
        DO  JWF=1,LFRM(ISPN)
         IWF=IWF+1
         do mwf = LFRM(1)*(ispn-1)+1,LFRM(ispn)
!YY. Changed the loop variable by +1. Or segmentation fault.
!<< do mwf = LFRM(1)*(nspn-1),LFRM(ispn)
!>> do mwf = LFRM(1)*(ispn-1)+1,LFRM(ispn)
          over1(iwf,mwf) = over1(iwf,mwf)+PSIN(IWF,I_PT)*
     &                     PSIN(MWF,I_PT)*wmsh(lpts)
          ek(iwf,mwf) = ek(iwf,mwf)+PSIL(IWF,I_PT)*
     &                  PSIN(MWF,I_PT)*wmsh(lpts)
          dftV(mwf,iwf) = dftV(mwf,iwf) +  (PSIN(mwf,i_pt)
     &               *pot(jpts))*PSIN(IWF,I_PT)*WMSH(LPTS)
          allH(iwf,mwf) = allH(iwf,mwf)+(PSIL(mWF,I_PT)+
     &                    PSIN(mwf,I_PT)*POT(JPTS))
     &                   *PSIN(iWF,I_PT)*wmsh(lpts)
 
          if(LIBXC1.OR.ISMGGA) then
            term_ibp = 
     &                (PSIR(MWF,I_PT,1)*MIXINS(1,I_PT+IOFS) +
     &                 PSIR(MWF,I_PT,2)*MIXINS(2,I_PT+IOFS) +
     &                 PSIR(MWF,I_PT,3)*MIXINS(3,I_PT+IOFS))
     &                *PSIN(IWF,I_PT)*WMSH(LPTS)
     &               +
     &                (PSIR(IWF,I_PT,1)*MIXINS(1,I_PT+IOFS) +
     &                 PSIR(IWF,I_PT,2)*MIXINS(2,I_PT+IOFS) +
     &                 PSIR(IWF,I_PT,3)*MIXINS(3,I_PT+IOFS))
     &                *PSIN(MWF,I_PT)*WMSH(LPTS)        
     &               +
     &                (PSIR(MWF,I_PT,1)*PSIR(IWF,I_PT,1)+
     &                 PSIR(MWF,I_PT,2)*PSIR(IWF,I_PT,2)+
     &                 PSIR(MWF,I_PT,3)*PSIR(IWF,I_PT,3))*
     &                 MIXINS(4,I_PT+IOFS)*WMSH(LPTS)

             dftV(mwf,iwf) = dftV(mwf,iwf) +  term_ibp
             allH(iwf,mwf) = allH(iwf,mwf) +  term_ibp
           end if
         end do
        end do
       end do
      end do
C
C DEALLOCATE LOCAL ARRAYS
C
      IF(LIBXC1.OR.ISMGGA)THEN
       deallocate(FRC4)
       deallocate(FRC3)
       deallocate(PSIR)
       deallocate(PPSIGZ)
       deallocate(PPSIGY)
       deallocate(PPSIGX)
       deallocate(PPSIG)
      ENDIF
      DEALLOCATE(PSIN,STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR DEALLOCATING PSIN'
      DEALLOCATE(PSIG,STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR DEALLOCATING PSIG'
      DEALLOCATE(PSIL,STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR DEALLOCATING PSIL'
      RETURN
      END
