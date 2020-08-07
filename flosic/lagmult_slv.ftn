C UTEP Electronic Structure Lab (2020)
       SUBROUTINE LAGMULT_SLV(MPTS,LPTS)
c
c  <phi_i|V_i|phi_j> stored in SIC(i,j); need to update this to full
c   lagrange multiplier matrix (lambda_ij)
c   This routine performs the numerical
c  integrals required to add <phi_i|H^DFT|phi_j> to get the full
c  lambda_ij 
c  This routine assumes the DFT potential is stored in array POT
c    
      use mixpot1,only : VXC=>POTIN,POT=>POTOUT
      use mesh1,only : wmsh,rmsh,nmsh
      use common2,only : RIDT, N_CON, LSYMMAX, N_POS,
     &  NFNCT, ISPN, NSPN, FRC1, FRC2, over1, ek, dftV,allH
       use common3,only : RMAT, NGRP
       use common5,only : N_OCC, PSI, NWF, NWFS, EVLOCC
       use common8,only : REP, N_REP, NDMREP
       use SICMAT,only : SIC
       use FRM,only : LFRM
! Converion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MPTS, LPTS, I, I_CON, I_LOC, I_POS, I_PT, ICON, IERR,
     & IFNCT, IGRP, IID, ILOC, IOCC, IPT, IPTS, IROW, ISHDUM, ISIZE,
     & ISPFAC, IWF, IX, J, J_POS, JBEG, JEND, JPTS, JROW, JWF, JWF_BEG,
     & K_REP, L_NUC, LI, LMAX1, LPTS_BEG, LPV, M_NUC, MU, NMAX,npv,
     & kspn, kwf, mdum, mdx, mwf,nwftot
       REAL*8 :: SYMBOL , FACTOR, GNAD, GRAD, GSUB, PTS, RR, RVECA, X,
     & Y, Z
      SAVE
      PARAMETER (NMAX=MPBLOCK)
      LOGICAL ICOUNT,IUPDAT
      REAL*8,ALLOCATABLE :: PSIN(:,:),PSIG(:,:,:),PSIL(:,:)
c
      DIMENSION GSUB(3,MAX_OCC)
      DIMENSION GRAD(3,MAX_CON,10,3) 
      DIMENSION RVECA(3,MX_GRP)
     & ,ISIZE(3),RR(3),JBEG(3),JEND(3)
      DIMENSION PTS(NSPEED,3),GNAD(NSPEED,10,6,MAX_CON,3) 
      DIMENSION ICOUNT(MAX_CON,3) 
      DATA JBEG,JEND/1,2,5,1,4,10/
      DATA ISIZE/1,3,6/
C
      IF (MPTS .GT. NMAX) THEN
       write(6,*)'LAGMULT_SLV: NMAX MUST BE AT LEAST: ',MPTS
       CALL STOPIT
      END IF
C
C ALLOCATE LOCAL ARRAYS
C
      ALLOCATE(PSIN(MAX_OCC,NMAX),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR ALLOCATING PSIN'
C     ALLOCATE(PSIG(3,MAX_OCC,NMAX),STAT=IERR)
C     IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR ALLOCATING PSIG'
      ALLOCATE(PSIL(MAX_OCC,NMAX),STAT=IERR) 
      IF(IERR/=0)WRITE(6,*)'FRCSLV:ERROR ALLOCATING PSIL'
      ISPFAC=2/NSPN
      LPTS_BEG=LPTS
c
C  SIC PART OF LAGRANGE MULTIPLIER MATRIX stored in SIC
C

C
C  CALCULATE WAVEFUNCTIONS, LAPLACIANS on current block of points
C  STORE WAVEFUNCTIONS IN PSIN AND NEGATIVE LAPLACIAN IN PSIL
C
      nwftot = lfrm(1) + lfrm(2)
      DO IPT=1,MPTS
       DO IWF=1,NWF    
        PSIN(IWF,IPT)=0.0D0
        PSIL(IWF,IPT)=0.0D0
       END DO
      END DO
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
          MDX=1
          KSPN=NSPN
          MDUM=1
          CALL UNRAVEL2(MDX,KSPN,MDUM,IFNCT,IID,J_POS,RIDT(1,IID),
     &      RVECA,L_NUC,1)
         IF(NWF.GT.MAX_OCC)THEN
          write(6,*)'FRCSLV: MAX_OCC MUST BE AT LEAST:',NWF
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
               DO  IWF=1,NWFTOT
                FACTOR=PSI(ILOC,IWF,1)
                DO LPV=1,NPV
                 PSIN(IWF,IPTS+LPV)=PSIN(IWF,IPTS+LPV) 
     &           +      FACTOR* GNAD(LPV,1,MU,ICON,LI)   
                 PSIL(IWF,IPTS+LPV)=PSIL(IWF,IPTS+LPV)
     &           -0.5D0*FACTOR*(GNAD(LPV,5,MU,ICON,LI)+  
     &                       GNAD(LPV,6,MU,ICON,LI)+  
     &                       GNAD(LPV,7,MU,ICON,LI))  
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
c
c  Compute <phi_i|H^DFT|phi_j>
c
       DO ISPN = 1,NSPN
        JPTS = LPTS_BEG + NMSH*(ISPN-1)   !JPTS is pts index for pot
        LPTS = LPTS_BEG 
        DO I_PT = 1,MPTS
          JPTS = JPTS + 1
          LPTS = LPTS + 1
c         IWF=NWFS(1)*(ISPN-1)
          IWF=LFRM(1)*(ISPN-1)          !IWF is local orbital index
c
c  Compute matrix elements; use test integrals as needed...
c
          DO JWF=1,LFRM(ISPN)
            IWF=IWF+1
            mwf = LFRM(1)*(ISPN -1)    !mwf is local orbital index 
            do kwf = 1,LFRM(ISPN)
              mwf = mwf + 1
              sic(jwf,kwf,ispn) = sic(jwf,kwf,ispn)+ 
     &          psin(iwf,i_pt)
     &        *(psil(mwf,i_pt) + pot(jpts)*psin(mwf,i_pt))*wmsh(lpts)
c            over1(iwf,mwf) = over1(iwf,mwf)+PSIN(IWF,I_PT)*
c    &            PSIN(MWF,I_PT)*wmsh(lpts)
c            ek(iwf,mwf) = ek(iwf,mwf)+PSIL(IWF,I_PT)*
c    &             PSIN(MWF,I_PT)*wmsh(lpts)
c            dftV(mwf,iwf) = dftV(mwf,iwf) +  (PSIN(mwf,i_pt)
c    &       *pot(jpts))*PSIN(IWF,I_PT)*WMSH(LPTS)
c            allH(iwf,mwf) = allH(iwf,mwf)+(PSIL(mWF,I_PT)+
c    &       PSIN(mwf,I_PT)*POT(JPTS))
c    &      *PSIN(iWF,I_PT)*wmsh(lpts)
            end do     !kwf
           end do      !jwf
          end do       !i_pt
         end do        !ispn
c
c  end <phi_i|H^DFT|phi_j>
c
 
c  Over1 = overlap matrix in local orbital basis:  nocc x nocc
c  ek = kinetic energy matrix elements
c  dftV = dft potential matrix elements
c  allH = full DFT matrix elements
c
C DEALLOCATE LOCAL ARRAYS
C
      DEALLOCATE(PSIN,STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'LAGMULT_SLV:ERROR DEALLOCATING PSIN'
      DEALLOCATE(PSIL,STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'LAGMULT_SLV:ERROR DEALLOCATING PSIL'
      RETURN
      END
