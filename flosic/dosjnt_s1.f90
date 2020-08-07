! UTEP Electronic Structure Lab (2020)
SUBROUTINE DOSJNT_S1(IPILE)
  use dosjnt_mod,only : H, PSIG, RVECA, PTS, GRAD, ICOUNT, &
        ISIZE, P, Q, V, FCGRP, CHARGE
  use mpidat1,only : NPROC,NCALLED
  use zero1
  use debug1
  use mesh1,only : WMSH,RMSH,NMSH
  use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT, E_UP, &
        E_DN, ISPN, NSPN, DIPOLE
  use common3,only : RMAT, NGRP
  use common5,only : PSI, NWF, NWFS, EFERMI, EVLOCC
! Conversion to implicit none.  Raja Zope Sun Oct 23 22:35:29 MDT 2016

!       INCLUDE  'PARAMAS'  
        INCLUDE  'PARAMA2'  
  INTEGER,INTENT(IN) :: IPILE
  LOGICAL            :: IUPDAT
  INTEGER,PARAMETER  :: NMAX=MPBLOCK
       INTEGER :: I_POS, ICON, IFNCT, IGP, ILOC, IPTS, ISHDUMMY, &
     & ISHELLA, IWF, IX, J, J_POS, JPTS, JWF, JWMAX, K, L_NUC, LI, &
     & LMAX1, LPV, M_NUC, MPTS, MU, NOFS,NPV
       REAL*8 :: FACTOR 

        NOFS=IPILE*NMAX
        MPTS=MIN(NMAX,NMSH-NOFS)
        DO IPTS=1,MPTS
         Q(IPTS,1)=RMSH(1,IPTS+NOFS)
         Q(IPTS,2)=RMSH(2,IPTS+NOFS)
         Q(IPTS,3)=RMSH(3,IPTS+NOFS)
         V(  IPTS)=WMSH(IPTS+NOFS)*FCGRP
        END DO
        DO 800 IGP=1,NGRP
         DO J=1,3
          DO IPTS=1,MPTS
           P(IPTS,J)=0.0D0
          END DO
          DO K=1,3
           DO IPTS=1,MPTS
            P(IPTS,J)=P(IPTS,J)+RMAT(K,J,IGP)*Q(IPTS,K)
           END DO
          END DO
         END DO
!
! INITIALIZE PSIG 
!
         CALL ZEROREALMATRIX(PSIG,MPTS,NWF)
         ISHELLA=0
         DO 86 IFNCT=1,NFNCT
          LMAX1=LSYMMAX(IFNCT)+1
          DO 84 I_POS=1,N_POS(IFNCT)
           ISHELLA=ISHELLA+1
           CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
           DO 82 J_POS=1,M_NUC
            CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA), &
                        RVECA,L_NUC,1)
            IF(L_NUC.NE.M_NUC)THEN
             PRINT *,'DOSJNT: PROBLEM IN UNRAVEL'
             CALL STOPIT
            END IF
            DO 80 JPTS=1,MPTS,NSPEED
             NPV=MIN(NSPEED,MPTS-JPTS+1)
             IPTS=JPTS-1
             DO LPV=1,NPV
              PTS(LPV,1)=P(IPTS+LPV,1)-RVECA(1,J_POS)
              PTS(LPV,2)=P(IPTS+LPV,2)-RVECA(2,J_POS)
              PTS(LPV,3)=P(IPTS+LPV,3)-RVECA(3,J_POS)
             END DO
             CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
             IF (IUPDAT) THEN
              ILOC=0
              DO 78 LI=1,LMAX1
               DO MU=1,ISIZE(LI)
                DO ICON=1,N_CON(LI,IFNCT)
                 ILOC=ILOC+1
                 IF (ICOUNT(ICON,LI)) THEN
                  DO IWF=1,NWF
                   FACTOR=PSI(ILOC,IWF,1)
                   DO LPV=1,NPV
                    PSIG(IPTS+LPV,IWF)=PSIG(IPTS+LPV,IWF) &
                   +FACTOR*GRAD(LPV,1,MU,ICON,LI)
                   END DO  
                  END DO  
                 END IF
                END DO  
               END DO  
   78         CONTINUE
             END IF
   80       CONTINUE
   82      CONTINUE
   84     CONTINUE
   86    CONTINUE
!
! UPDATE CHARGE AND DIPOLE MATRICES
!
         DO IWF=1,NWF
          DO IPTS=1,MPTS
           CHARGE=CHARGE+V(IPTS)*PSIG(IPTS,IWF)**2
          END DO
         END DO
         JWMAX=NWF
         DO IX=1,3
          DO IWF=1,NWF
!           JWMAX=NWFS(1)
!           IF (IWF.GT.NWFS(1)) JWMAX=NWF
           DO JWF=IWF,JWMAX
            DO IPTS=1,MPTS
             H(JWF,IWF,IX)=H(JWF,IWF,IX) &
            +PSIG(IPTS,IWF)*PSIG(IPTS,JWF)*V(IPTS)*P(IPTS,IX)
            END DO
           END DO
          END DO
         END DO
  800   CONTINUE

      RETURN
END SUBROUTINE DOSJNT_S1
