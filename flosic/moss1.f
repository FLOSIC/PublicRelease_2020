C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE  MOSS1(RHOFERMI,EFGRAD,CENTER,NSPHERE)
C
C COMPUTE EXPECTATION VALUES OF X**2/R**5, X*Y/R**5, ETC. FOR USE
C IN FINDING ELECTRIC FIELD GRADIENTS AT THE NUCLEAR SITES THIS IS USED
C FOR CALCULATING MOSSBAUER PARAMETERS.
C
       use mesh1,only : wmsh,rmsh,nmsh
       use debug1
       use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT, ISPN, NSPN
       use common3,only : RMAT, NGRP
       use common5,only : PSI, NWF, NWFS
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LMXX, NSPHERE, MAXANG, MAXSPH, I, I_POS, ICON, IDT,
     & IFNCT, ILOC, IPTS, IROT, ISHDUMMY, ISHELLA, ISIZE, IWF, J,
     & J_POS, JPTS, JWF, L_NUC, LI, LMX1, LPBEG, LPTS, LPV, LSIZ,
     & M_NUC, MAXRAD, MPTS, MSITES, MU, NMAX, NPV, NXTRA
       REAL*8 :: SYMBOL , RHOFERMI, EFGRAD, CENTER, ANGLE, AU2ANG,
     & CCCCCC, DNS, DOMEGA, DOS, EMN, EMX, FACTOR, GRAD, PSIG, PSIL,
     & PT, PTREL, PTS, QL, QLDS, QTOT, R1, R2, R2MAT, R5, RANG, RHOG,
     & RMN, RVECA, SPN, VTEST, WTRAD, XMAT, XRAD, YLM
       SAVE
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (LMXX=20)
       PARAMETER (LSIZ=(LMXX+1)**2)
       PARAMETER (MAXANG=((2*LMXX+1)*(LMXX+1)))
       PARAMETER (MAXSPH=500)
       PARAMETER (MAXRAD=1000)
       PARAMETER (MXCHG=56)
       CHARACTER*80 LINE
C
       LOGICAL LMKFIL,EXIST
       LOGICAL IUPDAT,ICOUNT
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       COMMON/TMP1/PTS(NSPEED,3),PSIL(MAX_OCC,LSIZ,2)
     &  ,RANG(3,MAXANG),QL(LMXX+2),QTOT(LMXX+2),DOS(LMXX+2)
     &  ,QLDS(LMXX+2,MAX_OCC),YLM(MAXANG,LSIZ)
     &  ,EMN(MAXSPH),EMX(MAXSPH),XRAD(MAXRAD),WTRAD(MAXRAD)
     &  ,CCCCCC(6,MAXSPH),ANGLE(3,MAXANG),DOMEGA(MAXANG)
     &  ,RVECA(3,MX_GRP),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3)
       COMMON/TMP2/PSIG(NMAX,MAX_OCC),RHOG(MAX_PTS,KRHOG,MXSPN)
C
       DIMENSION  CENTER(3,NSPHERE)
C
       DIMENSION ISIZE(3),MSITES(1),SPN(2)
       DIMENSION XMAT(3,3,MAX_IDENT), PT(3), PTREL(3)
       DIMENSION EFGRAD(3,3,MAX_IDENT),RHOFERMI(MAX_IDENT)
       DATA ISIZE/1,3,6/
       DATA AU2ANG/0.529177D0/

C
C CALCULATE CHARGE DENSITIES AT EACH POINT IN SPACE:
C
        IF (DEBUG) THEN
           PRINT *, 'IN NEW MOSS1'
           PRINT *, 'NMSH NSPN', NMSH, NSPN
        END IF
           NXTRA = NSPHERE
           DO J = 1, NSPHERE
              RMSH(1,NMSH+J) = CENTER(1,J)
              RMSH(2,NMSH+J) = CENTER(2,J)
              RMSH(3,NMSH+J) = CENTER(3,J)
           END DO
           NMSH = NMSH+NXTRA
           DO ISPN=1,NSPN
            SPN(ISPN) = 0.0D0
            DO IPTS=1,NMSH
             RHOG(IPTS,ISPN,1)=0.0D0
            END DO
           END DO
C
C POINTS LOOP
C
         LPTS=0
   40    CONTINUE
          MPTS=MIN(NMAX,NMSH-LPTS)
          LPBEG=LPTS
C
C INITIALIZE PSIG 
C
          DO IWF=1,NWF
           DO IPTS=1,MPTS
            PSIG(IPTS,IWF)=0.0D0
           END DO  
          END DO  
          ISHELLA=0
          DO 86 IFNCT=1,NFNCT
           LMX1=LSYMMAX(IFNCT)+1
           DO 84 I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
            DO 82 J_POS=1,M_NUC
             CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                    RVECA,L_NUC,1)
             IF(L_NUC.NE.M_NUC)THEN
              PRINT *,'DECOMP: PROBLEM IN UNRAVEL'
              CALL STOPIT 
             END IF
             LPTS=LPBEG
             DO 80 JPTS=1,MPTS,NSPEED
              NPV=MIN(NSPEED,MPTS-JPTS+1)
              DO LPV=1,NPV
               PTS(LPV,1)=RMSH(1,LPTS+LPV)-RVECA(1,J_POS)
               PTS(LPV,2)=RMSH(2,LPTS+LPV)-RVECA(2,J_POS)
               PTS(LPV,3)=RMSH(3,LPTS+LPV)-RVECA(3,J_POS)
              END DO
C
C GET ORBITS AND DERIVATIVES
C
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
C
C UPDATING ARRAY PSIG
C
              IF (IUPDAT) THEN
               IPTS=JPTS-1
               ILOC=0
               DO 78 LI=1,LMX1
                DO MU=1,ISIZE(LI)
                 DO ICON=1,N_CON(LI,IFNCT)
                  ILOC=ILOC+1
                  IF (ICOUNT(ICON,LI)) THEN
                   DO IWF=1,NWF
                    FACTOR=PSI(ILOC,IWF,1)
                    DO LPV=1,NPV
                     PSIG(IPTS+LPV,IWF)=PSIG(IPTS+LPV,IWF)
     &               +FACTOR*GRAD(LPV,1,MU,ICON,LI)
                    END DO
                   END DO  
                  END IF
                 END DO  
                END DO  
   78          CONTINUE
              END IF
              LPTS=LPTS+NPV
   80        CONTINUE
   82       CONTINUE
   84      CONTINUE
   86     CONTINUE
C
C UPDATE CHARGE DENSITY:
C
           IWF=0
           DO ISPN=1,NSPN
           DO JWF=1,NWFS(ISPN)
           IWF=IWF+1
            DO IPTS=1,MPTS
             RHOG(IPTS+LPBEG,ISPN,1)=
     &       RHOG(IPTS+LPBEG,ISPN,1)+PSIG(IPTS,IWF)**2
            END DO
           END DO
           END DO
C
C CHECK IF ALL POINTS DONE
C
          LPTS=LPBEG+MPTS
          IF(LPTS.GT.MAX_PTS)THEN
           PRINT *,'DECOMP: ERROR: LPTS >',MAX_PTS
           CALL STOPIT
          ELSE
           IF(LPTS.LT.NMSH) GOTO 40
          END IF
  500    CONTINUE
C
C  REMOVE EXTRA POINTS FROM MESH.  FILL UP RHOFERMI ARRAY.  
C
          NMSH = NMSH - NXTRA
          DO J=1,NSPHERE
            IF(NSPN.EQ.1) THEN
              RHOFERMI(J) = 2.0D0*RHOG(NMSH+J,1,1) 
            ELSE
              RHOFERMI(J) = RHOG(NMSH+J,1,1) + RHOG(NMSH+J,2,1)
            END IF
          END DO
C 
C COMPUTE EXPECTATION VALUES
C 
          DO IDT=1,NSPHERE
          DO I=1,3
          DO J=1,3
            EFGRAD(I,J,IDT)=0.0D0
          END DO
          END DO
          END DO
          DO ISPN=1,NSPN
          DO IPTS=1,NMSH
            DO IROT = 1,NGRP
              DO I=1,3
                PT(I)=0.0
                DO J=1,3
                   PT(I) = PT(I) + RMSH(J,IPTS)*RMAT(J,I,IROT)
                END DO
              END DO
C
C  EVALUATE MATRIX ELEMENTS FOR ALL IDENTITY MEMBERS
C
              DO IDT = 1,NSPHERE
                PTREL(1) = PT(1) - CENTER(1,IDT)
                PTREL(2) = PT(2) - CENTER(2,IDT)
                PTREL(3) = PT(3) - CENTER(3,IDT)
                R2 = PTREL(1)**2 + PTREL(2)**2 + PTREL(3)**2
                R1 = SQRT(R2)
                R5 = R2**2.5D0
                DO I=1,3
                DO J=1,3
                  EFGRAD(I,J,IDT) = EFGRAD(I,J,IDT) + 
     &              3.0D0*PTREL(I)*PTREL(J)/R5*2.0D0/FLOAT(NSPN)*
     &              RHOG(IPTS,ISPN,1)*WMSH(IPTS)/FLOAT(NGRP)
                END DO
                END DO
                IF(IDT.EQ.1) VTEST = VTEST + RHOG(IPTS,ISPN,1)/R1
     &                *WMSH(IPTS)/FLOAT(NGRP)*2.0D0/FLOAT(NSPN)
              END DO
              SPN(ISPN)=SPN(ISPN)+RHOG(IPTS,ISPN,1)*WMSH(IPTS)
     &             /FLOAT(NGRP)
            END DO
          END DO
          END DO
         IF (DEBUG) THEN
          PRINT *, 'VTEST', VTEST
          DNS=SPN(1)+SPN(NSPN)
          RMN=SPN(NSPN)-SPN(1) 
          SPN(1)=DNS
          SPN(2)=RMN
         PRINT *, 'CHARGES'
         PRINT *, 'NSPN', NSPN
         PRINT *, DNS, RMN
         END IF 
         DO IDT = 1,NSPHERE
           R2MAT = EFGRAD(1,1,IDT) + EFGRAD(2,2,IDT) + EFGRAD(3,3,IDT)
           R2MAT = R2MAT/3.0D0
           DO I=1,3
              EFGRAD(I,I,IDT) = EFGRAD(I,I,IDT) - R2MAT
           END DO
          END DO
        RETURN
       END
