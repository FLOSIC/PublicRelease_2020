C UTEP Electronic Structure Lab (2020)
C
C *****************************************************
C *****************************************************
C DENSOLD BASED ON OLD VERSION OF APOTNL BY M. PEDERSON AND D. POREZAG
C ATTENTION: FIRST TWO ARRAYS OF COMMON BLOCK TMP1 MUST BE IDENTICAL IN 
C DENSOLD AND APOTNL SINCE THEY ARE USED TO PASS DENSITY AND COULOMB POT
C
       SUBROUTINE DENSOLD_old(TIMEGORB)
       use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT, IGGA, ISPN, NSPN
       use common5,only : PSI, NWF, NWFS
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:48 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I_POS, ICON, IFNCT, IGR, ILOC, IPTS, ISHDUM, ISHELLA,
     & ISIZE, IWF, J_POS, JBEG, JGR, JLOC, JPTS, JWF, KGR, KPTS, L_NUC,
     & LI, LMAX1, LPTS, LPV, M_NUC, MPTS, MU, NDERV, NGRAD, NMAX, NMSH,
     & NPV, NSPN
       REAL*8 :: SYMBOL , TIMEGORB, APT1, COULOMB, FACTOR, GRAD, PSIG,
     & PTS, RHOG, RVECA, TIME3, TIME4
       SAVE
       PARAMETER (NMAX=MPBLOCK)
C
C RETURN:
C RHOG(IPTS,1, 1)= rho_up   
C RHOG(IPTS,2, 1)= d rho_up/dx
C RHOG(IPTS,3, 1)= d rho_up/dy
C RHOG(IPTS,4, 1)= d rho_up/dz
C RHOG(IPTS,5, 1)= d^2 rho_up/dx^2
C RHOG(IPTS,6, 1)= d^2 rho_up/dy^2
C RHOG(IPTS,7, 1)= d^2 rho_up/dz^2
C RHOG(IPTS,8, 1)= d^2 rho_up/dxdy
C RHOG(IPTS,9, 1)= d^2 rho_up/dxdz
C RHOG(IPTS,10,1)= d^2 rho_up/dydz
C RHOG(IPTS,1, 2)= rho_dn   
C RHOG(IPTS,2, 2)= d rho_dn/dx
C RHOG(IPTS,3, 2)= d rho_dn/dy
C RHOG(IPTS,4, 2)= d rho_dn/dz
C RHOG(IPTS,5, 2)= d^2 rho_dn/dx^2
C RHOG(IPTS,6, 2)= d^2 rho_dn/dy^2
C RHOG(IPTS,7, 2)= d^2 rho_dn/dz^2
C RHOG(IPTS,8, 2)= d^2 rho_dn/dxdy
C RHOG(IPTS,9, 2)= d^2 rho_dn/dxdz
C RHOG(IPTS,10,2)= d^2 rho_dn/dydz
C
       LOGICAL ICOUNT
       COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,KRHOG,MXSPN)
       COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       LOGICAL LGGA,IUPDAT
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/
	double precision ftime
	external ftime
C
       TIMEGORB=0.0D0
       CALL GTTIME(APT1)
       LGGA= .FALSE.
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) THEN
        LGGA= .TRUE.
        NGRAD=10
       END IF
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
C
C INITIALIZE PSIG AND RHOB
C
        DO IWF=1,NWF
         DO IGR=1,NGRAD
          DO IPTS=1,MPTS
           PSIG(IPTS,IGR,IWF)=0.0D0
          END DO
         END DO  
        END DO  
        DO ISPN=1,NSPN
         DO IGR=1,NGRAD
          DO IPTS=1,MPTS
           RHOG(LPTS+IPTS,IGR,ISPN)=0.0D0
          END DO
         END DO  
        END DO  
        ISHELLA=0
C
C FOR ALL CENTER TYPES
C
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
           write(6,*)'APTSLV: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
          DO 82 J_POS=1,M_NUC
C
C UNSYMMETRIZE 
C
           CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                  RVECA,L_NUC,1)
           IF(L_NUC.NE.M_NUC)THEN
            write(6,*)'APTSLV: PROBLEM IN UNRAVEL'
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
C
C GET ORBITS AND DERIVATIVES
C
            NDERV=0
            IF (LGGA) NDERV=2
            CALL GTTIME(TIME3)
            CALL GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
            CALL GTTIME(TIME4)
            TIMEGORB=TIMEGORB+TIME4-TIME3
C
C UPDATING ARRAY PSIG
C
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
                  DO IGR=1,NGRAD
                   DO LPV=1,NPV
                    PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF)
     &              +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
                   END DO
                  END DO  
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
C
C UPDATING RHOG, START WITH DENSITY 
C
        DO ISPN=1,NSPN
         JBEG= (ISPN-1)*NWFS(1) 
         DO JWF=1,NWFS(ISPN)
          JLOC=JWF+JBEG
          DO IPTS=1,MPTS
           RHOG(LPTS+IPTS,1,ISPN)=RHOG(LPTS+IPTS,1,ISPN)
     &     +PSIG(IPTS,1,JLOC)**2
          END DO
         END DO
        END DO
C
C UPDATE DERIVATIVES IF GGA CALCULATION
C         
        IF (LGGA) THEN
         DO 96 ISPN=1,NSPN
          JBEG= (ISPN-1)*NWFS(1)
          DO 94 JWF=1,NWFS(ISPN)
           JLOC=JWF+JBEG
C
C GRADIENT 
C
           DO IGR=2,4
            DO IPTS=1,MPTS
             RHOG(LPTS+IPTS,IGR,ISPN)=RHOG(LPTS+IPTS,IGR,ISPN)
     &       +2*PSIG(IPTS,1,JLOC)*PSIG(IPTS,IGR,JLOC)
            END DO
           END DO
C
C SECOND DERIVATIVES (XX,YY,ZZ)
C
           DO IGR=5,7
            JGR=IGR-3
            DO IPTS=1,MPTS
             RHOG(LPTS+IPTS,IGR,ISPN)=RHOG(LPTS+IPTS,IGR,ISPN)
     &       +2*(PSIG(IPTS,JGR,JLOC)**2
     &          +PSIG(IPTS,IGR,JLOC)*PSIG(IPTS,1,JLOC))
            END DO
           END DO
C
C SECOND DERIVATIVES (XY,XZ,YZ)
C
           DO IGR=2,3
            DO JGR=IGR+1,4
             KGR=IGR+JGR+3
             DO IPTS=1,MPTS
              RHOG(LPTS+IPTS,KGR,ISPN)=RHOG(LPTS+IPTS,KGR,ISPN)
     &        +2*(PSIG(IPTS,IGR,JLOC)*PSIG(IPTS,JGR,JLOC)
     &           +PSIG(IPTS,KGR,JLOC)*PSIG(IPTS,1,JLOC))
             END DO
            END DO
           END DO
   94     CONTINUE
   96    CONTINUE
        END IF
        LPTS=LPTS+MPTS
        IF (LPTS .LT. NMSH) GOTO 10
       CONTINUE
       RETURN
       END
