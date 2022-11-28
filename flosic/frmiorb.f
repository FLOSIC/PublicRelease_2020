C UTEP Electronic Structure Lab (2020)
C> @file frmiorb


       SUBROUTINE FRMIORB(NXBG,NXED)
       use pot_dens,only : COULOMB,RHOG
       use mixpot1,only : POTIN,POT=>POTOUT
       use mesh1,only : NMSH,RMSH,WMSH
       use common2,only : IGGA,LSYMMAX,N_CON,N_POS,NFNCT,RIDT,NSPN
     &                   ,E_UP, E_DN
       use common5,only : NWF, NWFS,PSI
       use for_diag1
! lsymmax, n_pos, ridt, n_con, psi, nfrm
!SIC module
       use ORBENG,only : EVALOCC
       use FRM,only    : BFRM,RESULTS,LFRM,DEBDAX
       use HMATSIC,only : OVTM,HMTM
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:45 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NXBG, NXED, I, I_POS, ICON, IERR, IFNCT, IGR, ILOC,
     & IORB, IPTS, ISHDUM, ISHELLA, ISIZE, ISPN, ITRI, IWF, J, J_POS,
     & JBEG, JGR, JLOC, JORB, JPTS, JWF, KGR, KORB, KPTS, KS, L_NUC,
     & LI, LMAX1, LORB, LPTS, LPV, M_NUC, MMSH, MPTS, MU, NDERV, NGRAD,
     & NMAX, NORB, NPV, NTRI
       REAL*8 :: SYMBOL , APT1, CHG, FACTOR, TIME3, TIME4,
     & TIMEGORB
C      INCLUDE 'commons.inc'
!      COMMON/ORBENG/EVALOCC(MAX_OCC)
       PARAMETER (NMAX=MPBLOCK)
       LOGICAL EXIST
C>
C> RETURN:
C> RHOG(IPTS,1, 1)= rho_up
C> RHOG(IPTS,2, 1)= d rho_up/dx
C> RHOG(IPTS,3, 1)= d rho_up/dy
C> RHOG(IPTS,4, 1)= d rho_up/dz
C> RHOG(IPTS,5, 1)= d^2 rho_up/dx^2
C> RHOG(IPTS,6, 1)= d^2 rho_up/dy^2
C> RHOG(IPTS,7, 1)= d^2 rho_up/dz^2
C> RHOG(IPTS,8, 1)= d^2 rho_up/dxdy
C> RHOG(IPTS,9, 1)= d^2 rho_up/dxdz
C> RHOG(IPTS,10,1)= d^2 rho_up/dydz
C> RHOG(IPTS,1, 2)= rho_dn
C> RHOG(IPTS,2, 2)= d rho_dn/dx
C> RHOG(IPTS,3, 2)= d rho_dn/dy
C> RHOG(IPTS,4, 2)= d rho_dn/dz
C> RHOG(IPTS,5, 2)= d^2 rho_dn/dx^2
C> RHOG(IPTS,6, 2)= d^2 rho_dn/dy^2
C> RHOG(IPTS,7, 2)= d^2 rho_dn/dz^2
C> RHOG(IPTS,8, 2)= d^2 rho_dn/dxdy
C> RHOG(IPTS,9, 2)= d^2 rho_dn/dxdz
C> RHOG(IPTS,10,2)= d^2 rho_dn/dydz
C>
!      LOGICAL ICOUNT,IFSIC
       LOGICAL IFSIC
!      COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C      COMMON/LOCORB/TMAT(MAX_OCC,MAX_OCC,MXSPN),MORB(2),IRBSIC,ZSIC
!      COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
!      COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
!    &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
!    &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
!      COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!    &   LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)
!      DIMENSION RORB(3,MAX_OCC),CHAM(NDH,NDH,2)
C      DIMENSION OVTM(NDH,NDH,2),HMTM(NDH,NDH,2)
!      COMMON/HMATSIC/ OVTM(MAX_OCC,MAX_OCC,2),HMTM(MAX_OCC,MAX_OCC,2)
!YY copied from common.inc
!      COMMON/FERMIONS/E_UP,E_DN,WFFILE
!      COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
!    &  EVAL(NDH),SC1(NDH),SC2(NDH)
!      COMMON/LOCORB/TMAT(MAX_OCC,MAX_OCC,MXSPN),MORB(2),ZSIC,IRBSIC
!      COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),ZTZL,JJJJJJ
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
       REAL*8,ALLOCATABLE :: RORB(:,:),CHAM(:,:,:)
       REAL*8,ALLOCATABLE :: PSIG(:,:,:),PTS(:,:)
     &                      ,GRAD(:,:,:,:,:),RVECA(:,:)
       LOGICAL,ALLOCATABLE :: ICOUNT(:,:)

       LOGICAL LGGA,IUPDAT,SCPT
       DIMENSION ISIZE(3)
       DATA ISIZE/1,3,6/

       MORB(1)=NINT(E_UP)
       MORB(2)=NINT(E_DN)
       PRINT*,'MORB(1), MORB(2):',MORB(1),MORB(2)
C      INQUIRE(FILE='SICPOT',EXIST=EXIST)
C      SCPT=EXIST
C      IF(SCPT)THEN
C       OPEN(91,FILE='SICPOT',FORM='UNFORMATTED')
C      END IF
C
C      WRITE(INITOUT,*)'EXECUTING FRMIORB',IRANK,INITOUT
       TIMEGORB=0.0D0
       CALL GTTIME(APT1)
!YY allocate for_diag1
!      if(.not.allocated(HAM)) ALLOCATE(HAM(NDH,NDH))
!      if(.not.allocated(OVER)) ALLOCATE(OVER(NDH,NDH))
!      if(.not.allocated(EVAL)) ALLOCATE(EVAL(NDH))
!      if(.not.allocated(SC1)) ALLOCATE(SC1(NDH))
!      if(.not.allocated(FILO))ALLOCATE(FILO(NDH,NDH))

       ALLOCATE(RORB(3,MAX_OCC))
       ALLOCATE(CHAM(NDH,NDH,2))

       ALLOCATE(PSIG(NMAX,10,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING PSIG'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING PTS'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING GRAD'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING RVECA'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR ALLOCATING ICOUNT'
       LGGA= .FALSE.
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) THEN
        LGGA= .TRUE.
        NGRAD=10
       END IF
C      WRITE(INITOUT,*)'IRANK LGGA',IRANK,LGGA
       NORB=0
       DO ISPN=1,NSPN
        DO JORB=1,MORB(ISPN)
         NORB=NORB+1
         DO J=1,3
          RORB(J,NORB)=BFRM(J,JORB,ISPN)
         END DO

C        READ(10,*,END=6)(RORB(J,NORB),J=1,3)
        END DO
       END DO
6      CONTINUE
       NTRI=2
       IF(.NOT.SCPT)NTRI=1
       DO 1000 ITRI=NXBG,NXED
        IF(ITRI.EQ.2)THEN
         DO ISPN=1,NSPN
          DO I=1,MAX_OCC
           DO J=1,MAX_OCC
            OVTM(J,I,ISPN)=0.0D0
            HMTM(J,I,ISPN)=0.0D0
           END DO
          END DO
         END DO
        END IF
        IF(ITRI.EQ.1)MMSH=NORB
        IF(ITRI.EQ.2)MMSH=NMSH
        IF(ITRI.EQ.1)THEN
         DO ISPN=1,NSPN
          DO IORB=1,NORB
           DO IWF =1,NWF
            TMAT(IORB,IWF,ISPN)=0.0D0
           END DO
          END DO
         END DO
        END IF


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
           PRINT *,'APTSLV: MAX_OCC MUST BE AT LEAST:',NWF
           CALL STOPIT
          END IF
C
C FOR ALL EQUIVALENT POSITIONS OF THIS ATOM
C
          DO 82 J_POS=1,M_NUC
C
C UNSYMMETRIZE
C
           IF(ITRI.EQ.1)THEN
            CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                   RVECA,L_NUC,1)
           ELSE
            CALL LORAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                    RVECA,L_NUC,1)
           END IF
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
             IF(ITRI.EQ.1)THEN
              PTS(LPV,1)=RORB(1,LPTS+KPTS)-RVECA(1,J_POS)
              PTS(LPV,2)=RORB(2,LPTS+KPTS)-RVECA(2,J_POS)
              PTS(LPV,3)=RORB(3,LPTS+KPTS)-RVECA(3,J_POS)
             ELSE
              PTS(LPV,1)=RMSH(1,LPTS+KPTS)-RVECA(1,J_POS)
              PTS(LPV,2)=RMSH(2,LPTS+KPTS)-RVECA(2,J_POS)
              PTS(LPV,3)=RMSH(3,LPTS+KPTS)-RVECA(3,J_POS)
             END IF
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
C            PRINT *,'IRANK,LMAX1',IRANK,LMAX1
             DO 78 LI=1,LMAX1
              DO MU=1,ISIZE(LI)
               DO ICON=1,N_CON(LI,IFNCT)
                ILOC=ILOC+1
                IF (ICOUNT(ICON,LI)) THEN
                 DO IWF=1,NWF
                  FACTOR=PSI(ILOC,IWF,1)
                  if(abs(FACTOR) .GT. 1.0d-10)then
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
        IF(ITRI.EQ.2)THEN
         DO ISPN=1,NSPN
          IF(ISPN.EQ.1)KS=0
          IF(ISPN.EQ.2)KS=MORB(1)
          DO IORB=1,MORB(ISPN)
C          REWIND(91)
C          DO IRD=1,KS+IORB
C           READ(91)(POT(IMSH),IMSH=1,NMSH)
C          END DO
           DO JORB=IORB,MORB(ISPN)
            DO IPTS=1,MPTS
             OVTM(JORB,IORB,ISPN)=OVTM(JORB,IORB,ISPN)
     &      -PSIG(IPTS,1,IORB+KS)*PSIG(IPTS,1,JORB+KS)*WMSH(IPTS+LPTS)
            END DO
            OVTM(IORB,JORB,ISPN)=OVTM(JORB,IORB,ISPN)
            DO IPTS=1,MPTS
             HMTM(JORB,IORB,ISPN)=HMTM(JORB,IORB,ISPN)
     &      -PSIG(IPTS,1,IORB+KS)*PSIG(IPTS,1,JORB+KS)
     &      *WMSH(IPTS+LPTS)*POT(IPTS+LPTS)
            END DO
           END DO
           DO JORB=IORB+1,MORB(ISPN)
C           READ(91)(POT(IMSH),IMSH=1,NMSH)
            DO IPTS=1,MPTS
             HMTM(IORB,JORB,ISPN)=HMTM(IORB,JORB,ISPN)
     &      -PSIG(IPTS,1,IORB+KS)*PSIG(IPTS,1,JORB+KS)
     &      *WMSH(IPTS+LPTS)*POT(IPTS+LPTS)
            END DO
           END DO
          END DO
         END DO
        END IF
C
C UPDATING RHOG, START WITH DENSITY
C
        DO ISPN=1,NSPN
         JBEG= (ISPN-1)*NWFS(1)
         DO JWF=1,NWFS(ISPN)
          JLOC=JWF+JBEG
          DO IPTS=1,MPTS
           RHOG(LPTS+IPTS,1,ISPN)=RHOG(LPTS+IPTS,1,ISPN)
     &    +PSIG(IPTS,1,JLOC)**2
          END DO
         END DO
        END DO
        IF(ITRI.EQ.1)THEN
C CALCULATE NON-ORTHOGONAL TRANSFORMATION
         PRINT*,'RHOG:',(RHOG(IPTS,1,1),IPTS=1,MORB(1)+MORB(NSPN))
         PRINT*,'MPTS:',MPTS
         DO ISPN=1,NSPN
          JBEG=(ISPN-1)*NWFS(1)
          DO JWF=1,NWFS(ISPN)
           JLOC=JWF+JBEG
           DO IPTS=1,MPTS
            TMAT(LPTS+IPTS,JWF,ISPN)=
     &      PSIG(IPTS,1,JLOC)/SQRT(RHOG(LPTS+IPTS,1,ISPN))
           END DO
          END DO
         END DO
        END IF
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
     &      +2*PSIG(IPTS,1,JLOC)*PSIG(IPTS,IGR,JLOC)
            END DO
           END DO
C
C SECOND DERIVATIVES (XX,YY,ZZ)
C
           DO IGR=5,7
            JGR=IGR-3
            DO IPTS=1,MPTS
             RHOG(LPTS+IPTS,IGR,ISPN)=RHOG(LPTS+IPTS,IGR,ISPN)
     &        +2*(PSIG(IPTS,JGR,JLOC)**2
     &           +PSIG(IPTS,IGR,JLOC)*PSIG(IPTS,1,JLOC))
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
        IF (LPTS .LT. MMSH) GOTO 10
        CONTINUE
        IF(ITRI.EQ.1)THEN
C        PRINT*,'TMAT:'
C        DO ISPN=1,NSPN
C         PRINT*,'ISPN:',ISPN
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(TMAT(IORB,JORB,ISPN),JORB=1,MORB(ISPN))
C         END DO
C        END DO
         DO ISPN=1,NSPN
          IF(ISPN.EQ.2)THEN
           DO IORB=1,MORB(ISPN)
            DO JWF=1,NWFS(ISPN)
             TMAT(IORB,JWF,NSPN)=TMAT(IORB+MORB(1),JWF,NSPN)
            END DO
           END DO
          END IF
C         PRINT*,'ORBITAL TRANSFORMATION FOR SPIN:',ISPN,MORB(ISPN)
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(TMAT(IORB,JWF,ISPN),JWF=1,NWFS(ISPN))
C         END DO
C calculate overlap  matrix:
          PRINT*,"NON-ORTHONORMAL OVERLAP MATRIX FOR SPIN:",ISPN
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            OVER(IORB,JORB)=0.0D0
            DO JWF=1,NWFS(ISPN)
             OVER(IORB,JORB)=OVER(IORB,JORB)
     &      +TMAT(IORB,JWF,ISPN)*TMAT(JORB,JWF,ISPN)
            END DO
           END DO
C          PRINT 510,(OVER(IORB,JORB),JORB=1,MORB(ISPN))
          END DO
          PRINT*,'CALLING LOWDEN'
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            HAM(JORB,IORB)=0.0D0
           END DO
           HAM(JORB,JORB)=1.0D0
          END DO
          IF(MORB(ISPN).NE.NFRM(ISPN))THEN
           PRINT*,'> SIC ERROR, MORB .NE. NFRM                 <'
           PRINT*,'> Check if you have fract. electron numbers <'
           PRINT*,'> in your SYMBOL / INPUT files              <'
           PRINT*,'> (we need to fix this somehow)             <'
           CALL STOPIT
          END IF
          CALL LOWSIC(ISPN,NDH,MORB(ISPN),OVER,HAM,FILO,EVAL,SC1)
          CALL SWAPFGI(ISPN,MORB(ISPN),0)
C         IF(ISPN.EQ.1)OPEN(47,FILE='FILO1',FORM='UNFORMATTED')
C         IF(ISPN.EQ.2)OPEN(47,FILE='FILO2',FORM='UNFORMATTED')
C         REWIND(47)
C         WRITE(47)((FILO(JFM,IFM),JFM=1,NFRM(ISPN)),IFM=1,NFRM(ISPN))
C         WRITE(47)( EVAL(JFM)    ,JFM=1,NFRM(ISPN))
C         CLOSE(47)
          PRINT*,'EVALOCC:'
          PRINT 510,(EVAL(IORB),IORB=1,MORB(ISPN))
C         PRINT*,'TRANFORMATION FOR SPIN:',ISPN
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(HAM(JORB,IORB),JORB=1,MORB(ISPN))
C         END DO
C NOW ROTATE TMAT:
          DO IORB=1,MORB(ISPN)
           DO JWF=1,NWFS(ISPN)
            OVER(IORB,JWF)=0.0D0
           END DO
          END DO
          DO IORB=1,MORB(ISPN)
           DO KORB=1,MORB(ISPN)
            DO JWF=1,NWFS(ISPN)
             OVER(IORB,JWF)=OVER(IORB,JWF)
     &      +HAM(KORB,IORB)*TMAT(KORB,JWF,ISPN)
            END DO
           END DO
          END DO
          DO IORB=1,MORB(ISPN)
           DO JWF=1,NWFS(ISPN)
            TMAT(IORB,JWF,ISPN)=OVER(IORB,JWF)
           END DO
          END DO
          PRINT*,'ORTHOGONAL TRANSFORMATION'
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(TMAT(IORB,JWF,ISPN),JWF=1,NWFS(ISPN))
C         END DO
          PRINT*,"ORTHONORMAL OVERLAP MATRIX FOR SPIN:",ISPN
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            OVER(IORB,JORB)=0.0D0
            CHAM(IORB,JORB,ISPN)=0.0D0
            DO JWF=1,NWFS(ISPN)
             JLOC=JWF
             IF(ISPN.EQ.2)JLOC=JWF+NWFS(1)
             OVER(IORB,JORB)=OVER(IORB,JORB)
     &       +TMAT(IORB,JWF,ISPN)*TMAT(JORB,JWF,ISPN)
             CHAM(IORB,JORB,ISPN)=CHAM(IORB,JORB,ISPN)
     &       +TMAT(IORB,JWF,ISPN)*TMAT(JORB,JWF,ISPN)*EVALOCC(JLOC)
            END DO
           END DO
C          PRINT 510,(OVER(IORB,JORB),JORB=1,MORB(ISPN))
          END DO
!         PRINT*,"ORTHONORMAL HAMILT FOR SPIN:",ISPN
!         DO IORB=1,MORB(ISPN)
!          PRINT 510,(CHAM(IORB,JORB,ISPN),JORB=1,MORB(ISPN))
!         END DO
         END DO
        ELSE
C CHECK TO SEE IF THE CHARGE INTEGRATES CORRECTLY:
         DO ISPN=1,NSPN
          CHG=0.0d0
          DO IPTS=1,MMSH
           CHG=CHG+RHOG(IPTS,1,ISPN)*WMSH(IPTS)
          END DO
          PRINT*,ISPN,CHG,',SPIN & CHG','ITRI:',ITRI,MORB(ISPN)
C OVERLAP MATRIX:
          PRINT*,'OVERLAP:'
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(OVTM(IORB,JORB,ISPN),JORB=1,MORB(ISPN))
C         END DO
C SIC      MATRIX:
          PRINT*,'SIC MATRIX:'
          DO IORB=1,MORB(ISPN)
           PRINT 510,(HMTM(IORB,JORB,ISPN),JORB=1,MORB(ISPN))
          END DO
          PRINT*,'LSD MATRIX:'
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            HAM (JORB,IORB)=CHAM(JORB,IORB,ISPN)
            OVER(JORB,IORB)=0.0D0
           END DO
           OVER(IORB,IORB)=1.0D0
C          PRINT 510,(CHAM(JORB,IORB,ISPN),JORB=1,MORB(ISPN))
          END DO
          CALL DIAGGE(NDH,MORB(ISPN),HAM,OVER,EVALOCC,SC1,1)
C NOW TRANSFORM THE SIC MATRIX TO THE CANONICAL REPRESENTATION
          DO IORB=1   ,MORB(ISPN)
           DO JORB=IORB,MORB(ISPN)
            HMTM(JORB,IORB,ISPN)=
     &     (HMTM(JORB,IORB,ISPN)+HMTM(IORB,JORB,ISPN))/2.
            HAM(JORB,IORB)=
     &      HMTM(JORB,IORB,ISPN)+CHAM(IORB,JORB,ISPN)
            HAM(IORB,JORB)=HAM(JORB,IORB)
           END DO
          END DO
C
C CONVERT THIS TO AN N^3 LOOP:
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            OVER(JORB,IORB)=0.0D0
            DO KORB=1,MORB(ISPN)
             DO LORB=1,MORB(ISPN)
              OVER(JORB,IORB)=OVER(JORB,IORB)+
     &        TMAT(KORB,JORB,ISPN)*TMAT(LORB,IORB,ISPN)*HAM(KORB,LORB)
             END DO
            END DO
           END DO
          END DO
          PRINT *,'EIGENVALUES'
          PRINT 510,(EVALOCC(I),I=1,MORB(ISPN))
          PRINT*,'CANONICAL SIC-LSD MATRIX FOR SPIN:',ISPN
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            HAM(JORB,IORB)=OVER(JORB,IORB)
           END DO
           OVER(IORB,IORB)=OVER(IORB,IORB)-EVALOCC(IORB) !SHOULD HAVE OFFSET
          END DO
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(HAM(JORB,IORB),JORB=1,MORB(ISPN))
C         END DO
          PRINT*,'CANONICAL SIC SCISSOR OPERATOR:',ISPN
          IF(ISPN.EQ.1)OPEN(95,FILE='SCISSOR1')
          IF(ISPN.EQ.2)OPEN(95,FILE='SCISSOR2')
          REWIND(95)
          WRITE(95,*)MORB(ISPN)
          DO IORB=1,MORB(ISPN)
           WRITE(95,510)(OVER(JORB,IORB),JORB=1,MORB(ISPN))
          END DO
          CLOSE(95)
C         DO IORB=1,MORB(ISPN)
C          PRINT 510,(OVER(JORB,IORB),JORB=1,MORB(ISPN))
C         END DO
          DO IORB=1,MORB(ISPN)
           DO JORB=1,MORB(ISPN)
            OVER(JORB,IORB)=0.0D0
           END DO
           OVER(IORB,IORB)=1.0D0
          END DO
          CALL DIAGGE(NDH,MORB(ISPN),HAM,OVER,EVAL,SC1,1)
          IF(MORB(ISPN).LE.10)THEN
           PRINT 511,(EVAL(I),I=1,MORB(ISPN))
          ELSE
           PRINT 511,(EVAL(I),I=1,10)
           PRINT 510,(EVAL(I),I=11,MORB(ISPN))
          END IF
!         DO IORB=1,MORB(ISPN)
!          PRINT*,'EIGENVECTOR:',IORB
C          PRINT 510,(HAM(J,IORB),J=1,MORB(ISPN))
!         END DO
         END DO
        END IF
 1000  CONTINUE
C      CLOSE(91)
 500   FORMAT(I5,3F15.6,' Fermi Orbital & Position')
 510   FORMAT(20F12.6)
 511   FORMAT('SIC-LSD EIGENVALUES:',20F12.6)



       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING ICOUNT'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING RVECA'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING GRAD'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING PTS'
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DENSOLD:ERROR DEALLOCATING PSIG'
       DEALLOCATE(PTS,STAT=IERR)

       DEALLOCATE(CHAM)
       DEALLOCATE(RORB)

       RETURN !returning from Fermiorb'
       END
