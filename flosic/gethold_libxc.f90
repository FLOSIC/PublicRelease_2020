! UTEP Electronic Structure Lab (2020)
!
! ***************************************************************
!
  SUBROUTINE GETHOLD_LIBXC(NEWIT,IFNCT,JFNCT,AI,AJ,HSUB)
!
! WRITTEN BY MARK R PEDERSON (1998)
! modified from gethold 4/2016
!
!       use hstor1,only : HTIM1,HTIM2,HTIM3,HTIM4
       use mixpot1,only : POTIN,DVPOT=>POTOUT
!       use mixpot,only : POTDV,POTIN,DVPOT=>POTOUT
       use XTMP2A,only : mixin
       use mesh1,only : rmsh,nmsh,wmsh
       use debug1
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : HOLD
       use common6,only : TSPH, LIMSPH, NSPHERES
!       INCLUDE 'PARAMS'
       INCLUDE 'PARAMA2'
       SAVE
!The input will be sandwiched with <phi|input|phi>
!       real(8), intent(in) :: mixin(5,NMSH) 
!The output of this subroutine is <phi|input|Grad phi>+<Grad phi|input|phi>

       LOGICAL :: DOIT(2,MX_SPH),SKIP(MAX_BARE,MAX_BARE)
       INTEGER :: IBEG(3),IEND(3)
       INTEGER :: NIMAX(MAX_BARE),NJMAX(MAX_BARE)
       REAL*8  :: DOTL(MX_SPH)
       REAL*8  :: AI(3),AJ(3),AIMJ(3),AL(3),RI(3),RJ(3)
       REAL*8  :: ALPHA2RI(3), ALPHA2RJ(3)
       REAL*8  :: EXPI(MPBLOCK,MAX_BARE), EXPJ(MPBLOCK,MAX_BARE)
       REAL*8  :: ENV(MPBLOCK), VOL(MPBLOCK)
       REAL*8  :: P(11,MPBLOCK,2)
       REAL*8  :: PRODB(10,10,MPBLOCK)
       REAL*8  :: HSUB(MAXUNSYM,MAXUNSYM)
       REAL*8  :: SS(10,10,MAX_BARE,MAX_BARE)
!       REAL*8,ALLOCATABLE :: PISV(:,:)
       real*8  :: VOL1(5,MPBLOCK)
       real*8  :: GRADP(10,3,MPBLOCK,2)
       real*8  :: R( 3,MPBLOCK,2)

       LOGICAL :: FULLSKIPI,FULLSKIPJ,ONSITE
       integer :: I,J,K,L
       integer :: NEW, NEWIT, N1, N3, ICLC, ICALL
       integer :: IDOIT, IERR, IPTR
       integer :: MUI, MUJ, MIMAX, MJMAX, MAXI, MAXJ, LI, LJ, LMAX
       integer :: NIM, NJM
       integer :: IPTS, MPTS, KPTS, ISPH, IFNCS
       integer :: IC, IBASE, IALP, INDEX, IFNCT
       integer :: JC, JBASE, JALP, JNDEX, JFNCT
       real*8  :: ALPHAI, ALPHAJ
       REAL*8  :: RISQ,RJSQ,RARG, RSQ
       real*8  :: R1, R2
       real*8  :: RMN, RMX, RMINI, RMINJ, RMAXI, RMAXJ, RMAX2
       real*8  :: SIZ, SIZT, PTOT, PROD, DIFF, ARG
       real*8  :: VLM, PI4D3, ZERO
       real*8  :: ALPIONE, ALPITWO, ALPITHR
       real*8  :: ALPI2, ALPJ2
       REAL*8  :: TM2, TM1, TIM1
       real*8  :: ENVI, VOLMGGA, GPI(3), GPJ(3), PI, PJ, GPJS

       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA ICALL/0/
       DATA ZERO/1.0D-8/ 
!

         IPTR = (ISPN-1)*NMSH
       IF (DEBUG) PRINT 1010,NEWIT,IFNCT,JFNCT,NSPHERES,NMSH
 1010  FORMAT(' GETHOLD: NIJSP:',5I8)

       CALL GTTIME(TIM1)

!       ALLOCATE(PISV(13,MAX_PTS),STAT=IERR)
!       IF(IERR.NE.0)WRITE(6,*)'GETHOLD:ERROR ALLOCATING PISV'

       IF (NEWIT.EQ.1) THEN
        IFNCS= -1
        AL(1)=1.0D30
        AL(2)=1.0D30
        AL(3)=1.0D30
       END IF
       NEWIT=0
       DIFF=ABS(IFNCS-IFNCT)
       IFNCS=IFNCT
       DO J=1,3
        DIFF=DIFF+ABS(AL(J)-AI(J))
        AL(J)=AI(J)
       END DO
       NEW=0
       IF (DIFF.GT.0.00001) NEW=1
       DIFF=ABS(AI(1)-AJ(1))+ABS(AI(2)-AJ(2))+ABS(AI(3)-AJ(3)) 
       AIMJ(1)=AI(1)-AJ(1)
       AIMJ(2)=AI(2)-AJ(2)
       AIMJ(3)=AI(3)-AJ(3)
       ONSITE=.FALSE.
       IF (DIFF .LT. 0.001D0) THEN
        ONSITE=.TRUE.
        ICALL=ICALL+1     
        ICLC =0
       END IF
       IF (NEW.EQ.1) THEN
        ALPIONE=BFALP(N_BARE(IFNCT),IFNCT)
        ALPITWO=BFALP(N_BARE(IFNCT)-1,IFNCT)-ALPIONE
        ALPITHR=BFALP(N_BARE(IFNCT)-2,IFNCT)-ALPIONE
        RMAX2=20.0D0/ALPIONE 
        PI4D3=16.0D0*ATAN(1.0D0)/3.0D0
        CALL GTTIME(TM1)
!
! FIGURE OUT WHICH SPHERES SHOULD BE SKIPPED
!
        DO ISPH=1,NSPHERES
         RSQ=0.0D0
         DO J=1,3
          RSQ=(AI(J)-TSPH(J,ISPH))**2+RSQ
         END DO
         VLM=PI4D3*TSPH(4,ISPH)**3
         RSQ=SQRT(RSQ)
         RMN=RSQ-TSPH(4,ISPH)
         RMX=RSQ+TSPH(4,ISPH)
         RMX=MAX(RMX,1.0D0)
         DOIT(1,ISPH)=.FALSE.
         IF (RMN .GT. 0.0D0) THEN
          SIZ=EXP(-BFALP(N_BARE(IFNCT),IFNCT)*RMN*RMN)
          SIZ=SIZ*VLM*RMX*RMX
          DOTL(ISPH)=SIZ
          IF (SIZ.GT.ZERO) DOIT(1,ISPH)=.TRUE.
         ELSE 
          DOTL(ISPH)=VLM   
          DOIT(1,ISPH)=.TRUE.
         END IF
        END DO
! PISV removed to save memory -CMD
!
!        DO 20 ISPH=1,NSPHERES
!         IF (DOIT(1,ISPH)) THEN
!          DO 10 KPTS=LIMSPH(1,ISPH),LIMSPH(2,ISPH)
!           PISV( 2,KPTS)=RMSH(1,KPTS)-AI(1)
!           PISV( 3,KPTS)=RMSH(2,KPTS)-AI(2)
!           PISV( 4,KPTS)=RMSH(3,KPTS)-AI(3)
!           PISV( 5,KPTS)=PISV(2,KPTS)*PISV(2,KPTS)
!           PISV( 6,KPTS)=PISV(3,KPTS)*PISV(3,KPTS)
!           PISV( 7,KPTS)=PISV(4,KPTS)*PISV(4,KPTS)
!           PISV( 1,KPTS)=PISV(5,KPTS)+PISV(6,KPTS)+PISV(7,KPTS)
!           PISV( 8,KPTS)=PISV(2,KPTS)*PISV(3,KPTS)
!           PISV( 9,KPTS)=PISV(2,KPTS)*PISV(4,KPTS)
!           PISV(10,KPTS)=PISV(3,KPTS)*PISV(4,KPTS)
!           PISV(11,KPTS)=DVPOT(KPTS)*EXP(-ALPIONE*PISV(1,KPTS))
!           PISV(12,KPTS)=            EXP(-ALPITWO*PISV(1,KPTS))
!           PISV(13,KPTS)=            EXP(-ALPITHR*PISV(1,KPTS))
!   10     CONTINUE
!         END IF
!   20   CONTINUE
        CALL GTTIME(TM2)
!
        IF (DEBUG) write(6,*)'NEW:',TM2-TM1,NMSH,ALPIONE,RMAX2
       END IF
!       CALL GTTIME(TIM2)
!       HTIM1= HTIM1+TIM2-TIM1
!       TIM1= TIM2
!
! FIGURE OUT WHICH SPHERES TO SKIP FOR OTHER ATOM
!
       N1=0
       N3=0
       SIZT=0.0D0
       DO ISPH=1,NSPHERES
        RSQ=0.0D0
        DO J=1,3
         RSQ=(AJ(J)-TSPH(J,ISPH))**2+RSQ
        END DO
        RSQ=SQRT(RSQ)
        RMN=RSQ-TSPH(4,ISPH)
        RMX=RSQ+TSPH(4,ISPH)
        VLM=PI4D3*TSPH(4,ISPH)**3
        RMX=MAX(RMX,1.0D0)
        IF (RMN .GT. 0.0D0) THEN
         SIZ=EXP(-BFALP(N_BARE(JFNCT),JFNCT)*RMN*RMN)
         SIZ=SIZ*VLM*RMX*RMX
        ELSE 
         SIZ=VLM   
        END IF 
        DOIT(2,ISPH)=DOIT(1,ISPH)
        IF (SIZ.LT.ZERO) DOIT(2,ISPH)=.FALSE.
        SIZT=SIZT+SIZ
        SIZ=SIZ*DOTL(ISPH)/VLM
        IF (SIZ .LT. ZERO**1.5D0) DOIT(2,ISPH)=.FALSE.
        IF (DOIT(1,ISPH)) N1=N1+1
        IF (DOIT(2,ISPH)) N3=N3+1
       END DO
       IF (DEBUG) write(6,*)'DOIT:',N1,N3,SIZT,PI4D3,RMX
!
       MAXI=N_CON(1,IFNCT)+3*N_CON(2,IFNCT)+6*N_CON(3,IFNCT)
       MAXJ=N_CON(1,JFNCT)+3*N_CON(2,JFNCT)+6*N_CON(3,JFNCT)
       MIMAX=1
       DO IALP=1,N_BARE(IFNCT)
        LMAX=1
        DO L=2,3
         DO IC=1,N_CON(L,IFNCT)
          IF (ABS(BFCON(IALP,IC,L,IFNCT)) .GT. 0.0D0) THEN
           LMAX=L 
          END IF
         END DO
        END DO
        IF (LMAX.EQ.1) NIMAX(IALP)=1
        IF (LMAX.EQ.2) NIMAX(IALP)=4
        IF (LMAX.EQ.3) NIMAX(IALP)=10
        MIMAX=MAX(NIMAX(IALP),MIMAX)
       END DO
!
       MJMAX=1
       DO JALP=1,N_BARE(JFNCT)
        LMAX=1
        DO L=2,3
         DO JC=1,N_CON(L,JFNCT)
          IF (ABS(BFCON(JALP,JC,L,JFNCT)) .GT. 0.0D0) THEN
           LMAX=L 
          END IF
         END DO
        END DO
        IF (LMAX.EQ.1) NJMAX(JALP)=1
        IF (LMAX.EQ.2) NJMAX(JALP)=4
        IF (LMAX.EQ.3) NJMAX(JALP)=10
        MJMAX=MAX(NJMAX(JALP),MJMAX)
       END DO
!
       DO I=1,MAXI
        DO J=1,MAXJ
         HSUB(J,I)=0.0D0
        END DO
       END DO
       RARG=(AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2
       DO IALP=1,N_BARE(IFNCT)
        ALPHAI=BFALP(IALP,IFNCT)
        DO JALP=1,N_BARE(JFNCT)
         ALPHAJ=BFALP(JALP,JFNCT)
         SKIP(IALP,JALP)=.FALSE.
         ARG=(ALPHAI*ALPHAJ/(ALPHAI+ALPHAJ)) &
     &       * RARG
!     &      *((AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2)
         IF (ARG.GT.CUTEXP) SKIP(IALP,JALP)=.TRUE.
         DO I=1,10
          DO J=1,10
           SS(I,J,IALP,JALP)=0.0D0
          END DO
         END DO
        END DO
       END DO

!       CALL GTTIME(TIM2)
!       HTIM2= HTIM2+TIM2-TIM1
!       TIM1= TIM2

       GRADP(:,:,:,:)=0.0D0 
!
! END ZEROING, CALCULATE OVERLAP OR KINETIC ENERGY MATRIX ELEMENTS
!
       MPTS=0
       DO 100 ISPH=1,NSPHERES
        IF (.NOT.DOIT(2,ISPH)) GOTO 100
        DO 90 KPTS=LIMSPH(1,ISPH),LIMSPH(2,ISPH)
         MPTS=MPTS+1
         RI(1)=RMSH(1,KPTS)-AI(1)
         RI(2)=RMSH(2,KPTS)-AI(2)
         RI(3)=RMSH(3,KPTS)-AI(3)
         RISQ=RI(1)**2 + RI(2)**2 + RI(3)**2
         P( 1,MPTS,1)=1.0D0
         P( 2,MPTS,1)=RI(1)
         P( 3,MPTS,1)=RI(2)
         P( 4,MPTS,1)=RI(3)
         P( 5,MPTS,1)=RI(1)**2
         P( 6,MPTS,1)=RI(2)**2
         P( 7,MPTS,1)=RI(3)**2
         P( 8,MPTS,1)=RI(1)*RI(2)
         P( 9,MPTS,1)=RI(1)*RI(3)
         P(10,MPTS,1)=RI(2)*RI(3)
         P(11,MPTS,1)=RISQ

!         VOL(MPTS)= DVPOT(KPTS)    *EXP(-ALPIONE*P(11,MPTS,1))
   
         EXPI(MPTS,N_BARE(IFNCT))=  1.0D0
         EXPI(MPTS,N_BARE(IFNCT)-1)=EXP(-ALPITWO*RISQ)
         EXPI(MPTS,N_BARE(IFNCT)-2)=EXP(-ALPITHR*RISQ)
!
! 4/8 if ispn=2, POTDV filled with spin 2 V in overnum_libxc
!         VOL1(1,MPTS)=POTDV(1,KPTS) *EXP(-ALPIONE*RISQ)
!         VOL1(2,MPTS)=POTDV(2,KPTS) *EXP(-ALPIONE*RISQ)
!         VOL1(3,MPTS)=POTDV(3,KPTS) *EXP(-ALPIONE*RISQ)
!         VOL1(4,MPTS)=POTDV(4,KPTS) *EXP(-ALPIONE*RISQ)

! mixin has both spins. multiplied by WMSH in overnum_libxc
         IPTR = (ISPN-1)*NMSH
!         VOL1(1,MPTS)=DVPOT(KPTS) *EXP(-ALPIONE*RISQ)
         VOL1(2,MPTS)=mixin(1,KPTS+IPTR) *EXP(-ALPIONE*RISQ)
         VOL1(3,MPTS)=mixin(2,KPTS+IPTR) *EXP(-ALPIONE*RISQ)
         VOL1(4,MPTS)=mixin(3,KPTS+IPTR) *EXP(-ALPIONE*RISQ)
!
         GRADP( 2 ,1,MPTS,1)=1.0D0   *VOL1(2,MPTS)!d/dx (x-Rx)
         GRADP( 3 ,2,MPTS,1)=1.0D0   *VOL1(3,MPTS)!d/dy (y-Ry)
         GRADP( 4 ,3,MPTS,1)=1.0D0   *VOL1(4,MPTS)!d/dz (z-Rz)
                    
         GRADP( 5 ,1,MPTS,1)=2*RI(1) *VOL1(2,MPTS)!d/dx (x-Rx)^2
         GRADP( 6 ,2,MPTS,1)=2*RI(2) *VOL1(3,MPTS)!d/dy (y-Ry)^2
         GRADP( 7 ,3,MPTS,1)=2*RI(3) *VOL1(4,MPTS)!d/dz (z-Rz)^2
                    
         GRADP( 8 ,1,MPTS,1)=RI(2)   *VOL1(2,MPTS)!d/dx (x-Rx)(y-Ry)
         GRADP( 8 ,2,MPTS,1)=RI(1)   *VOL1(3,MPTS)!d/dy (x-Rx)(y-Ry)
         GRADP( 9 ,1,MPTS,1)=RI(3)   *VOL1(2,MPTS)!d/dx (x-Rx)(z-Rz)
         GRADP( 9 ,3,MPTS,1)=RI(1)   *VOL1(4,MPTS)!d/dz (x-Rx)(z-Rz)
         GRADP(10 ,2,MPTS,1)=RI(3)   *VOL1(3,MPTS)!d/dy (y-Ry)(z-Rz)
         GRADP(10 ,3,MPTS,1)=RI(2)   *VOL1(4,MPTS)!d/dz (y-Ry)(z-Rz)
  
         R(1, MPTS, 1)=RI(1) *VOL1(2,MPTS)
         R(2, MPTS, 1)=RI(2) *VOL1(3,MPTS)
         R(3, MPTS, 1)=RI(3) *VOL1(4,MPTS)

         IDOIT=0
         IF (MPTS .EQ. MPBLOCK)       IDOIT=1
         IF (KPTS .EQ. LIMSPH(2,ISPH)) IDOIT=1
         IF (IDOIT.EQ.1) THEN     
          RMINI=1.0D30
          RMAXI=0.0D0
          DO IPTS=1,MPTS
           RMINI=MIN(RMINI,P(11,IPTS,1))
           RMAXI=MAX(RMAXI,P(11,IPTS,1))
          END DO
          IF (.NOT.ONSITE) THEN
           RMINJ=1.0D30
           RMAXJ=0.0D0
           DO IPTS=1,MPTS
            RJ(1)=P(2,IPTS,1)+AIMJ(1) !RJ(1)=(RMSH(1,IPTS)-AI(1))+(AI(1)-AJ(1)) = RMSH(1,IPTS)-AJ(1)
            RJ(2)=P(3,IPTS,1)+AIMJ(2)
            RJ(3)=P(4,IPTS,1)+AIMJ(3)
            P( 1,IPTS,2)=1.0D0
            P( 2,IPTS,2)=RJ(1)
            P( 3,IPTS,2)=RJ(2)
            P( 4,IPTS,2)=RJ(3)
            P( 5,IPTS,2)=RJ(1)**2 
            P( 6,IPTS,2)=RJ(2)**2 
            P( 7,IPTS,2)=RJ(3)**2 
            P( 8,IPTS,2)=RJ(1)*RJ(2) 
            P( 9,IPTS,2)=RJ(1)*RJ(3) 
            P(10,IPTS,2)=RJ(2)*RJ(3) 
!
            P(11,IPTS,2)=P(5,IPTS,2)+P(6,IPTS,2)+P(7,IPTS,2)
            RMINJ=MIN(RMINJ,P(11,IPTS,2))
            RMAXJ=MAX(RMAXJ,P(11,IPTS,2))
!       
            GRADP( 2 ,1,IPTS,2)=1.0D0   *VOL1(2,IPTS) !d/dx (x-Rx)
            GRADP( 3 ,2,IPTS,2)=1.0D0   *VOL1(3,IPTS) !d/dy (y-Ry)
            GRADP( 4 ,3,IPTS,2)=1.0D0   *VOL1(4,IPTS) !d/dz (z-Rz)
                       
            GRADP( 5 ,1,IPTS,2)=2*RJ(1) *VOL1(2,IPTS) !d/dx (x-Rx)^2
            GRADP( 6 ,2,IPTS,2)=2*RJ(2) *VOL1(3,IPTS) !d/dy (y-Ry)^2
            GRADP( 7 ,3,IPTS,2)=2*RJ(3) *VOL1(4,IPTS) !d/dz (z-Rz)^2
                       
            GRADP( 8 ,1,IPTS,2)=RJ(2)   *VOL1(2,IPTS) !d/dx (x-Rx)(y-Ry)
            GRADP( 8 ,2,IPTS,2)=RJ(1)   *VOL1(3,IPTS) !d/dy (x-Rx)(y-Ry)
            GRADP( 9 ,1,IPTS,2)=RJ(3)   *VOL1(2,IPTS) !d/dx (x-Rx)(z-Rz)
            GRADP( 9 ,3,IPTS,2)=RJ(1)   *VOL1(4,IPTS) !d/dz (x-Rx)(z-Rz)
            GRADP(10 ,2,IPTS,2)=RJ(3)   *VOL1(3,IPTS) !d/dy (y-Ry)(z-Rz)
            GRADP(10 ,3,IPTS,2)=RJ(2)   *VOL1(4,IPTS) !d/dz (y-Ry)(z-Rz)

            R(1, IPTS, 2)=RJ(1) *VOL1(2,IPTS)
            R(2, IPTS, 2)=RJ(2) *VOL1(3,IPTS)
            R(3, IPTS, 2)=RJ(3) *VOL1(4,IPTS)
           END DO
          ELSE
!           DO IPTS=1,MPTS
!            DO J=1,11
!             P(J,IPTS,2)=P(J,IPTS,1)
!            END DO
!           END DO
           P(:,:,2)=P(:,:,1)
!
!           DO K=1,3
!            DO IPTS=1,MPTS
!             DO J=1,10
!              GRADP(J,K,IPTS,2)=GRADP(J,K,IPTS,1)
!             END DO
!            END DO
!           END DO
           GRADP(:,:,:,2)=GRADP(:,:,:,1) 
!
           R(:,:,2)=R(:,:,1)
!
           RMINJ=RMINI
           RMAXJ=RMAXI
          END IF
!
! NOTE THAT EXPI.NE.EXPJ EVEN FOR ONSITE CALCS
!
          FULLSKIPI=.FALSE.
          DO IALP=1,N_BARE(IFNCT)-3
           ALPHAI=BFALP(IALP,IFNCT)-ALPIONE
           IF (ALPHAI*RMINI .LT. CUTEXP) THEN
            DO IPTS=1,MPTS
             EXPI(IPTS,IALP)=EXP(-ALPHAI*P(11,IPTS,1))
            END DO
           END IF
          END DO
          FULLSKIPJ=.TRUE.
          DO JALP=1,N_BARE(JFNCT)
           ALPHAJ=BFALP(JALP,JFNCT)
           IF (ALPHAJ*RMINJ .LT. CUTEXP) THEN
            FULLSKIPJ=.FALSE.
            DO IPTS=1,MPTS
             EXPJ(IPTS,JALP)=EXP(-ALPHAJ*P(11,IPTS,2))
            END DO
           END IF
          END DO
          IF (FULLSKIPJ.OR.FULLSKIPI) GOTO 60 
!
! VOL can't be attached to P anymore
! attached to GRADP and R
!           DO IPTS=1,MPTS
!            DO I=1,MIMAX
!             P(I,IPTS,1)=P(I,IPTS,1)*VOL(IPTS)
!            END DO
!           END DO
!
          DO IPTS=1,MPTS
           DO J=1,MJMAX
            DO I=1,MIMAX
!Pi*Pj only needed to calculate PTOT
             PRODB(I,J,IPTS)=P(I,IPTS,1)*P(J,IPTS,2) !*VOL1(1,IPTS)
            END DO
           END DO
          END DO
          PTOT=0.0D0
          DO IPTS=1,MPTS
           DO J=1,MJMAX
            DO I=1,MIMAX
             PTOT=PTOT+ABS(PRODB(I,J,IPTS))
            END DO
           END DO
          END DO
          DO 50 IALP=1,N_BARE(IFNCT)
           ALPHAI=BFALP(IALP,IFNCT)
           ALPI2=2*ALPHAI
           IF (ALPHAI*RMINI .GT. CUTEXP) GOTO 50
           DO 40 JALP=1,N_BARE(JFNCT)
            ALPHAJ=BFALP(JALP,JFNCT)
            ALPJ2=2*ALPHAJ
            IF (ALPHAJ*RMINJ .GT. CUTEXP) GOTO 40
            IF (.NOT.SKIP(IALP,JALP)) THEN
!
! UPDATE 10X10 MATRICES
!
             DO IPTS=1,MPTS
              ENV(IPTS)=EXPI(IPTS,IALP)*EXPJ(IPTS,JALP)
             END DO
             NJM=NJMAX(JALP)
             NIM=NIMAX(IALP)
             DO IPTS=1,MPTS
              IF (PTOT*ENV(IPTS) .GT. ZERO/100.0D0) THEN
               IF (ONSITE) THEN
                IF ((IALP .EQ. N_BARE(IFNCT)) .AND. &
                    (JALP .EQ. N_BARE(IFNCT))) ICLC=ICLC+1
               END IF

               ALPHA2RI(1)=R(1,IPTS,1)*ALPI2
               ALPHA2RI(2)=R(2,IPTS,1)*ALPI2 
               ALPHA2RI(3)=R(3,IPTS,1)*ALPI2 

               ALPHA2RJ(1)=R(1,IPTS,2)*ALPJ2
               ALPHA2RJ(2)=R(2,IPTS,2)*ALPJ2
               ALPHA2RJ(3)=R(3,IPTS,2)*ALPJ2

               ENVI=ENV(IPTS)
!unrolling k=1,3 (gradients) loop speeds up by at least 2x.
                DO J=1,NJM   !NJMAX(JALP)
                 PJ=P(J,IPTS,2) 
                 GPJ(1)=(GRADP(J,1,IPTS,2) -ALPHA2RJ(1)*PJ) !d/dx 
                 GPJ(2)=(GRADP(J,2,IPTS,2) -ALPHA2RJ(2)*PJ) !d/dy
                 GPJ(3)=(GRADP(J,3,IPTS,2) -ALPHA2RJ(3)*PJ) !d/dz

                 GPJS=(GPJ(1)+GPJ(2)+GPJ(3))*ENVI
                 PJ=PJ*ENVI !
                 DO I=1,NIM  !JIMAX(IALP)
                  PI=P(I,IPTS,1)
                  GPI(1)=GRADP(I,1,IPTS,1) -ALPHA2RI(1)*PI
                  GPI(2)=GRADP(I,2,IPTS,1) -ALPHA2RI(2)*PI
                  GPI(3)=GRADP(I,3,IPTS,1) -ALPHA2RI(3)*PI

! Grad(phi_i)*phi_j + Grad(phi_j)*phi_i
                  SS(I,J,IALP,JALP)=SS(I,J,IALP,JALP)  &
                    + (GPI(1)+GPI(2)+GPI(3))*PJ  &
                     + GPJS*PI
!                     + PRODB(I,J,IPTS)*ENVI !standard
                    
! First try: combined terms when possible to optimize
!                    SS(I,J,IALP,JALP)=SS(I,J,IALP,JALP)  &
!                      + ((GRADP(I,IPTS,K,1) -ALPHAI*2*RI(K)*P(I,IPTS,1))* P(J,IPTS,2) &
!                      +  (GRADP(J,IPTS,K,2) -ALPHAJ*2*RJ(K)*P(J,IPTS,2))* P(I,IPTS,1) ) &
!                       *VOL1(K+1,IPTS) *ENVI
!
                 END DO !I
                END DO !J
              END IF !PTOT
             END DO !IPTS
            END IF
   40       CONTINUE
   50      CONTINUE       
   60     CONTINUE
          MPTS=0
         END IF !IDOIT.EQ.1
   90   CONTINUE
  100  CONTINUE

!       CALL GTTIME(TIM2)
!       HTIM3= HTIM3+TIM2-TIM1
!       TIM1= TIM2
!
       DO 240 IALP=1,N_BARE(IFNCT)
        ALPHAI=BFALP(IALP,IFNCT)
        DO 220 JALP=1,N_BARE(JFNCT)
         ALPHAJ=BFALP(JALP,JFNCT)
!
!        CALL OVMXSF(ALPHAI,ALPHAJ,AI,AJ,SS(1,1,IALP,JALP))
!
         IF (.NOT.SKIP(IALP,JALP)) THEN
          INDEX=0
          DO LI=0,LSYMMAX(IFNCT)
           DO IBASE=1,N_CON(LI+1,IFNCT)
!           BI=BFCON(IALP,IBASE,LI+1,IFNCT)
            DO MUI=IBEG(LI+1),IEND(LI+1)
             INDEX=INDEX+1
             JNDEX=0
             DO LJ=0,LSYMMAX(JFNCT)
              DO JBASE=1,N_CON(LJ+1,JFNCT)
               PROD=BFCON(IALP,IBASE,LI+1,IFNCT)  &
                   *BFCON(JALP,JBASE,LJ+1,JFNCT)
               DO MUJ=IBEG(LJ+1),IEND(LJ+1)
                JNDEX=JNDEX+1
                HSUB(JNDEX,INDEX)=HSUB(JNDEX,INDEX)  &
                                 +PROD*SS(MUI,MUJ,IALP,JALP)
               END DO
              END DO
             END DO
            END DO
           END DO
          END DO
         END IF
  220   CONTINUE
  240  CONTINUE
!       CALL GTTIME(TIM2)
!       HTIM4= HTIM4+TIM2-TIM1
!       TIM1= TIM2
       IF (DEBUG.AND.ONSITE) THEN
        write(6,*)'OVERNUM ICALL, ICLC:',ICALL,ICLC,FLOAT(ICLC)/NMSH
       END IF
!       DEALLOCATE(PISV,STAT=IERR)
!       IF(IERR.NE.0)WRITE(6,*)'GETHOLD:ERROR DEALLOCATING PISV'

       RETURN
       END 
