C UTEP Electronic Structure Lab (2020)
C
C ***************************************************************
C
       SUBROUTINE GETHOLD(NEWIT,IFNCT,JFNCT,AI,AJ,HSUB)
C
C WRITTEN BY MARK R PEDERSON (1998)
C
       use mixpot1,only : POTIN,DVPOT=>POTOUT
       use mesh1,only : rmsh,nmsh
       use debug1
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
       use common3,only : RMAT
       use common5,only : HOLD
       use common6,only : TSPH, LIMSPH, NSPHERES
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NEWIT, IFNCT, JFNCT, I, IALP, IBASE, IBEG, IC, ICALL,
     & ICLC, IDOIT, IEND, IERR, IFNCS, INDEX, IPTS, ISPH, J, JALP,
     & JBASE, JC, JNDEX, KPTS, L, LI, LJ, LMAX, MAXI, MAXJ, MIMAX,
     & MJMAX, MPTS, MUI, MUJ, N1, N3, NEW, NIMAX, NJMAX
       REAL*8 :: AI , AJ, HSUB, AIMJ, AL, ALPHAI, ALPHAJ, ALPIONE,
     & ALPITHR, ALPITWO, ARG, DIFF, DOTL, ENV, EXPI, EXPJ, HTIM1,
     & HTIM2, HTIM3, HTIM4, P, PI4D3, PROD, PRODB, PTOT, RMAX2, RMAXI,
     & RMAXJ, RMINI, RMINJ, RMN, RMX, RSQ, SIZ, SIZT, SS, TIM1, TIM2,
     & TM1, TM2, VLM, VOL, ZERO
       SAVE
       LOGICAL SKIP , FULLSKIPI,FULLSKIPJ,ONSITE
       LOGICAL DOIT
       COMMON/HAMTIM/HTIM1,HTIM2,HTIM3,HTIM4
C       COMMON/TMP1/PISV(13,MAX_PTS)
C       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),DVPOT(MAX_PTS*MXSPN)
       DIMENSION AI(3),AJ(3),IBEG(3),IEND(3),AIMJ(3)
       DIMENSION SS(10,10,MAX_BARE,MAX_BARE)
       DIMENSION P(11,MPBLOCK,2),EXPI(MPBLOCK,MAX_BARE)
     &          ,ENV(MPBLOCK)   ,EXPJ(MPBLOCK,MAX_BARE)
       DIMENSION DOIT(2,MX_SPH),DOTL(MX_SPH)
       DIMENSION PRODB(10,10,MPBLOCK),VOL(MPBLOCK)
       DIMENSION HSUB(MAXUNSYM,MAXUNSYM)
       DIMENSION SKIP(MAX_BARE,MAX_BARE)
       DIMENSION NIMAX(MAX_BARE),NJMAX(MAX_BARE)
       DIMENSION AL(3)
       REAL*8,ALLOCATABLE :: PISV(:,:)
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA ICALL/0/
       DATA ZERO/1.0D-8/ 
C
       IF (DEBUG) PRINT 1010,NEWIT,IFNCT,JFNCT,NSPHERES,NMSH
 1010  FORMAT(' GETHOLD: NIJSP:',5I8)

       CALL GTTIME(TIM1)

       ALLOCATE(PISV(13,MAX_PTS),STAT=IERR)
       IF(IERR.NE.0)WRITE(6,*)'GETHOLD:ERROR ALLOCATING PISV'

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
        IF (N_BARE(IFNCT) .NE. 1) THEN
          ALPITWO=BFALP(N_BARE(IFNCT)-1,IFNCT)-ALPIONE
          ALPITHR=BFALP(N_BARE(IFNCT)-2,IFNCT)-ALPIONE
        END IF
        RMAX2=20.0D0/ALPIONE 
        PI4D3=16.0D0*ATAN(1.0D0)/3.0D0
        CALL GTTIME(TM1)
C
C FIGURE OUT WHICH SPHERES SHOULD BE SKIPPED
C
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
        DO 20 ISPH=1,NSPHERES
         IF (DOIT(1,ISPH)) THEN
          DO 10 KPTS=LIMSPH(1,ISPH),LIMSPH(2,ISPH)
           PISV( 2,KPTS)=RMSH(1,KPTS)-AI(1)
           PISV( 3,KPTS)=RMSH(2,KPTS)-AI(2)
           PISV( 4,KPTS)=RMSH(3,KPTS)-AI(3)
           PISV( 5,KPTS)=PISV(2,KPTS)*PISV(2,KPTS)
           PISV( 6,KPTS)=PISV(3,KPTS)*PISV(3,KPTS)
           PISV( 7,KPTS)=PISV(4,KPTS)*PISV(4,KPTS)
           PISV( 1,KPTS)=PISV(5,KPTS)+PISV(6,KPTS)+PISV(7,KPTS)
           PISV( 8,KPTS)=PISV(2,KPTS)*PISV(3,KPTS)
           PISV( 9,KPTS)=PISV(2,KPTS)*PISV(4,KPTS)
           PISV(10,KPTS)=PISV(3,KPTS)*PISV(4,KPTS)
           PISV(11,KPTS)=DVPOT(KPTS)*EXP(-ALPIONE*PISV(1,KPTS))
           PISV(12,KPTS)=            EXP(-ALPITWO*PISV(1,KPTS))
           PISV(13,KPTS)=            EXP(-ALPITHR*PISV(1,KPTS))
   10     CONTINUE
         END IF
   20   CONTINUE
        CALL GTTIME(TM2)
C
        IF (DEBUG) write(6,*)'NEW:',TM2-TM1,NMSH,ALPIONE,RMAX2
       END IF
       CALL GTTIME(TIM2)
       HTIM1= HTIM1+TIM2-TIM1
       TIM1= TIM2
C
C FIGURE OUT WHICH SPHERES TO SKIP FOR OTHER ATOM
C
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
C
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
C
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
C
       DO I=1,MAXI
        DO J=1,MAXJ
         HSUB(J,I)=0.0D0
        END DO
       END DO
       DO IALP=1,N_BARE(IFNCT)
        ALPHAI=BFALP(IALP,IFNCT)
        DO JALP=1,N_BARE(JFNCT)
         ALPHAJ=BFALP(JALP,JFNCT)
         SKIP(IALP,JALP)=.FALSE.
         ARG=(ALPHAI*ALPHAJ/(ALPHAI+ALPHAJ))
     &      *((AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2)
         IF (ARG.GT.CUTEXP) SKIP(IALP,JALP)=.TRUE.
         DO I=1,10
          DO J=1,10
           SS(I,J,IALP,JALP)=0.0D0
          END DO
         END DO
        END DO
       END DO

       CALL GTTIME(TIM2)
       HTIM2= HTIM2+TIM2-TIM1
       TIM1= TIM2
C
C END ZEROING, CALCULATE OVERLAP OR KINETIC ENERGY MATRIX ELEMENTS
C
       MPTS=0
       DO 100 ISPH=1,NSPHERES
        IF (.NOT.DOIT(2,ISPH)) GOTO 100
        DO 90 KPTS=LIMSPH(1,ISPH),LIMSPH(2,ISPH)
         MPTS=MPTS+1
         P( 1,MPTS,1)=1.0D0
         P( 2,MPTS,1)=PISV( 2,KPTS)
         P( 3,MPTS,1)=PISV( 3,KPTS)   
         P( 4,MPTS,1)=PISV( 4,KPTS)
         P( 5,MPTS,1)=PISV( 5,KPTS)
         P( 6,MPTS,1)=PISV( 6,KPTS)
         P( 7,MPTS,1)=PISV( 7,KPTS)
         P( 8,MPTS,1)=PISV( 8,KPTS)
         P( 9,MPTS,1)=PISV( 9,KPTS)
         P(10,MPTS,1)=PISV(10,KPTS)
         P(11,MPTS,1)=PISV( 1,KPTS)
         VOL(MPTS)=PISV(11,KPTS)
         EXPI(MPTS,N_BARE(IFNCT))=  1.0D0
         IF (N_BARE(IFNCT) .NE. 1) THEN
           EXPI(MPTS,N_BARE(IFNCT)-1)=PISV(12,KPTS)
           EXPI(MPTS,N_BARE(IFNCT)-2)=PISV(13,KPTS)
         END IF
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
           DO IPTS=1,MPTS
            P( 1,IPTS,2)=1.0D0
            P( 2,IPTS,2)=P(2,IPTS,1)+AIMJ(1)
            P( 3,IPTS,2)=P(3,IPTS,1)+AIMJ(2)
            P( 4,IPTS,2)=P(4,IPTS,1)+AIMJ(3)
            P( 5,IPTS,2)=P(2,IPTS,2)*P(2,IPTS,2)
            P( 6,IPTS,2)=P(3,IPTS,2)*P(3,IPTS,2)
            P( 7,IPTS,2)=P(4,IPTS,2)*P(4,IPTS,2)
            P( 8,IPTS,2)=P(2,IPTS,2)*P(3,IPTS,2)
            P( 9,IPTS,2)=P(2,IPTS,2)*P(4,IPTS,2)
            P(10,IPTS,2)=P(3,IPTS,2)*P(4,IPTS,2)
           END DO
           RMINJ=1.0D30
           RMAXJ=0.0D0
           DO IPTS=1,MPTS
            P(11,IPTS,2)=P(5,IPTS,2)+P(6,IPTS,2)+P(7,IPTS,2)
            RMINJ=MIN(RMINJ,P(11,IPTS,2))
            RMAXJ=MAX(RMAXJ,P(11,IPTS,2))
           END DO
          ELSE
           DO IPTS=1,MPTS
            DO J=1,11
             P(J,IPTS,2)=P(J,IPTS,1)
            END DO
           END DO
           RMINJ=RMINI
           RMAXJ=RMAXI
          END IF
C
C NOTE THAT EXPI.NE.EXPJ EVEN FOR ONSITE CALCS
C
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
          DO IPTS=1,MPTS
           DO I=1,MIMAX
            P(I,IPTS,1)=P(I,IPTS,1)*VOL(IPTS)
           END DO
          END DO
          DO IPTS=1,MPTS
           DO J=1,MJMAX
            DO I=1,MIMAX
             PRODB(I,J,IPTS)=P(I,IPTS,1)*P(J,IPTS,2)
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
           IF (ALPHAI*RMINI .GT. CUTEXP) GOTO 50
           DO 40 JALP=1,N_BARE(JFNCT)
            ALPHAJ=BFALP(JALP,JFNCT)
            IF (ALPHAJ*RMINJ .GT. CUTEXP) GOTO 40
            IF (.NOT.SKIP(IALP,JALP)) THEN
C
C UPDATE 10X10 MATRICES
C
             DO IPTS=1,MPTS
              ENV(IPTS)=EXPI(IPTS,IALP)*EXPJ(IPTS,JALP)
             END DO
             DO IPTS=1,MPTS
              IF (PTOT*ENV(IPTS) .GT. ZERO/100.0D0) THEN
               IF (ONSITE) THEN
                IF ((IALP .EQ. N_BARE(IFNCT)) .AND.
     &              (JALP .EQ. N_BARE(IFNCT))) ICLC=ICLC+1
               END IF
               IF ((NJMAX(JALP) .EQ. 10) .AND. 
     &             (NIMAX(IALP) .EQ. 10)) THEN
                DO J=1,10
                 DO I=1,10
                  SS(I,J,IALP,JALP)=SS(I,J,IALP,JALP)
     &                             +PRODB(I,J,IPTS)*ENV(IPTS)
                 END DO
                END DO
               ELSE 
                DO J=1,NJMAX(JALP)
                 DO I=1,NIMAX(IALP)
                  SS(I,J,IALP,JALP)=SS(I,J,IALP,JALP)
     &                             +PRODB(I,J,IPTS)*ENV(IPTS)
                 END DO
                END DO
               END IF
              END IF
             END DO
            END IF
   40       CONTINUE
   50      CONTINUE       
   60     CONTINUE
          MPTS=0
         END IF
   90   CONTINUE
  100  CONTINUE

       CALL GTTIME(TIM2)
       HTIM3= HTIM3+TIM2-TIM1
       TIM1= TIM2
C
       DO 240 IALP=1,N_BARE(IFNCT)
        ALPHAI=BFALP(IALP,IFNCT)
        DO 220 JALP=1,N_BARE(JFNCT)
         ALPHAJ=BFALP(JALP,JFNCT)
C
C        CALL OVMXSF(ALPHAI,ALPHAJ,AI,AJ,SS(1,1,IALP,JALP))
C
         IF (.NOT.SKIP(IALP,JALP)) THEN
          INDEX=0
          DO LI=0,LSYMMAX(IFNCT)
           DO IBASE=1,N_CON(LI+1,IFNCT)
            DO MUI=IBEG(LI+1),IEND(LI+1)
             INDEX=INDEX+1
             JNDEX=0
             DO LJ=0,LSYMMAX(JFNCT)
              DO JBASE=1,N_CON(LJ+1,JFNCT)
               PROD=BFCON(IALP,IBASE,LI+1,IFNCT)
     &             *BFCON(JALP,JBASE,LJ+1,JFNCT)
               DO MUJ=IBEG(LJ+1),IEND(LJ+1)
                JNDEX=JNDEX+1
                HSUB(JNDEX,INDEX)=HSUB(JNDEX,INDEX)
     &                           +PROD*SS(MUI,MUJ,IALP,JALP)
               END DO
              END DO
             END DO
            END DO
           END DO
          END DO
         END IF
  220   CONTINUE
  240  CONTINUE
       CALL GTTIME(TIM2)
       HTIM4= HTIM4+TIM2-TIM1
       TIM1= TIM2
       IF (DEBUG.AND.ONSITE) THEN
        write(6,*)'OVERNUM ICALL, ICLC:',ICALL,ICLC,FLOAT(ICLC)/NMSH
       END IF

       DEALLOCATE(PISV,STAT=IERR)
       IF(IERR.NE.0)WRITE(6,*)'GETHOLD:ERROR DEALLOCATING PISV'

       RETURN
       END 
