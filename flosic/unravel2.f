C UTEP Electronic Structure Lab (2020)
C> @file unravel2.ftn
C> @brief YY. unravel for SIC
C> @param[in] MD=0  FIND OCCUPIED ORBITAL KORB OF SPIN KSPN
C>            MD=1  FIND ALL OCCUPIED ORBITALS.
C>            MD=2  FIND ALL OCCUPIED AND UNOCCUPIED ORBITALS
       SUBROUTINE UNRAVEL2
     &(MD,KSPN,KORB,IFNCT,ISHELLA,I_SITE,RVEC,RVECI,N_NUC,ILOC)
C> @author ORIGINALLY WRITTEN BY MARK R PEDERSON (1985)
       use debug1
       use common2,only : N_CON, LSYMMAX, ISPN, NSPN
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF, NWFS
       use common7,only : T1UNRV, T2UNRV
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC, INDBEG
     &                   ,NS_TOT
!SIC module
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MD, KSPN, KORB, IFNCT, ISHELLA, I_SITE, N_NUC, ILOC,
     & I_LOCAL, I_SALC, IBASE, ICALL1, ICALL2, IL, IMS, IND_SALC,
     & INDEX, IOCC, IORB, IORBB, IORBE, IQ, IQ_BEG, IROW, ISHELL,
     & ISHELLV, ITOT, IWF, J_LOCAL, JORB, JWF, K_REP, K_ROW, KSALC, LI,
     & LSPN, LSPNB, LSPNE, MLOCAL, MO, MSTOT, MU, NDEG
       REAL(8) :: RVEC , RVECI, ADD1, ADD2, ADD3, ADDTMAT, PHILOC,
     & TIMER1, TIMER2, TMATMAX
!       INCLUDE 'commons.inc'
       SAVE
!       COMMON/LOCORB/TMAT(MAX_OCC,MAX_OCC,MXSPN),MORB(2),ZSIC,IRBSIC
!       COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),ZTZL,JJJJJJ

       DIMENSION NDEG(3),IND_SALC(ISMAX,MAX_CON,3)
       DIMENSION RVECI(3,MX_GRP),RVEC(3)
       DIMENSION ISHELLV(2)
       DIMENSION PHILOC(MAX_OCC)
       DIMENSION MSTOT(2)
       DATA NDEG/1,3,6/
       DATA ICALL1,ICALL2,ISHELL/0,0,0/
C MD=0  FIND OCCUPIED ORBITAL KORB OF SPIN KSPN
C MD=1  FIND ALL OCCUPIED ORBITALS.
C MD=2  FIND ALL OCCUPIED AND UNOCCUPIED ORBITALS
C
!        ADDTMAT=0.0D0     
!        TMATMAX=0.0D0
!        DO JWF =1,MORB(KSPN)
!          DO IORB=1,NWFS(KSPN)
!            ADDTMAT=ADDTMAT+ABS(TMAT(IORB,JWF,KSPN))
!            TMATMAX=MAX(TMAT(IORB,JWF,KSPN),TMATMAX)
!          END DO
!        END DO
!        IF(ADDTMAT.LE.0.001.OR.MD.EQ.0)THEN
!C       PRINT 40,((ABS(TMAT(I,J,KSPN)),I=1,MORB(KSPN)),J=1,MORB(KSPN))
!C       PRINT*,'TMATMAX,IR,KSP,MD:',TMATMAX,IRANK,KSPN,MD
!        IF(ADDTMAT.LE.0.001)CALL STOPIT
! 40     FORMAT(10F10.2)
!        END IF

C       if(kspn.eq.2) then
C         print*,'KSPN!!!'
C       END IF


C Note: MD=0 puts result in PSI(:,1,ILOC)

       IF(I_SITE.EQ.1)THEN
        ICALL1=ICALL1+1
        CALL GTTIME(TIMER1)
        CALL OBINFO(1,RVEC,RVECI,N_NUC,ISHELLV(ILOC))
!        if(MD==2)return !testing timing of next section. reduced. time in gsmat
        CALL GSMAT(ISHELLV(ILOC),ILOC)
!        if(MD==2)return !testing timing of next section. same. 
        CALL GTTIME(TIMER2)
        T1UNRV=T1UNRV+TIMER2-TIMER1
        IF (DEBUG.AND.(1000*(ICALL1/1000).EQ.ICALL1)) THEN
         PRINT*,'WASTED1=',T1UNRV,ISHELL
         PRINT*,'ICALL2,AVERAGE:',ICALL1,T1UNRV/ICALL2
        END IF
       END IF
       CALL GTTIME(TIMER1)
       ISHELL=ISHELLV(ILOC)
       IF(MD.EQ.0)THEN
         LSPNB=KSPN
         LSPNE=KSPN
       ELSE IF(ABS(MD).EQ.1)THEN
         LSPNB=1
         LSPNE=NSPN
       ELSE IF(ABS(MD).EQ.2)THEN
         LSPNB=1
         LSPNE=NSPN
       END IF
C ZERO PSI:
       DO MO=1,MAX_OCC
         DO IL=1,MAXUNSYM
           PSI(IL,MO,ILOC)=0.0D0
         END DO
       END DO
       MLOCAL=0

C
C UNSYMMETRIZE THE WAVEFUNCTIONS....
C
!       ITOT=0
       IWF=0
       DO 1020 LSPN=LSPNB,LSPNE
        KSALC=0
        DO 1010 K_REP=1,N_REP
C
C CALCULATE ARRAY LOCATIONS:
C
         DO 5 K_ROW=1,NDMREP(K_REP)
          KSALC=KSALC+1
    5    CONTINUE
         INDEX=INDBEG(ISHELLA,K_REP)
         DO 20 LI =0,LSYMMAX(IFNCT)
          DO 15 IBASE=1,N_CON(LI+1,IFNCT)
           DO 10 IQ=1,N_SALC(KSALC,LI+1,ISHELL)
            INDEX=INDEX+1
            IND_SALC(IQ,IBASE,LI+1)=INDEX
   10      CONTINUE
   15     CONTINUE
   20    CONTINUE
         IF(MD.EQ.0)THEN
           MSTOT(1)=NS_TOT(K_REP)
           MSTOT(LSPN)=NFRM(LSPN)
         ELSE IF(MD.EQ.1)THEN
!<<<<< Force0.0
! If MD = 1, get coefs for occupied orbitals;
! start spin down orbital indices at N(occupied spin up)+1
!
           MSTOT(1)=NFRM(1) !NS_TOT(K_REP)
           MSTOT(LSPN)=NFRM(LSPN)
         ELSE IF(MD.EQ.2)THEN 
           MSTOT(1)   =NS_TOT(K_REP)
           MSTOT(NSPN)=NS_TOT(K_REP)
         END IF
         DO 1000 IOCC=1,MSTOT(LSPN)
!           ITOT=ITOT+1
           I_SALC=KSALC-NDMREP(K_REP)
           DO 950 IROW=1,NDMREP(K_REP)
             I_SALC=I_SALC+1
             IWF=IWF+1
             I_LOCAL=0
             DO 900 LI=0,LSYMMAX(IFNCT)
              DO 890 MU=1,NDEG(LI+1)
               IMS=MU+NDEG(LI+1)*(I_SITE-1)
               DO 880 IBASE=1,N_CON(LI+1,IFNCT)
                I_LOCAL=I_LOCAL+1
                MLOCAL=MAX(I_LOCAL,MLOCAL) !always equals to total orbitals!
                PSI(I_LOCAL,IWF,ILOC)=0.0D0
                IQ_BEG=IND_SALC(1,IBASE,LI+1)-1
                DO 800 IQ=1,N_SALC(KSALC,LI+1,ISHELL)
                 PSI(I_LOCAL,IWF,ILOC)=PSI(I_LOCAL,IWF,ILOC)+
     &           PSI_COEF(IQ+IQ_BEG,IOCC,K_REP,LSPN)*
     &           U_MAT(IMS,IQ,I_SALC,LI+1,ILOC)
  800           CONTINUE
  880          CONTINUE
  890         CONTINUE
  900        CONTINUE
             IF(I_LOCAL.GT.MAXUNSYM)THEN
              PRINT*,'UNRAVEL: MAXUNSYM MUST BE AT LEAST:',I_LOCAL
              CALL STOPIT
             END IF
  950      CONTINUE
C NOW ROTATE INTO THE FERMI ORBITALS
 1000    CONTINUE
 1010   CONTINUE
 1020  CONTINUE

!       ADD1=0.0D0
!       ADD2=0.0D0
!       ADD3=0.0D0
      !PRINT*,'BEFORE 2010',MD,KORB
      
      DO 2010 LSPN=LSPNB,LSPNE
       IF(MD.EQ.0) THEN
        IORBB=KORB
        IORBE=KORB
        if(kspn.eq.2) then
          if(mstot(1).eq.0)then
            print*,'mstot(1) is zero'
          end if
          iorbb=korb-mstot(1)
          iorbe=korb-mstot(1)
          IF(iorbb.le.0)then
            print*,'iorbb.le.0',irank
            call stopit
          end if
        end if
       ELSE
        IORBB=1
        IORBE=MORB(LSPN)
       END IF
       DO 2000 J_LOCAL=1,MLOCAL
        DO IORB=IORBB,IORBE
          PHILOC(IORB)=0.0D0
          IWF=0
          IF((LSPN.EQ.2).and.(MD.ne.0))IWF=MSTOT(1)
C         PRINT*,'IWF,MD:',IWF,MD
          DO JWF=1,MORB(LSPN)
            IWF=IWF+1
C           PRINT*,IORB,JWF,LSPN,IWF,TMAT(IORB,JWF,LSPN)
!            ADD1=ADD1+ABS(TMAT(IORB,JWF,LSPN))
C            IF(MD.EQ.0)THEN
C             PRINT*,'MD EQ 0', IORB,JWF,LSPN,TMAT(IORB,JWF,LSPN)
C            END IF
!            ADD2=ADD2+ABS(PSI(J_LOCAL,IWF,ILOC))
!            ADD3=ADD3+1.
            PHILOC(IORB)=PHILOC(IORB)
     &              +PSI(J_LOCAL,IWF,ILOC)*TMAT(IORB,JWF,LSPN)
          END DO
        END DO
        JORB=0
        IF((LSPNB.NE.LSPNE).AND.(LSPN.EQ.2)) THEN
          JORB=MSTOT(1)
          IF(MD.EQ.1) JORB=MORB(1) !<<<< Added from Force0.0 merge
        END IF
        DO IORB=IORBB,IORBE
          JORB=JORB+1
          PSI(J_LOCAL,JORB,ILOC)=PHILOC(IORB)  !COMMENT OUT ROTATION HERE
        END DO
 2000  CONTINUE
 2010 CONTINUE
      IF(.FALSE.) THEN !cmd testing timing of ADD statements
!      IF(ADD1*ADD2*ADD3.LE.0.001)THEN
C     PRINT*,IRANK,'UNRVL:',MD,LSPN,ADD1,ADD2,ADD3
        ADDTMAT=0.0D0     
        TMATMAX=0.0D0
        DO JWF =1,MORB(KSPN)
          DO IORB=1,MORB(KSPN)
            ADDTMAT=ADDTMAT+ABS(TMAT(IORB,JWF,KSPN))
            TMATMAX=MAX(TMAT(IORB,JWF,2),TMATMAX)
          END DO
        END DO
C       PRINT*,'ADDTMAT:',ADDTMAT,ADDUNRAVEL1
C       PRINT*,'TMATMAX:',TMATMAX
        PRINT*,'WHY SHOULD I BE STOPPING HERE?',IRANK,LSPN,KSPN,MD
        CALL STOPIT
      END IF
      IF(MD.EQ.0)THEN
        DO JWF=2,NWF
          DO J_LOCAL=1,MLOCAL
            PSI(J_LOCAL,JWF,ILOC)=0.0D0
          END DO
        END DO
C       DO J_LOCAL=1,MLOCAL
C         PSI(J_LOCAL,KORB,ILOC)=PSI(J_LOCAL,1,ILOC)
C         PSI(J_LOCAL, 1  ,ILOC)=0.0D0
C       END DO
      END IF
       RETURN
       END
