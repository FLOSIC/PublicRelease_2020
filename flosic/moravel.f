C UTEP Electronic Structure Lab (2020)
       SUBROUTINE UNRAVEL3
     &(MD,KSPN,KORB,IFNCT,ISHELLA,I_SITE,RVEC,RVECI,N_NUC,ILOC)
C ORIGINALLY WRITTEN BY MARK R PEDERSON (1985)
       use debug1
       use common2,only : N_CON, LSYMMAX, ISPN, NSPN
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF
       use common7,only : T1UNRV, T2UNRV
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC, INDBEG
     &                   ,NS_TOT
!SIC module
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MD, KSPN, KORB, IFNCT, ISHELLA, I_SITE, N_NUC, ILOC,
     & I_LOCAL, I_SALC, IBASE, ICALL1, ICALL2, IL, IMS, IND_SALC,
     & INDEX, IOCC, IORB, IORBB, IORBE, IQ, IQ_BEG, IROW, ISHELL,
     & ISHELLV, ISPNB, ISPNE, ITOT, IWF, J_LOCAL, JORB, JWF, K_REP,
     & K_ROW, KSALC, LI, MLOCAL, MO, MSTOT, MU, NDEG
       REAL*8 :: RVEC , RVECI, PHILOC, TIMER1, TIMER2
!       COMMON/LOCORB/TMAT(MAX_OCC,MAX_OCC,MXSPN),MORB(2),ZSIC,IRBSIC
!       COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),ZTZL,JJJJJJ
       save
C       INCLUDE 'commons.inc'
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
C       PRINT *,'IN UNRAVEL',IRANK

C       if(kspn.eq.2) then
C        print*,'KSPN!!!'
C       END IF



C Note: MD=0 puts result in PSI(:,1,ILOC)


       IF(I_SITE.EQ.1)THEN
        ICALL1=ICALL1+1
        CALL GTTIME(TIMER1)
        CALL OBINFO(1,RVEC,RVECI,N_NUC,ISHELLV(ILOC))
        CALL GSMAT(ISHELLV(ILOC),ILOC)
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
              ISPNB=KSPN
              ISPNE=KSPN
       ELSE IF(ABS(MD).EQ.1)THEN
              ISPNB=1
              ISPNE=NSPN
       ELSE IF(ABS(MD).EQ.2)THEN
              ISPNB=1
              ISPNE=NSPN
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
       ITOT=0
       IWF=0
       DO 1020 ISPN=ISPNB,ISPNE
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
C
C END CALCULATION OF SALC INDICES FOR REPRESENTATION K_REP
C               IF(MD.EQ.0)MSTOT(ISPN)=N_OCC(K_REP,ISPN) !N_OCC is also the number of all orbitals!!!!!!
C               IF(MD.EQ.1)MSTOT(ISPN)=N_OCC(K_REP,ISPN)

               IF(MD.EQ.0)MSTOT(ISPN)=NFRM(ISPN) !Only works if no symmetry
               IF(MD.EQ.1)MSTOT(ISPN)=NFRM(ISPN)



               IF(MD.EQ.2)MSTOT(ISPN)=NS_TOT(K_REP)
         DO 1000 IOCC=1,MSTOT(ISPN)
          ITOT=ITOT+1
          I_SALC=KSALC-NDMREP(K_REP)
          DO 950 IROW=1,NDMREP(K_REP)
           I_SALC=I_SALC+1
           IWF=IWF+1
C           PRINT *,'IRANK',IRANK,IWF
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
     &         PSI_COEF(IQ+IQ_BEG,IOCC,K_REP,ISPN)*
     &         U_MAT(IMS,IQ,I_SALC,LI+1,ILOC)
  800         CONTINUE
  880        CONTINUE
  890       CONTINUE
  900      CONTINUE
           IF(I_LOCAL.GT.MAXUNSYM)THEN
            PRINT*,'UNRAVEL: MAXUNSYM MUST BE AT LEAST:',I_LOCAL
            CALL STOPIT
           END IF
  950     CONTINUE
C NOW ROTATE INTO THE FERMI ORBITALS
 1000    CONTINUE
 1010   CONTINUE
 1020  CONTINUE

      DO 2000 ISPN=ISPNB,ISPNE
          IF(MD.EQ.0) THEN
            IORBB=KORB
            IORBE=KORB
            if(kspn.eq.2) then
             iorbb=korb-mstot(1)
             iorbe=korb-mstot(1)
            end if
          ELSE
            IORBB=1
            IORBE=MORB(ISPN)
          END IF
      DO 2000 J_LOCAL=1,MLOCAL
      DO IORB=IORBB,IORBE
        PHILOC(IORB)=0.0D0
             IWF=0
             IF((ISPN.EQ.2).and.(MD.ne.0))IWF=MSTOT(1)
        DO JWF=1,MORB(ISPN)
             IWF=IWF+1
        PHILOC(IORB)=PHILOC(IORB)
     &              +PSI(J_LOCAL,IWF,ILOC)*TMAT(IORB,JWF,ISPN)
        END DO
      END DO
          JORB=0
          IF((ISPNB.NE.ISPNE).AND.(ISPN.EQ.2))JORB=MSTOT(1)
      DO IORB=IORBB,IORBE
          JORB=JORB+1
      PSI(J_LOCAL,JORB,ILOC)=PHILOC(IORB)  !COMMENT OUT ROTATION HERE
      END DO
 2000 CONTINUE
      IF(MD.EQ.0)THEN
             DO JWF=2,NWF
             DO J_LOCAL=1,MLOCAL
             PSI(J_LOCAL,JWF,ILOC)=0.0D0
             END DO
             END DO
C            DO J_LOCAL=1,MLOCAL
C            PSI(J_LOCAL,KORB,ILOC)=PSI(J_LOCAL,1,ILOC)
C            PSI(J_LOCAL, 1  ,ILOC)=0.0D0
C            END DO
      END IF


*      if((md.eq.0).and.(kspn.eq.2)) then
*       print*,'md=',md,'kspn=',kspn,'korb=',korb
*       do iwf=1,NWF
*        print'(100F10.5)',(psi(ibas,iwf,iloc),ibas=1,MLOCAL)
*       end do
*       stop
*      end if


       RETURN
       END
