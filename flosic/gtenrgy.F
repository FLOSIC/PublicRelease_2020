C UTEP Electronic Structure Lab (2020)
C
C *****************************************************************
C
       SUBROUTINE GTENRGY
C
C WRITTEN BY MARK R PEDERSON (1986-1989), MODIFIED BY DIRK POREZAG
C
C       use debug1
       use blas_module,only : TRIPLE_PROD2
       use global_inputs,only : SOLVENT1,PCM1,DMAT1,EFP1
       use for_diag1
       use hstor1,only : hstor
       use common2,only : RIDT, IFUIDT, NIDENT, ZELC, ISPN, NSPN,
     & ETOTAL, ENNUC, ELOCAL, ECOUL, EKINONL, ERGXL, EKINONL2,
     & ERGXN, ERGCL, ERGCN, ERGFLD, ESOLC, ENUCSOL, ESOLTOT,EPCM
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI
       use common8,only : REP, N_REP, NDMREP, NS_TOT
       use locorb,only  : ZSIC
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I_BASIS, I_NUC, I_REP, IERR, IERROR, IFNCT, IOUT,
     & IREC, IWF, J, J_NUC, JFNCT, K, K_NUC, K_VIRT, M_NUC, MM, MREC,
     & MSITES, MSPN, N_BASIS, NBAS, NREC
       REAL*8 :: ADD , ENERGY, EREP, ERGBEST, ESIC, ESICTOTAL, EXCCOR,
     & OCC_K, R_NUCA, RABS, SHIFT, TOLER, VEC
       SAVE
       DIMENSION VEC(3),R_NUCA(3,MX_GRP)
       DIMENSION MSITES(1)
       REAL*8,ALLOCATABLE :: COEF(:,:)
       COMMON/SICENRGY/ESIC(11,MAX_CON,LDIM,MAX_IDENT),
     &                   SHIFT(11,Max_CON,LDIM,MAX_IDENT),ESICTOTAL
       DATA TOLER/1.0D-5/
C
C CALL OVERLAP TO GET KINETIC+NONLOCAL ENERGY 
C
       IERROR=0
       CALL OVERLAP(3)

      !This block is moved down
      !IF(DMAT1)THEN
      !  CALL GTKINRGY_DMAT
      !  EKINONL=EKINONL2
      !ELSE
C
C CALCULATE NUCLEAR REPULSION ENERGY
C
       ENERGY=0.0D0
       DO 90 I_NUC=1,NIDENT
        IFNCT=IFUIDT(I_NUC)
        CALL GASITES(1,RIDT(1,I_NUC),M_NUC,R_NUCA,MSITES)
        DO 60 J_NUC=1,M_NUC
         DO 50 K_NUC=J_NUC+1,M_NUC
          VEC(1)=R_NUCA(1,K_NUC)-R_NUCA(1,J_NUC)
          VEC(2)=R_NUCA(2,K_NUC)-R_NUCA(2,J_NUC)
          VEC(3)=R_NUCA(3,K_NUC)-R_NUCA(3,J_NUC)
          RABS=SQRT((VEC(1)**2+VEC(2)**2+VEC(3)**2))
          IF (RABS.LT.TOLER) THEN
           IERROR=1
           write(6,*)'PROBABLE ERROR(1) IN GTENRGY'
           write(6,*)'RABS=',RABS,' M_NUC=',M_NUC
           write(6,*)'I_NUC,R:',I_NUC,(R_NUCA(J,I_NUC),J=1,3)
           write(6,*)'J_NUC,R:',J_NUC,(R_NUCA(J,J_NUC),J=1,3)
           write(6,*)'K_NUC,R:',K_NUC,(R_NUCA(J,K_NUC),J=1,3)
          END IF
          ENERGY=ENERGY+ZELC(IFNCT)*ZELC(IFNCT)/RABS
   50    CONTINUE
   60   CONTINUE
        MM=M_NUC
        DO 80 J_NUC=I_NUC+1,NIDENT
         JFNCT=IFUIDT(J_NUC)
         CALL GASITES(1,RIDT(1,J_NUC),M_NUC,R_NUCA,MSITES)
         DO 70 K_NUC=1,M_NUC
          VEC(1)=RIDT(1,I_NUC)-R_NUCA(1,K_NUC)
          VEC(2)=RIDT(2,I_NUC)-R_NUCA(2,K_NUC)
          VEC(3)=RIDT(3,I_NUC)-R_NUCA(3,K_NUC)
          RABS=SQRT((VEC(1)**2+VEC(2)**2+VEC(3)**2))
          IF (RABS.LT.TOLER) THEN
           IERROR=1
           write(6,*)'PROBABLE ERROR(2) IN GTENRGY'
           write(6,*)'RABS=',RABS,' M_NUC,MM=',M_NUC,MM
           write(6,*)'I_NUC,R:',I_NUC,(RIDT(J,I_NUC),J=1,3)
           write(6,*)'J_NUC,R:',K_NUC,(RIDT(J,J_NUC),J=1,3)
           write(6,*)'K_NUC,R:',K_NUC,(R_NUCA(J,K_NUC),J=1,3)
          END IF
          ENERGY=ENERGY+ZELC(IFNCT)*ZELC(JFNCT)*MM/RABS
   70    CONTINUE
   80   CONTINUE
   90  CONTINUE
       ENNUC=ENERGY 

       IF(DMAT1)THEN
         CALL GTKINRGY_DMAT
         EKINONL=EKINONL2
       ELSE       
C
C CALCULATE KINETIC+NONLOCAL ENERGY 
C
  100  EKINONL=0.0D0
       K_VIRT=0
       DO 200 ISPN=1,NSPN
        K=0
        DO 190 I_REP=1,N_REP
         NBAS=NS_TOT(I_REP)
         ALLOCATE(COEF(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'gtenrgy:Error allocating COEF'
         ENDIF
!         ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
!         IF(IERR.NE.0)THEN
!           WRITE(6,*)'gtenrgy:Error allocating Ham'
!         ENDIF
         ALLOCATE(AOVER(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'gtenrgy:Error allocating Overlap'
         ENDIF
         DO 120 I=1,NS_TOT(I_REP)
          DO 110 J=I,NS_TOT(I_REP)
           K=K+1
!           AHAM(J,I)=HSTOR(K,1)
           AOVER(J,I)=HSTOR(K,1)
           AOVER(I,J)=AOVER(J,I)
  110     CONTINUE
  120    CONTINUE
         N_BASIS=NS_TOT(I_REP)
         DO 140 IWF=1,N_OCC(I_REP,ISPN)
          K_VIRT=K_VIRT+1
          OCC_K=sqrt(OCCUPANCY(K_VIRT))
          DO 130 I_BASIS=1,N_BASIS
!           AOVER(IWF,I_BASIS)=PSI_COEF(I_BASIS,IWF,I_REP,ISPN)
!     &                      *SQRT(OCCUPANCY(K_VIRT))
           COEF(I_BASIS,IWF)=PSI_COEF(I_BASIS,IWF,I_REP,ISPN)
     &                      *OCC_K
!     &                      *SQRT(OCCUPANCY(K_VIRT))
  130     CONTINUE
  140    CONTINUE
         EREP=0.0D0
!         DO 170 I_BASIS=1,N_BASIS
!          DO 160 J_BASIS=I_BASIS,N_BASIS
!           ADD=0.0D0
!           DO 150 IWF=1,N_OCC(I_REP,ISPN)
!            ADD=ADD+AOVER(IWF,I_BASIS)*AOVER(IWF,J_BASIS)
!  150      CONTINUE
!           IF (J_BASIS.NE.I_BASIS) ADD=ADD*2
!           EREP=EREP+ADD*AHAM(J_BASIS,I_BASIS)
!  160     CONTINUE
!  170    CONTINUE

         CALL TRIPLE_PROD2(N_BASIS,N_OCC(I_REP,ISPN),COEF,AOVER,ADD)
         EREP=EREP+ADD

         EKINONL=EKINONL+EREP*NDMREP(I_REP)
!         DEALLOCATE(AHAM,STAT=IERR)
!         IF(IERR.NE.0)THEN
!           WRITE(6,*)'gtenrgy:Error deallocating Ham'
!         ENDIF
         DEALLOCATE(AOVER,STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'gtenrgy:Error deallocating Overlap'
         ENDIF
         DEALLOCATE(COEF,STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'gtenrgy:Error deallocating Overlap'
         ENDIF
  190   CONTINUE
  200  CONTINUE
       ENDIF
C
c       CALL SICSET(2,1)
       EKINONL=(2*EKINONL)/NSPN
       EXCCOR=ERGXL+ERGXN+ERGCL+ERGCN
       ETOTAL=ENNUC+ELOCAL+ECOUL+EKINONL+EXCCOR+ERGFLD
       ETOTAL=ETOTAL+ENUCSOL+ESOLC
       ESICTOTAL=ZSIC
       IF(PCM1) ETOTAL=ETOTAL+EPCM
       PRINT 1000,'================================================'
       PRINT 1000,'TOTAL ENERGY:               ',ETOTAL
       IF(SOLVENT1) THEN 
         PRINT 1000,'TOTAL ENERGY WITH SOLVENTS: ',ESOLTOT +ETOTAL
       ENDIF
       PRINT 1000,'CONTRIBUTIONS:'
       IF(SOLVENT1) THEN
         PRINT 1000,'SOLVENT-ELECTRON COULOMB:   ',ESOLC
         PRINT 1000,'SOLVENT-NUCLEAR COULOMB:    ',ENUCSOL
         PRINT 1000,'SOLVENT-SOLVENT COULOMB:    ',ESOLTOT
       ENDIF
       PRINT 1000,'NUCLEAR REPULSION:          ',ENNUC
       PRINT 1000,'LOCAL POTENTIAL:            ',ELOCAL 
       PRINT 1000,'MEAN-FIELD COULOMB:         ',ECOUL 
       PRINT 1000,'KINETIC+NONLOCAL POTENTIAL: ',EKINONL
       PRINT 1000,'EXCHANGE-CORRELATION:       ',EXCCOR
       PRINT 1000,'EXTERNAL ELECTRIC FIELD:    ',ERGFLD 
       IF(PCM1) THEN
         PRINT 1000,'PCM POLARIZATION ENERGY:    ',EPCM
       ENDIF
       PRINT 1000,'SELF-INTERACTION ENERGY:    ',ESICTOTAL
       PRINT 1000,'TOTAL ENERGY + SIC ENERGY:  ',ETOTAL+ESICTOTAL
       PRINT 1000,'================================================'
        call flush(iout)
1000   FORMAT(A,F20.6)
       IF (IERROR.EQ.1) CALL STOPIT
C
C COPY HAMOLD TO HAMBST. CHECK IF HAMOLD IS VALID.
C CHECK ALSO IF HAMOLD AND HAMBST ARE COMPATIBLE
C
       OPEN(81,FILE='HAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND(81)
       READ(81,END=280)NREC,MSPN
       IF (MSPN.NE.NSPN) GOTO 280
C
       ERGBEST=ETOTAL+1.0D30
       MREC=NREC
       MSPN=NSPN
       OPEN(82,FILE='HAMBST',FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND(82)
       READ(82,END=230)MREC,MSPN,ERGBEST
  230  IF ((MREC.NE.NREC).OR.(MSPN.NE.NSPN)
     &                   .OR.(ETOTAL.LT.ERGBEST)) THEN
        REWIND(82)
        WRITE(82)NREC,NSPN,ETOTAL
        DO 240 ISPN=1,NSPN
         READ(81,END=260)(HSTOR(IREC,1),IREC=1,NREC)
         WRITE(82)(HSTOR(IREC,1),IREC=1,NREC)
  240   CONTINUE
       END IF
       CLOSE(81)
       CLOSE(82)
       GOTO 300
C
C ERROR HANDLING
C
  260  CLOSE(82,STATUS='DELETE')
  280  CLOSE(81,STATUS='DELETE')
       write(6,*)'GTENRGY: HAMOLD IS INCOMPATIBLE OR NOT AVAILABLE'
  300  CONTINUE
       RETURN
       END
