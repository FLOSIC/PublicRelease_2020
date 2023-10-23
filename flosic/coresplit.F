C UTEP Electronic Structure Lab (2020)
C
C
C
       SUBROUTINE CORESPLIT(MODE,ECORE,EMIN,EMAX,CHGVAL)
C
C  GETS THE RIGHT NUMBER OF VALENCE STATES AND ENERGY WINDOW
C 
       use debug1
       use common2,only : E_UP, E_DN, NSPN
       use common5,only : OCCUPANCY, N_OCC, NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:47 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MODE, I, IOCC, IREP, IROW, IS, ISORT, IWF, JSORT,
     & JWF, NVIRT
       REAL*8 :: SYMBOL , ECORE, EMIN, EMAX, CHGVAL, CCORE, CHARGE,
     & CTOT, EFMAX, EFMIN, EGAP, EVALSAV, GAP
       SAVE
       PARAMETER (EGAP=1.0D0)
       DIMENSION EVALSAV(MAX_OCC),GAP(MAX_OCC)
C
       IF(MODE.EQ.1) THEN
       EFMIN=MIN(EFERMI(1),EFERMI(NSPN))
       EFMAX=MAX(EFERMI(1),EFERMI(NSPN))
       ECORE=-9999.D0
       IF(DEBUG) THEN
         write(6,*)'CS: ',E_UP,E_DN,EFMIN,EFMAX,ECORE
C         write(6,*)'CS: ',N_REP,NWF, NWFS(1),NWFS(2)
         write(6,*)'CS: ',N_REP,NWF
       ENDIF
        DO I=1,NWF
         EVALSAV(I)=EVLOCC(I)
        ENDDO
C SORT EIGENVALUES
C
        NVIRT=NWF+10
        DO ISORT=1,NVIRT
         DO JSORT=ISORT+1,NVIRT
           IF (EVALSAV(JSORT).LT.EVALSAV(ISORT)) THEN
            CALL SWAP(EVALSAV(ISORT),EVALSAV(JSORT))
          END IF
         ENDDO
        ENDDO
        DO I=1,NVIRT-1
         GAP(I)=EVALSAV(I+1)-EVALSAV(I)
         IF(DEBUG) write(6,*)'GAP : ',i,gap(i),EVALSAV(I)
        ENDDO
        DO I=NVIRT-1,1,-1
         IF(GAP(I).GT.EGAP) THEN
            ECORE=EVALSAV(I)+0.5D0*GAP(I)
            GOTO 200
         ENDIF
        ENDDO
  200   CONTINUE
        IF(ECORE.LT.EVALSAV(1)) THEN
         write(6,*)'ERROR IN CORESPLIT: NO GAP FOUND!'
         EMIN=EVALSAV(1)
         EMAX=EFMAX+EMAX
         CHGVAL=E_UP+E_DN
        ELSE
         EMIN=ECORE
         EMAX=EFMAX+EMAX
        ENDIF
        RETURN
       ENDIF    !MODE=1 
C
       PRINT'(A,F9.3)',
     &    'STATES BELOW ARE EXCLUDED FROM SPIN-ORBIT: ',ECORE
       CHARGE=0.0D0
       CCORE=0.0D0
       IWF=0
       JWF=0
       DO IS=1,NSPN
        DO IREP=1,N_REP
          DO IOCC=1,N_OCC(IREP,IS)
             DO IROW=1,NDMREP(IREP)
             JWF=JWF+1
             END DO
          IWF=IWF+1
          CHARGE=CHARGE+OCCUPANCY(IWF)*NDMREP(IREP)
C       write(6,*)'CHG: ',IWF,CHARGE,OCCUPANCY(IWF),NDMREP(IREP),
C     &                EVLOCC(JWF) 
          IF(EVLOCC(JWF).LT.ECORE) THEN
            CCORE=CCORE+OCCUPANCY(IWF)*NDMREP(IREP)
          ENDIF
          ENDDO
        ENDDO
       END DO                                          
       CTOT=E_UP+E_DN
       CHGVAL=CTOT-CCORE
       write(6,*)'NUMBER OF CORE/VAL ELECTRONS: ', CCORE,CHGVAL
       RETURN
       END 
