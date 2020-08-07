C UTEP Electronic Structure Lab (2020)

********************************************************************
C
C      THIS ROUTINE FINDS SPIN, REPRESENTATION AND THE INDEX OF THE BASIS
C      HOMO OF THE MOLECULE AND THE HOMO+INDEX STATES.
C                                          TB 04/03
C
      SUBROUTINE FINDSTATE(INDEX,N,IS,IR,IB,IERR)
       use common2,only : NSPN
       use common3,only : RMAT
       use common8,only : REP, N_REP, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:45 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: INDEX, N, IS, IR, IB, IERR, I, IBS, IREPR, IRP, ISP,
     & ISPIN, IV, IVIRT, IWF, IWFHOMO, J, JVIRT, KWF, NBASE, NDIM
       REAL*8 :: EV , EVL
      SAVE 
      LOGICAL EXIST
      CHARACTER*20 TRASH
      DIMENSION INDEX(N),IS(N),IR(N),IB(N),EV(MXSPN*MAX_OCC)
      DIMENSION ISP(MXSPN*MAX_OCC),IRP(MXSPN*MAX_OCC),
     &           IBS(MXSPN*MAX_OCC),EVL(MXSPN*MAX_OCC)
 
        IERR=0
        INQUIRE(FILE='EVALUES',EXIST=EXIST) 
        IF(EXIST) THEN
          OPEN(90,FILE='EVALUES',STATUS='OLD')
          REWIND(90)
        ELSE
          WRITE(6,*)'EVALUES DOES NOT EXIST, CALLING'
          WRITE(6,*)'OLD ROUTINE TO FIND HOMO AND OCCUPIED STATES'
          WRITE(6,*)'NEED EVALUES TO FIND LUMO'
          IERR=1
          RETURN
        END IF
        IWF=0

        DO 300 ISPIN=1,NSPN
          READ(90,240) TRASH
          DO 200 IREPR=1,N_REP
            IF(NS_TOT(IREPR).LE.0) GO TO 200
            READ(90,240) TRASH
            READ(90,*) NDIM, NBASE
             JVIRT=0
             READ(90,*) (EV(IVIRT),IVIRT=1,NBASE)
              DO IV=1,NBASE
                IWF=IWF+1
                JVIRT=JVIRT+1
                EVL(IWF)=EV(IV)
                ISP(IWF)=ISPIN
                IRP(IWF)=IREPR
                IBS(IWF)=JVIRT
              END DO
 200      CONTINUE
 300    CONTINUE
 240    FORMAT(A20)
        DO I=1,IWF
          DO J=1,IWF
           IF(EVL(J)-EVL(I).GT.1.0D-6) THEN  
             CALL SWAP(EVL(I),EVL(J))
             CALL ISWAP(ISP(I),ISP(J))
             CALL ISWAP(IRP(I),IRP(J))
             CALL ISWAP(IBS(I),IBS(J))
           END IF
          END DO
        END DO
        CLOSE(90)
        CALL FINDHOMO(IWF,EVL,IWFHOMO) 
        DO I=1,N
            KWF=IWFHOMO+INDEX(I)
            IS(I)=ISP(KWF)
            IR(I)=IRP(KWF)
            IB(I)=IBS(KWF)
        END DO
      RETURN
      END
