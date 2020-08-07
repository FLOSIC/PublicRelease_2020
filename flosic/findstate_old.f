C UTEP Electronic Structure Lab (2020)


********************************************************************
C
C      THIS ROUTINE FINDS SPIN, REPRESENTATION AND THE INDEX OF THE BASIS
C      HOMO OF THE MOLECULE AND THE HOMO+INDEX STATES.
C                                          TB 04/03
C
      SUBROUTINE FINDSTATE_old(INDEX,N,IS,IR,IB)
       use common2,only : NSPN
       use common5,only : N_OCC, EVLOCC
       use common8,only : REP, N_REP, NDMREP
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:45 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: INDEX, N, IS, IR, IB, I, IBS, IDEG, IREPR, IRP, ISP,
     & ISPIN, IVIRT, IWF, IWFHOMO, J, JVIRT, KWF
       REAL*8 :: EVL 
      SAVE 
      DIMENSION INDEX(N),IS(N),IR(N),IB(N)
      DIMENSION ISP(MAX_OCC),IRP(MAX_OCC),IBS(MAX_OCC),EVL(MAX_OCC)
      
        WRITE(6,*)'FINDSTATE:' ,NSPN,N_REP
        IWF=0
        DO ISPIN=1,NSPN
          DO IREPR=1,N_REP
            JVIRT=0
            WRITE(6,*)'N_OCC', IREPR,ISPIN,N_OCC(IREPR,ISPIN)
            DO IVIRT=1,N_OCC(IREPR,ISPIN)
               DO IDEG=1,NDMREP(IREPR)
                IWF=IWF+1
                JVIRT=JVIRT+1
                EVL(IWF)=EVLOCC(IWF)
                ISP(IWF)=ISPIN
                IRP(IWF)=IREPR
                IBS(IWF)=JVIRT
               END DO
            END DO
          END DO
        END DO
        WRITE(6,*) IWF
        DO I=1,IWF
          DO J=1,IWF
           IF(EVL(J).GT.EVL(I)) THEN  
             CALL SWAP(EVL(I),EVL(J))
             CALL ISWAP(ISP(I),ISP(J))
             CALL ISWAP(IRP(I),IRP(J))
             CALL ISWAP(IBS(I),IBS(J))
           END IF
          END DO
        END DO
        CALL FINDHOMO(IWF,EVL,IWFHOMO) 
        DO I=1,N
           KWF=IWFHOMO+INDEX(I)
           IS(I)=ISP(KWF)
           IR(I)=IRP(KWF)
           IB(I)=IBS(KWF)
        END DO
      RETURN
      END
