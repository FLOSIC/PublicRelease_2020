C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE WFOUT2
C ORIGINAL VERSION BY MARK R PEDERSON
C 06/23/97 Converted to subroutine. David C. Patton
C
C WRITE WAVEFC
C
       use common2,only : WFFILE, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF,
     &   NWFS, EVLOCC
       use common8,only : REP, N_REP, NS_TOT
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:02:00 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I_OCC, I_REP, ITOT, IWF
      SAVE
C
      WFFILE='WFOUT2'
      OPEN(98,FILE=WFFILE,FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND(98)
      WRITE(98,*) NSPN
      WRITE(98,*) NWF,(NWFS(ISPN),ISPN=1,NSPN)
      WRITE(98,*) (EVLOCC(IWF), IWF=1,NWF)
      ITOT=0
      WRITE(98,*) N_REP
      DO 870 ISPN=1,NSPN
        DO 860 I_REP=1,N_REP
c          ITOT=ITOT+1
          WRITE(98,*) N_OCC(I_REP,ISPN),NS_TOT(I_REP)
          WRITE(98,*) (OCCUPANCY(I_OCC+ITOT),
     &             I_OCC=1,N_OCC(I_REP,ISPN))
          ITOT=ITOT+N_OCC(I_REP,ISPN)
          DO 850 IWF=1,N_OCC(I_REP,ISPN)
            WRITE(98,*) (PSI_COEF(I,IWF,I_REP,ISPN),
     &              I=1,NS_TOT(I_REP))
850       CONTINUE
860     CONTINUE
870   CONTINUE
      CLOSE(98)
      RETURN
      END
