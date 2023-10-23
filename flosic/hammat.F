C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
      SUBROUTINE HAMMAT(LPTS,MPTS,MDCL,POTDV,INDX)
C
C 02/12/97 David Clay Patton
C      
      use hstor1,only : HSTOR
      use xtmp2a,only : PSIBR,PSIB,VOL,QR
      use mesh1,only : RMSH
      use common3,only : NGRP
      use common5,only : PSI
      use common8,only : REP, N_REP, NDMREP
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:52 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LPTS, MPTS, MDCL, INDX, IBG, IGP, IPTS, ISS, JGP,
     & JNDX, JREP, JSS, KREP, KRP, MGRP, NBAS, NRP
       REAL*8 :: POTDV , DOT
      SAVE
      LOGICAL LSETUP
      DIMENSION POTDV(*)
C      COMMON/TMP2/PSIBR(MPBLOCK,5*NDH),PSIB(MPBLOCK)
C     & ,VOL(MPBLOCK),QR(3,MPBLOCK)
C      REAL*8,ALLOCATABLE :: PSIBR(:,:),PSIB(:),VOL(:),QR(:,:)
C
      MGRP=1
      IF (MDCL.EQ.1) MGRP=NGRP
C
      LSETUP=.FALSE. 
      DO IPTS=1,MPTS
       VOL(IPTS)=POTDV(LPTS+IPTS)/MGRP 
       QR(1,IPTS)=RMSH(1,LPTS+IPTS)
       QR(2,IPTS)=RMSH(2,LPTS+IPTS)
       QR(3,IPTS)=RMSH(3,LPTS+IPTS)
      END DO
C
C UPDATE MATRIX ELEMENTS HERE:
C
      DO 100 IGP=1,MGRP
       JGP=IGP
       INDX=0
       DO 200 KREP=1,N_REP
        JREP=KREP
        IF (MDCL.EQ.1) THEN
         NRP=1
        ELSE
         NRP=NDMREP(KREP)
        END IF
        CALL GETBAS(LSETUP,MPTS,JGP,QR,PSIBR,JREP,NBAS)
        JNDX=INDX
        DO 300 KRP=1,NRP
         INDX=JNDX
         IBG=(KRP-1)*NBAS
         DO 400 ISS=1,NBAS
          DO IPTS=1,MPTS
           PSIB(IPTS)=PSIBR(IPTS,ISS+IBG)*VOL(IPTS)/NRP
          END DO
C
C NOTE: WE COULD MAKE USE OF SPARCITY HERE. IF PSIB=0 AT EACH
C       MESHPOINT, THE NEXT LOOP CAN BE SKIPPED - INDX MUST BE
C       RESET OF COURSE
C
          DO JSS=ISS,NBAS
           DOT=0.0D0
           DO IPTS=1,MPTS
            DOT=DOT+PSIB(IPTS)*PSIBR(IPTS,JSS+IBG)
           END DO
           INDX=INDX+1
           HSTOR(INDX,2)=HSTOR(INDX,2)+DOT
          END DO
400      CONTINUE
300     CONTINUE
200    CONTINUE
100   CONTINUE 
      RETURN
      END
