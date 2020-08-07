C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE OVLONS(IFNCT,OVLTAB)
C
C CREATE TABLE OF ONSITE OVERLAP INTEGRALS
C
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFNCT, IBARE, ICON, JBARE, L, L1
       REAL*8 :: OVLTAB , ALP, ALRC, FACO, P, PI
        SAVE
        DIMENSION OVLTAB(MAX_CON,LDIM)
C 
        PI= 4*ATAN(1.0D0) 
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          OVLTAB(ICON,L1)= 0.0D0
         END DO
        END DO
        DO IBARE= 1,N_BARE(IFNCT)
         DO JBARE= IBARE,N_BARE(IFNCT)
          ALP= BFALP(IBARE,IFNCT)+BFALP(JBARE,IFNCT)
          ALRC= 1.0D0/ALP
          FACO= 2*PI*SQRT(PI*ALRC)
          IF (IBARE .NE. JBARE) THEN
           FACO= 2*FACO
          END IF
          DO L= 0,LSYMMAX(IFNCT)
           FACO= 0.5D0*ALRC*(2*L+1)*FACO
           L1= L+1
           DO ICON= 1,N_CON(L1,IFNCT)
            P= BFCON(IBARE,ICON,L1,IFNCT)*BFCON(JBARE,ICON,L1,IFNCT)
            OVLTAB(ICON,L1)= OVLTAB(ICON,L1)+FACO*P
           END DO
          END DO
         END DO
        END DO
        RETURN
       END
