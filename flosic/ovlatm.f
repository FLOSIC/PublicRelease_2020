C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE OVLATM(IFNCT,OVLTAB)
C
C CREATE TABLE OF ATOMIC OVERLAP INTEGRALS
C
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFNCT, IBARE, ICON, JBARE, JCON, L, L1
       REAL*8 :: OVLTAB , ALP, ALRC, FACO, P, PI
        SAVE
        DIMENSION OVLTAB(MAX_CON,MAX_CON,LDIM)
C 
        PI= 4*ATAN(1.0D0) 
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          DO JCON= 1,N_CON(L1,IFNCT)
           OVLTAB(JCON,ICON,L1)= 0.0D0
          END DO
         END DO
        END DO
C
        DO IBARE= 1,N_BARE(IFNCT)
         DO JBARE= 1,N_BARE(IFNCT)
          ALP= BFALP(IBARE,IFNCT)+BFALP(JBARE,IFNCT)
          ALRC= 1.0D0/ALP
          FACO= 2*PI*SQRT(PI*ALRC)
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           FACO= 0.5D0*ALRC*(2*L+1)*FACO
           DO ICON= 1,N_CON(L1,IFNCT)
            DO JCON= 1,N_CON(L1,IFNCT)
             P= BFCON(IBARE,ICON,L1,IFNCT)*BFCON(JBARE,JCON,L1,IFNCT)
             OVLTAB(JCON,ICON,L1)= OVLTAB(JCON,ICON,L1)+FACO*P
            END DO
           END DO
          END DO
         END DO
        END DO
        RETURN
       END
