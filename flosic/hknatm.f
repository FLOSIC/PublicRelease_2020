C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE HKNATM(IFNCT,HKNTAB)
C
C CREATE TABLE OF ATOMIC KINETIC ENERGY INTEGRALS
C
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:49 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFNCT, IBARE, ICON, JBARE, JCON, L, L1
       REAL*8 :: HKNTAB , ABET, ASRC, ASUM, FAC1, FAC2, P, PI
        SAVE
        DIMENSION HKNTAB(MAX_CON,MAX_CON,LDIM)
C 
        PI= 4*ATAN(1.0D0) 
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          DO JCON= ICON,N_CON(L1,IFNCT)
           HKNTAB(JCON,ICON,L1)= 0.0D0
          END DO
         END DO
        END DO
C
        DO IBARE= 1,N_BARE(IFNCT)
         DO JBARE= 1,N_BARE(IFNCT)
          ABET= BFALP(JBARE,IFNCT)
          ASUM= BFALP(IBARE,IFNCT)+ABET
          ASRC= 1.0D0/ASUM
          FAC1= 2*PI*SQRT(PI*ASRC)
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           FAC1= 0.5D0*ASRC*(2*L+1)*FAC1
           FAC2= 0.5D0*ASRC*(2*L+3)*FAC1
           DO ICON= 1,N_CON(L1,IFNCT)
            DO JCON= 1,N_CON(L1,IFNCT)
             P= BFCON(IBARE,ICON,L1,IFNCT)*BFCON(JBARE,JCON,L1,IFNCT)
             HKNTAB(JCON,ICON,L1)= HKNTAB(JCON,ICON,L1)
     &                            +ABET*P*((2*L+3)*FAC1-2*ABET*FAC2)
            END DO
           END DO
          END DO
         END DO
        END DO
       RETURN
       END
