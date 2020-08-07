C UTEP Electronic Structure Lab (2020)
     
********************************************************************

      SUBROUTINE FINDHOMO(IWF,EVL,IWFHOMO)
       use common2,only : NSPN
       use common5,only : EFERMI
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:44 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IWF, IWFHOMO, I
       REAL*8 :: EVL , DIFF, DIFFNEW, EF
      SAVE
      DIMENSION EVL(IWF) 
      EF=EFERMI(1)
      IF(NSPN.EQ.2) EF=MAX(EF,EFERMI(NSPN))
      DIFF=1.0D30
       DO I=1,IWF
         DIFFNEW=EF-EVL(I)
         IF(DIFFNEW.GE.-1.0D-4) THEN
           IF(DIFFNEW.LT.DIFF) THEN
            DIFF=DIFFNEW
            IWFHOMO=I
           END IF 
         END IF 
       END DO 
       WRITE(6,*) 'IWFHOMO = ', IWFHOMO
      RETURN
      END
