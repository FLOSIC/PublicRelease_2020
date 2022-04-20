C UTEP Electronic Structure Lab (2020)
!> @file swapfgi.f
       SUBROUTINE SWAPFGI(LSPN,NBXX,IW)
       use for_diag1
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:03 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LSPN, NBXX, IW, IB, JB
       !REAL*8 :: EVSIC , GILO
       REAL*8,allocatable :: EVSIC(:),GILO(:,:)
!       COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
!     &  EVAL(NDH),SC1(NDH),SC2(NDH)
!       COMMON/SIC_DIAG/GILO(NDH,NDH,2),EVSIC(NDH,2)
C Copy values FILO -> GILO
       allocate(EVSIC(NBXX), GILO(NBXX,NBXX))
       IF(IW.EQ.0)THEN
         DO IB=1,NBXX
           !EVSIC(IB,LSPN)=EVAL(IB)
           EVSIC(IB)=EVAL(IB)
           DO JB=1,NBXX
             !GILO(JB,IB,LSPN)=FILO(JB,IB)
             GILO(JB,IB)=FILO(JB,IB)
           END DO
         END DO
       ELSE
C Copy values FILO <- GILO
       DO IB=1,NBXX
         !EVAL(IB)=EVSIC(IB,LSPN)
         EVAL(IB)=EVSIC(IB)
           DO JB=1,NBXX
             !FILO(JB,IB)=GILO(JB,IB,LSPN)
             FILO(JB,IB)=GILO(JB,IB)
           END DO
         END DO
       END IF
       deallocate(EVSIC, GILO)
       RETURN
       END
