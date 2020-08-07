C UTEP Electronic Structure Lab (2020)
C
C *********************************************************************
C
       SUBROUTINE RHCDRV(IID,NPTS,RTMP,DTMP,RHOCDR)
C
C DIRK POREZAG, MARCH 1998
C RHCDRV CALCULATES THE RADIAL DERIVATIVES OF THE CORE DENSITY 
C NEEDED FOR NONLINEAR CORE CORRECTIONS
C RETURN: RHOCDR= dRHO/dR * 1/R WHERE R= RTMP(IPTS)
C
       use common1,only : ISNLCC, RRADTAB, RHOCOR, NRADTAB, NLCC
       use common2,only : IFUIDT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:01 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IID, NPTS, IFU, IPTS
       REAL*8 :: RTMP , DTMP, RHOCDR
        SAVE
        LOGICAL     SKIP 
        DIMENSION   RTMP(*),DTMP(3,*),RHOCDR(*)
C
C SETUP
C
        IFU= IFUIDT(IID)
        SKIP= (ISNLCC .NE. 1) 
        IF (.NOT. SKIP) THEN
         SKIP= (NLCC(IFU) .NE. 1)
        END IF
        IF (SKIP) THEN
         DO IPTS=1,NPTS
          RHOCDR(IPTS)= 0.0D0
         END DO
         RETURN
        END IF
C
C INTERPOLATION
C
        CALL FINTPOL(8,NPTS,RTMP,0.1D0,NRADTAB(IFU),3,2,
     &               RRADTAB(1,IFU),RHOCOR(1,1,IFU),DTMP)
        DO IPTS=1,NPTS
         RHOCDR(IPTS)= DTMP(2,IPTS)/RTMP(IPTS)
        END DO
        RETURN
        END
