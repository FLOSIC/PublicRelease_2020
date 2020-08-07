C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE POISSON2(NWRD,NPAIR,ND,MD,ALPHAV,A,BETAV,B,RHO)
C WRITTEN BY MARK R. PEDERSON 1996
       use common2,only : NSPN
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:57 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NPAIR, ND, MD, MXPR, IP
       REAL*8 :: SYMBOL , ALPHAV, A, BETAV, B, RHO, CLNRM, CLWRD,
     & TIME1, TIME2, TMNRM, TMWRD
       SAVE
       PARAMETER (MXPR=MXPOISS)
       COMMON/NEWTIMES/TMWRD,TMNRM,CLWRD,CLNRM
       LOGICAL NWRD
       DIMENSION ALPHAV(MXPR),BETAV(MXPR)
       DIMENSION A(3,MXPR),B(3,MXPR)
       DIMENSION RHO(10,10,MXPR)
       IF(NWRD)THEN
        CALL GTTIME(TIME1)
        DO IP=1,NPAIR,NSPN
         CALL POISSON1(NSPN,ND,MD,ALPHAV(IP),A(1,IP),BETAV(IP),B(1,IP),
     &                 RHO(1,1,IP))
        END DO
        CALL GTTIME(TIME2)
        TMWRD=TMWRD+TIME2-TIME1
        CLWRD=CLWRD+NPAIR
       ELSE
        CALL GTTIME(TIME1)
        CALL POISSON1(NPAIR,ND,MD,ALPHAV,A,BETAV,B,RHO)
        CALL GTTIME(TIME2)
        CLNRM=CLNRM+NPAIR
        TMNRM=TMNRM+TIME2-TIME1
       END IF
       RETURN
       END 
