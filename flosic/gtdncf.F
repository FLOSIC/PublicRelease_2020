C UTEP Electronic Structure Lab (2020)
C
C HERE STARTS THE PART WHERE ALL BEAUTIFYING IS HOPELESS - DVP
C
       SUBROUTINE GTDNCF(GGA,ALP,A,BET,B,GAM,COEF,DOEF,DLT,C)
C WRITTEN BY MARK R PEDERSON 17-SEPTEMBER 1997
       IMPLICIT REAL*8 (A-H,O-Z)
       LOGICAL GGA
       DIMENSION A(3),B(3),C(3)
       DIMENSION COEF(35),GAM(10,10),DOEF(20,10)
           DLT=ALP+BET
           ARG=0.0D0
           DO I=1,3
           C(I)=(ALP*A(I)+BET*B(I))/DLT
           ARG=ARG+(A(I)-B(I))**2
           END DO
           ARG=ARG*ALP*BET/DLT

       CXMAX1=C(1)-A(1)
       CYMAY1=C(2)-A(2)
       CZMAZ1=C(3)-A(3)
       CXMAX2=CXMAX1*CXMAX1   
       CYMAY2=CYMAY1*CYMAY1   
       CZMAZ2=CZMAZ1*CZMAZ1   
       CXMBX1=C(1)-B(1)
       CYMBY1=C(2)-B(2)
       CZMBZ1=C(3)-B(3)
       CXMBX2=CXMBX1*CXMBX1   
       CYMBY2=CYMBY1*CYMBY1   
       CZMBZ2=CZMBZ1*CZMBZ1   
C      POLY( 1)=1.0D0                      
       COEF( 1)=0.0D0
     & +                                                GAM( 1, 1)
     & +             CXMBX1*                            GAM( 1, 2)
     & +                           CYMBY1*              GAM( 1, 3)
     & +                                         CZMBZ1*GAM( 1, 4)
     & +             CXMBX2*                            GAM( 1, 5)
     & +                           CYMBY2*              GAM( 1, 6)
     & +                                         CZMBZ2*GAM( 1, 7)
     & +             CXMBX1*       CYMBY1*              GAM( 1, 8)
     & +             CXMBX1*                     CZMBZ1*GAM( 1, 9)
     & +                           CYMBY1*       CZMBZ1*GAM( 1,10)
     & +      CXMAX1*                                   GAM( 2, 1)
     & +      CXMAX1*CXMBX1*                            GAM( 2, 2)
     & +      CXMAX1*              CYMBY1*              GAM( 2, 3)
     & +      CXMAX1*                            CZMBZ1*GAM( 2, 4)
     & +      CXMAX1*CXMBX2*                            GAM( 2, 5)
     & +      CXMAX1*              CYMBY2*              GAM( 2, 6)
     & +      CXMAX1*                            CZMBZ2*GAM( 2, 7)
     & +      CXMAX1*CXMBX1*       CYMBY1*              GAM( 2, 8)
     & +      CXMAX1*CXMBX1*                     CZMBZ1*GAM( 2, 9)
       COEF( 1)=COEF( 1)
     & +      CXMAX1*              CYMBY1*       CZMBZ1*GAM( 2,10)
     & +                    CYMAY1*                     GAM( 3, 1)
     & +             CXMBX1*CYMAY1*                     GAM( 3, 2)
     & +                    CYMAY1*CYMBY1*              GAM( 3, 3)
     & +                    CYMAY1*              CZMBZ1*GAM( 3, 4)
     & +             CXMBX2*CYMAY1*                     GAM( 3, 5)
     & +                    CYMAY1*CYMBY2*              GAM( 3, 6)
     & +                    CYMAY1*              CZMBZ2*GAM( 3, 7)
     & +             CXMBX1*CYMAY1*CYMBY1*              GAM( 3, 8)
     & +             CXMBX1*CYMAY1*              CZMBZ1*GAM( 3, 9)
     & +                    CYMAY1*CYMBY1*       CZMBZ1*GAM( 3,10)
     & +                                  CZMAZ1*       GAM( 4, 1)
     & +             CXMBX1*              CZMAZ1*       GAM( 4, 2)
     & +                           CYMBY1*CZMAZ1*       GAM( 4, 3)
     & +                                  CZMAZ1*CZMBZ1*GAM( 4, 4)
     & +             CXMBX2*              CZMAZ1*       GAM( 4, 5)
     & +                           CYMBY2*CZMAZ1*       GAM( 4, 6)
     & +                                  CZMAZ1*CZMBZ2*GAM( 4, 7)
     & +             CXMBX1*       CYMBY1*CZMAZ1*       GAM( 4, 8)
       COEF( 1)=COEF( 1)
     & +             CXMBX1*              CZMAZ1*CZMBZ1*GAM( 4, 9)
     & +                           CYMBY1*CZMAZ1*CZMBZ1*GAM( 4,10)
     & +      CXMAX2*                                   GAM( 5, 1)
     & +      CXMAX2*CXMBX1*                            GAM( 5, 2)
     & +      CXMAX2*              CYMBY1*              GAM( 5, 3)
     & +      CXMAX2*                            CZMBZ1*GAM( 5, 4)
     & +      CXMAX2*CXMBX2*                            GAM( 5, 5)
     & +      CXMAX2*              CYMBY2*              GAM( 5, 6)
     & +      CXMAX2*                            CZMBZ2*GAM( 5, 7)
     & +      CXMAX2*CXMBX1*       CYMBY1*              GAM( 5, 8)
     & +      CXMAX2*CXMBX1*                     CZMBZ1*GAM( 5, 9)
     & +      CXMAX2*              CYMBY1*       CZMBZ1*GAM( 5,10)
     & +                    CYMAY2*                     GAM( 6, 1)
     & +             CXMBX1*CYMAY2*                     GAM( 6, 2)
     & +                    CYMAY2*CYMBY1*              GAM( 6, 3)
     & +                    CYMAY2*              CZMBZ1*GAM( 6, 4)
     & +             CXMBX2*CYMAY2*                     GAM( 6, 5)
     & +                    CYMAY2*CYMBY2*              GAM( 6, 6)
     & +                    CYMAY2*              CZMBZ2*GAM( 6, 7)
       COEF( 1)=COEF( 1)
     & +             CXMBX1*CYMAY2*CYMBY1*              GAM( 6, 8)
     & +             CXMBX1*CYMAY2*              CZMBZ1*GAM( 6, 9)
     & +                    CYMAY2*CYMBY1*       CZMBZ1*GAM( 6,10)
     & +                                  CZMAZ2*       GAM( 7, 1)
     & +             CXMBX1*              CZMAZ2*       GAM( 7, 2)
     & +                           CYMBY1*CZMAZ2*       GAM( 7, 3)
     & +                                  CZMAZ2*CZMBZ1*GAM( 7, 4)
     & +             CXMBX2*              CZMAZ2*       GAM( 7, 5)
     & +                           CYMBY2*CZMAZ2*       GAM( 7, 6)
     & +                                  CZMAZ2*CZMBZ2*GAM( 7, 7)
     & +             CXMBX1*       CYMBY1*CZMAZ2*       GAM( 7, 8)
     & +             CXMBX1*              CZMAZ2*CZMBZ1*GAM( 7, 9)
     & +                           CYMBY1*CZMAZ2*CZMBZ1*GAM( 7,10)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 1)
     & +      CXMAX1*CXMBX1*CYMAY1*                     GAM( 8, 2)
     & +      CXMAX1*       CYMAY1*CYMBY1*              GAM( 8, 3)
     & +      CXMAX1*       CYMAY1*              CZMBZ1*GAM( 8, 4)
     & +      CXMAX1*CXMBX2*CYMAY1*                     GAM( 8, 5)
     & +      CXMAX1*       CYMAY1*CYMBY2*              GAM( 8, 6)
       COEF( 1)=COEF( 1)
     & +      CXMAX1*       CYMAY1*              CZMBZ2*GAM( 8, 7)
     & +      CXMAX1*CXMBX1*CYMAY1*CYMBY1*              GAM( 8, 8)
     & +      CXMAX1*CXMBX1*CYMAY1*              CZMBZ1*GAM( 8, 9)
     & +      CXMAX1*       CYMAY1*CYMBY1*       CZMBZ1*GAM( 8,10)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 1)
     & +      CXMAX1*CXMBX1*              CZMAZ1*       GAM( 9, 2)
     & +      CXMAX1*              CYMBY1*CZMAZ1*       GAM( 9, 3)
     & +      CXMAX1*                     CZMAZ1*CZMBZ1*GAM( 9, 4)
     & +      CXMAX1*CXMBX2*              CZMAZ1*       GAM( 9, 5)
     & +      CXMAX1*              CYMBY2*CZMAZ1*       GAM( 9, 6)
     & +      CXMAX1*                     CZMAZ1*CZMBZ2*GAM( 9, 7)
     & +      CXMAX1*CXMBX1*       CYMBY1*CZMAZ1*       GAM( 9, 8)
     & +      CXMAX1*CXMBX1*              CZMAZ1*CZMBZ1*GAM( 9, 9)
     & +      CXMAX1*              CYMBY1*CZMAZ1*CZMBZ1*GAM( 9,10)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 1)
     & +             CXMBX1*CYMAY1*       CZMAZ1*       GAM(10, 2)
     & +                    CYMAY1*CYMBY1*CZMAZ1*       GAM(10, 3)
     & +                    CYMAY1*       CZMAZ1*CZMBZ1*GAM(10, 4)
     & +             CXMBX2*CYMAY1*       CZMAZ1*       GAM(10, 5)
       COEF( 1)=COEF( 1)
     & +                    CYMAY1*CYMBY2*CZMAZ1*       GAM(10, 6)
     & +                    CYMAY1*       CZMAZ1*CZMBZ2*GAM(10, 7)
     & +             CXMBX1*CYMAY1*CYMBY1*CZMAZ1*       GAM(10, 8)
     & +             CXMBX1*CYMAY1*       CZMAZ1*CZMBZ1*GAM(10, 9)
     & +                    CYMAY1*CYMBY1*CZMAZ1*CZMBZ1*GAM(10,10)
  
C      POLY( 2)=1.0D0              *ZZMCZ1 
       COEF( 2)=0.0D0
     & +                                                GAM( 1, 4)
     & +2.0D0*                                   CZMBZ1*GAM( 1, 7)
     & +             CXMBX1*                            GAM( 1, 9)
     & +                           CYMBY1*              GAM( 1,10)
     & +      CXMAX1*                                   GAM( 2, 4)
     & +2.0D0*CXMAX1*                            CZMBZ1*GAM( 2, 7)
     & +      CXMAX1*CXMBX1*                            GAM( 2, 9)
     & +      CXMAX1*              CYMBY1*              GAM( 2,10)
     & +                    CYMAY1*                     GAM( 3, 4)
     & +2.0D0*              CYMAY1*              CZMBZ1*GAM( 3, 7)
     & +             CXMBX1*CYMAY1*                     GAM( 3, 9)
     & +                    CYMAY1*CYMBY1*              GAM( 3,10)
     & +                                                GAM( 4, 1)
     & +             CXMBX1*                            GAM( 4, 2)
     & +                           CYMBY1*              GAM( 4, 3)
     & +                                  CZMAZ1*       GAM( 4, 4)
     & +                                         CZMBZ1*GAM( 4, 4)
     & +             CXMBX2*                            GAM( 4, 5)
     & +                           CYMBY2*              GAM( 4, 6)
       COEF( 2)=COEF( 2)
     & +2.0D0*                            CZMAZ1*CZMBZ1*GAM( 4, 7)
     & +                                         CZMBZ2*GAM( 4, 7)
     & +             CXMBX1*       CYMBY1*              GAM( 4, 8)
     & +             CXMBX1*              CZMAZ1*       GAM( 4, 9)
     & +             CXMBX1*                     CZMBZ1*GAM( 4, 9)
     & +                           CYMBY1*CZMAZ1*       GAM( 4,10)
     & +                           CYMBY1*       CZMBZ1*GAM( 4,10)
     & +      CXMAX2*                                   GAM( 5, 4)
     & +2.0D0*CXMAX2*                            CZMBZ1*GAM( 5, 7)
     & +      CXMAX2*CXMBX1*                            GAM( 5, 9)
     & +      CXMAX2*              CYMBY1*              GAM( 5,10)
     & +                    CYMAY2*                     GAM( 6, 4)
     & +2.0D0*              CYMAY2*              CZMBZ1*GAM( 6, 7)
     & +             CXMBX1*CYMAY2*                     GAM( 6, 9)
     & +                    CYMAY2*CYMBY1*              GAM( 6,10)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 1)
     & +2.0D0*       CXMBX1*              CZMAZ1*       GAM( 7, 2)
     & +2.0D0*                     CYMBY1*CZMAZ1*       GAM( 7, 3)
     & +                                  CZMAZ2*       GAM( 7, 4)
       COEF( 2)=COEF( 2)
     & +2.0D0*                            CZMAZ1*CZMBZ1*GAM( 7, 4)
     & +2.0D0*       CXMBX2*              CZMAZ1*       GAM( 7, 5)
     & +2.0D0*                     CYMBY2*CZMAZ1*       GAM( 7, 6)
     & +2.0D0*                            CZMAZ2*CZMBZ1*GAM( 7, 7)
     & +2.0D0*                            CZMAZ1*CZMBZ2*GAM( 7, 7)
     & +2.0D0*       CXMBX1*       CYMBY1*CZMAZ1*       GAM( 7, 8)
     & +             CXMBX1*              CZMAZ2*       GAM( 7, 9)
     & +2.0D0*       CXMBX1*              CZMAZ1*CZMBZ1*GAM( 7, 9)
     & +                           CYMBY1*CZMAZ2*       GAM( 7,10)
     & +2.0D0*                     CYMBY1*CZMAZ1*CZMBZ1*GAM( 7,10)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 4)
     & +2.0D0*CXMAX1*       CYMAY1*              CZMBZ1*GAM( 8, 7)
     & +      CXMAX1*CXMBX1*CYMAY1*                     GAM( 8, 9)
     & +      CXMAX1*       CYMAY1*CYMBY1*              GAM( 8,10)
     & +      CXMAX1*                                   GAM( 9, 1)
     & +      CXMAX1*CXMBX1*                            GAM( 9, 2)
     & +      CXMAX1*              CYMBY1*              GAM( 9, 3)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 4)
     & +      CXMAX1*                            CZMBZ1*GAM( 9, 4)
       COEF( 2)=COEF( 2)
     & +      CXMAX1*CXMBX2*                            GAM( 9, 5)
     & +      CXMAX1*              CYMBY2*              GAM( 9, 6)
     & +2.0D0*CXMAX1*                     CZMAZ1*CZMBZ1*GAM( 9, 7)
     & +      CXMAX1*                            CZMBZ2*GAM( 9, 7)
     & +      CXMAX1*CXMBX1*       CYMBY1*              GAM( 9, 8)
     & +      CXMAX1*CXMBX1*              CZMAZ1*       GAM( 9, 9)
     & +      CXMAX1*CXMBX1*                     CZMBZ1*GAM( 9, 9)
     & +      CXMAX1*              CYMBY1*CZMAZ1*       GAM( 9,10)
     & +      CXMAX1*              CYMBY1*       CZMBZ1*GAM( 9,10)
     & +                    CYMAY1*                     GAM(10, 1)
     & +             CXMBX1*CYMAY1*                     GAM(10, 2)
     & +                    CYMAY1*CYMBY1*              GAM(10, 3)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 4)
     & +                    CYMAY1*              CZMBZ1*GAM(10, 4)
     & +             CXMBX2*CYMAY1*                     GAM(10, 5)
     & +                    CYMAY1*CYMBY2*              GAM(10, 6)
     & +2.0D0*              CYMAY1*       CZMAZ1*CZMBZ1*GAM(10, 7)
     & +                    CYMAY1*              CZMBZ2*GAM(10, 7)
     & +             CXMBX1*CYMAY1*CYMBY1*              GAM(10, 8)
       COEF( 2)=COEF( 2)
     & +             CXMBX1*CYMAY1*       CZMAZ1*       GAM(10, 9)
     & +             CXMBX1*CYMAY1*              CZMBZ1*GAM(10, 9)
     & +                    CYMAY1*CYMBY1*CZMAZ1*       GAM(10,10)
     & +                    CYMAY1*CYMBY1*       CZMBZ1*GAM(10,10)
  
C      POLY( 3)=1.0D0       *YYMCY1        
       COEF( 3)=0.0D0
     & +                                                GAM( 1, 3)
     & +2.0D0*                     CYMBY1*              GAM( 1, 6)
     & +             CXMBX1*                            GAM( 1, 8)
     & +                                         CZMBZ1*GAM( 1,10)
     & +      CXMAX1*                                   GAM( 2, 3)
     & +2.0D0*CXMAX1*              CYMBY1*              GAM( 2, 6)
     & +      CXMAX1*CXMBX1*                            GAM( 2, 8)
     & +      CXMAX1*                            CZMBZ1*GAM( 2,10)
     & +                                                GAM( 3, 1)
     & +             CXMBX1*                            GAM( 3, 2)
     & +                    CYMAY1*                     GAM( 3, 3)
     & +                           CYMBY1*              GAM( 3, 3)
     & +                                         CZMBZ1*GAM( 3, 4)
     & +             CXMBX2*                            GAM( 3, 5)
     & +2.0D0*              CYMAY1*CYMBY1*              GAM( 3, 6)
     & +                           CYMBY2*              GAM( 3, 6)
     & +                                         CZMBZ2*GAM( 3, 7)
     & +             CXMBX1*CYMAY1*                     GAM( 3, 8)
     & +             CXMBX1*       CYMBY1*              GAM( 3, 8)
       COEF( 3)=COEF( 3)
     & +             CXMBX1*                     CZMBZ1*GAM( 3, 9)
     & +                    CYMAY1*              CZMBZ1*GAM( 3,10)
     & +                           CYMBY1*       CZMBZ1*GAM( 3,10)
     & +                                  CZMAZ1*       GAM( 4, 3)
     & +2.0D0*                     CYMBY1*CZMAZ1*       GAM( 4, 6)
     & +             CXMBX1*              CZMAZ1*       GAM( 4, 8)
     & +                                  CZMAZ1*CZMBZ1*GAM( 4,10)
     & +      CXMAX2*                                   GAM( 5, 3)
     & +2.0D0*CXMAX2*              CYMBY1*              GAM( 5, 6)
     & +      CXMAX2*CXMBX1*                            GAM( 5, 8)
     & +      CXMAX2*                            CZMBZ1*GAM( 5,10)
     & +2.0D0*              CYMAY1*                     GAM( 6, 1)
     & +2.0D0*       CXMBX1*CYMAY1*                     GAM( 6, 2)
     & +                    CYMAY2*                     GAM( 6, 3)
     & +2.0D0*              CYMAY1*CYMBY1*              GAM( 6, 3)
     & +2.0D0*              CYMAY1*              CZMBZ1*GAM( 6, 4)
     & +2.0D0*       CXMBX2*CYMAY1*                     GAM( 6, 5)
     & +2.0D0*              CYMAY2*CYMBY1*              GAM( 6, 6)
     & +2.0D0*              CYMAY1*CYMBY2*              GAM( 6, 6)
       COEF( 3)=COEF( 3)
     & +2.0D0*              CYMAY1*              CZMBZ2*GAM( 6, 7)
     & +             CXMBX1*CYMAY2*                     GAM( 6, 8)
     & +2.0D0*       CXMBX1*CYMAY1*CYMBY1*              GAM( 6, 8)
     & +2.0D0*       CXMBX1*CYMAY1*              CZMBZ1*GAM( 6, 9)
     & +                    CYMAY2*              CZMBZ1*GAM( 6,10)
     & +2.0D0*              CYMAY1*CYMBY1*       CZMBZ1*GAM( 6,10)
     & +                                  CZMAZ2*       GAM( 7, 3)
     & +2.0D0*                     CYMBY1*CZMAZ2*       GAM( 7, 6)
     & +             CXMBX1*              CZMAZ2*       GAM( 7, 8)
     & +                                  CZMAZ2*CZMBZ1*GAM( 7,10)
     & +      CXMAX1*                                   GAM( 8, 1)
     & +      CXMAX1*CXMBX1*                            GAM( 8, 2)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 3)
     & +      CXMAX1*              CYMBY1*              GAM( 8, 3)
     & +      CXMAX1*                            CZMBZ1*GAM( 8, 4)
     & +      CXMAX1*CXMBX2*                            GAM( 8, 5)
     & +2.0D0*CXMAX1*       CYMAY1*CYMBY1*              GAM( 8, 6)
     & +      CXMAX1*              CYMBY2*              GAM( 8, 6)
     & +      CXMAX1*                            CZMBZ2*GAM( 8, 7)
       COEF( 3)=COEF( 3)
     & +      CXMAX1*CXMBX1*CYMAY1*                     GAM( 8, 8)
     & +      CXMAX1*CXMBX1*       CYMBY1*              GAM( 8, 8)
     & +      CXMAX1*CXMBX1*                     CZMBZ1*GAM( 8, 9)
     & +      CXMAX1*       CYMAY1*              CZMBZ1*GAM( 8,10)
     & +      CXMAX1*              CYMBY1*       CZMBZ1*GAM( 8,10)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 3)
     & +2.0D0*CXMAX1*              CYMBY1*CZMAZ1*       GAM( 9, 6)
     & +      CXMAX1*CXMBX1*              CZMAZ1*       GAM( 9, 8)
     & +      CXMAX1*                     CZMAZ1*CZMBZ1*GAM( 9,10)
     & +                                  CZMAZ1*       GAM(10, 1)
     & +             CXMBX1*              CZMAZ1*       GAM(10, 2)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 3)
     & +                           CYMBY1*CZMAZ1*       GAM(10, 3)
     & +                                  CZMAZ1*CZMBZ1*GAM(10, 4)
     & +             CXMBX2*              CZMAZ1*       GAM(10, 5)
     & +2.0D0*              CYMAY1*CYMBY1*CZMAZ1*       GAM(10, 6)
     & +                           CYMBY2*CZMAZ1*       GAM(10, 6)
     & +                                  CZMAZ1*CZMBZ2*GAM(10, 7)
     & +             CXMBX1*CYMAY1*       CZMAZ1*       GAM(10, 8)
       COEF( 3)=COEF( 3)
     & +             CXMBX1*       CYMBY1*CZMAZ1*       GAM(10, 8)
     & +             CXMBX1*              CZMAZ1*CZMBZ1*GAM(10, 9)
     & +                    CYMAY1*       CZMAZ1*CZMBZ1*GAM(10,10)
     & +                           CYMBY1*CZMAZ1*CZMBZ1*GAM(10,10)
  
C      POLY( 4)=      XXMCX1               
       COEF( 4)=0.0D0
     & +                                                GAM( 1, 2)
     & +2.0D0*       CXMBX1*                            GAM( 1, 5)
     & +                           CYMBY1*              GAM( 1, 8)
     & +                                         CZMBZ1*GAM( 1, 9)
     & +                                                GAM( 2, 1)
     & +      CXMAX1*                                   GAM( 2, 2)
     & +             CXMBX1*                            GAM( 2, 2)
     & +                           CYMBY1*              GAM( 2, 3)
     & +                                         CZMBZ1*GAM( 2, 4)
     & +2.0D0*CXMAX1*CXMBX1*                            GAM( 2, 5)
     & +             CXMBX2*                            GAM( 2, 5)
     & +                           CYMBY2*              GAM( 2, 6)
     & +                                         CZMBZ2*GAM( 2, 7)
     & +      CXMAX1*              CYMBY1*              GAM( 2, 8)
     & +             CXMBX1*       CYMBY1*              GAM( 2, 8)
     & +      CXMAX1*                            CZMBZ1*GAM( 2, 9)
     & +             CXMBX1*                     CZMBZ1*GAM( 2, 9)
     & +                           CYMBY1*       CZMBZ1*GAM( 2,10)
     & +                    CYMAY1*                     GAM( 3, 2)
       COEF( 4)=COEF( 4)
     & +2.0D0*       CXMBX1*CYMAY1*                     GAM( 3, 5)
     & +                    CYMAY1*CYMBY1*              GAM( 3, 8)
     & +                    CYMAY1*              CZMBZ1*GAM( 3, 9)
     & +                                  CZMAZ1*       GAM( 4, 2)
     & +2.0D0*       CXMBX1*              CZMAZ1*       GAM( 4, 5)
     & +                           CYMBY1*CZMAZ1*       GAM( 4, 8)
     & +                                  CZMAZ1*CZMBZ1*GAM( 4, 9)
     & +2.0D0*CXMAX1*                                   GAM( 5, 1)
     & +      CXMAX2*                                   GAM( 5, 2)
     & +2.0D0*CXMAX1*CXMBX1*                            GAM( 5, 2)
     & +2.0D0*CXMAX1*              CYMBY1*              GAM( 5, 3)
     & +2.0D0*CXMAX1*                            CZMBZ1*GAM( 5, 4)
     & +2.0D0*CXMAX2*CXMBX1*                            GAM( 5, 5)
     & +2.0D0*CXMAX1*CXMBX2*                            GAM( 5, 5)
     & +2.0D0*CXMAX1*              CYMBY2*              GAM( 5, 6)
     & +2.0D0*CXMAX1*                            CZMBZ2*GAM( 5, 7)
     & +      CXMAX2*              CYMBY1*              GAM( 5, 8)
     & +2.0D0*CXMAX1*CXMBX1*       CYMBY1*              GAM( 5, 8)
     & +      CXMAX2*                            CZMBZ1*GAM( 5, 9)
       COEF( 4)=COEF( 4)
     & +2.0D0*CXMAX1*CXMBX1*                     CZMBZ1*GAM( 5, 9)
     & +2.0D0*CXMAX1*              CYMBY1*       CZMBZ1*GAM( 5,10)
     & +                    CYMAY2*                     GAM( 6, 2)
     & +2.0D0*       CXMBX1*CYMAY2*                     GAM( 6, 5)
     & +                    CYMAY2*CYMBY1*              GAM( 6, 8)
     & +                    CYMAY2*              CZMBZ1*GAM( 6, 9)
     & +                                  CZMAZ2*       GAM( 7, 2)
     & +2.0D0*       CXMBX1*              CZMAZ2*       GAM( 7, 5)
     & +                           CYMBY1*CZMAZ2*       GAM( 7, 8)
     & +                                  CZMAZ2*CZMBZ1*GAM( 7, 9)
     & +                    CYMAY1*                     GAM( 8, 1)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 2)
     & +             CXMBX1*CYMAY1*                     GAM( 8, 2)
     & +                    CYMAY1*CYMBY1*              GAM( 8, 3)
     & +                    CYMAY1*              CZMBZ1*GAM( 8, 4)
     & +2.0D0*CXMAX1*CXMBX1*CYMAY1*                     GAM( 8, 5)
     & +             CXMBX2*CYMAY1*                     GAM( 8, 5)
     & +                    CYMAY1*CYMBY2*              GAM( 8, 6)
     & +                    CYMAY1*              CZMBZ2*GAM( 8, 7)
       COEF( 4)=COEF( 4)
     & +      CXMAX1*       CYMAY1*CYMBY1*              GAM( 8, 8)
     & +             CXMBX1*CYMAY1*CYMBY1*              GAM( 8, 8)
     & +      CXMAX1*       CYMAY1*              CZMBZ1*GAM( 8, 9)
     & +             CXMBX1*CYMAY1*              CZMBZ1*GAM( 8, 9)
     & +                    CYMAY1*CYMBY1*       CZMBZ1*GAM( 8,10)
     & +                                  CZMAZ1*       GAM( 9, 1)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 2)
     & +             CXMBX1*              CZMAZ1*       GAM( 9, 2)
     & +                           CYMBY1*CZMAZ1*       GAM( 9, 3)
     & +                                  CZMAZ1*CZMBZ1*GAM( 9, 4)
     & +2.0D0*CXMAX1*CXMBX1*              CZMAZ1*       GAM( 9, 5)
     & +             CXMBX2*              CZMAZ1*       GAM( 9, 5)
     & +                           CYMBY2*CZMAZ1*       GAM( 9, 6)
     & +                                  CZMAZ1*CZMBZ2*GAM( 9, 7)
     & +      CXMAX1*              CYMBY1*CZMAZ1*       GAM( 9, 8)
     & +             CXMBX1*       CYMBY1*CZMAZ1*       GAM( 9, 8)
     & +      CXMAX1*                     CZMAZ1*CZMBZ1*GAM( 9, 9)
     & +             CXMBX1*              CZMAZ1*CZMBZ1*GAM( 9, 9)
     & +                           CYMBY1*CZMAZ1*CZMBZ1*GAM( 9,10)
       COEF( 4)=COEF( 4)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 2)
     & +2.0D0*       CXMBX1*CYMAY1*       CZMAZ1*       GAM(10, 5)
     & +                    CYMAY1*CYMBY1*CZMAZ1*       GAM(10, 8)
     & +                    CYMAY1*       CZMAZ1*CZMBZ1*GAM(10, 9)
  
C      POLY( 5)=1.0D0              *ZZMCZ2 
       COEF( 5)=0.0D0
     & +                                                GAM( 1, 7)
     & +      CXMAX1*                                   GAM( 2, 7)
     & +                    CYMAY1*                     GAM( 3, 7)
     & +                                                GAM( 4, 4)
     & +                                  CZMAZ1*       GAM( 4, 7)
     & +2.0D0*                                   CZMBZ1*GAM( 4, 7)
     & +             CXMBX1*                            GAM( 4, 9)
     & +                           CYMBY1*              GAM( 4,10)
     & +      CXMAX2*                                   GAM( 5, 7)
     & +                    CYMAY2*                     GAM( 6, 7)
     & +                                                GAM( 7, 1)
     & +             CXMBX1*                            GAM( 7, 2)
     & +                           CYMBY1*              GAM( 7, 3)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 4)
     & +                                         CZMBZ1*GAM( 7, 4)
     & +             CXMBX2*                            GAM( 7, 5)
     & +                           CYMBY2*              GAM( 7, 6)
     & +                                  CZMAZ2*       GAM( 7, 7)
     & +4.0D0*                            CZMAZ1*CZMBZ1*GAM( 7, 7)
       COEF( 5)=COEF( 5)
     & +                                         CZMBZ2*GAM( 7, 7)
     & +             CXMBX1*       CYMBY1*              GAM( 7, 8)
     & +2.0D0*       CXMBX1*              CZMAZ1*       GAM( 7, 9)
     & +             CXMBX1*                     CZMBZ1*GAM( 7, 9)
     & +2.0D0*                     CYMBY1*CZMAZ1*       GAM( 7,10)
     & +                           CYMBY1*       CZMBZ1*GAM( 7,10)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 7)
     & +      CXMAX1*                                   GAM( 9, 4)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 7)
     & +2.0D0*CXMAX1*                            CZMBZ1*GAM( 9, 7)
     & +      CXMAX1*CXMBX1*                            GAM( 9, 9)
     & +      CXMAX1*              CYMBY1*              GAM( 9,10)
     & +                    CYMAY1*                     GAM(10, 4)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 7)
     & +2.0D0*              CYMAY1*              CZMBZ1*GAM(10, 7)
     & +             CXMBX1*CYMAY1*                     GAM(10, 9)
     & +                    CYMAY1*CYMBY1*              GAM(10,10)
  
C      POLY( 6)=1.0D0       *YYMCY1*ZZMCZ1 
       COEF( 6)=0.0D0
     & +                                                GAM( 1,10)
     & +      CXMAX1*                                   GAM( 2,10)
     & +                                                GAM( 3, 4)
     & +2.0D0*                                   CZMBZ1*GAM( 3, 7)
     & +             CXMBX1*                            GAM( 3, 9)
     & +                    CYMAY1*                     GAM( 3,10)
     & +                           CYMBY1*              GAM( 3,10)
     & +                                                GAM( 4, 3)
     & +2.0D0*                     CYMBY1*              GAM( 4, 6)
     & +             CXMBX1*                            GAM( 4, 8)
     & +                                  CZMAZ1*       GAM( 4,10)
     & +                                         CZMBZ1*GAM( 4,10)
     & +      CXMAX2*                                   GAM( 5,10)
     & +2.0D0*              CYMAY1*                     GAM( 6, 4)
     & +4.0D0*              CYMAY1*              CZMBZ1*GAM( 6, 7)
     & +2.0D0*       CXMBX1*CYMAY1*                     GAM( 6, 9)
     & +                    CYMAY2*                     GAM( 6,10)
     & +2.0D0*              CYMAY1*CYMBY1*              GAM( 6,10)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 3)
       COEF( 6)=COEF( 6)
     & +4.0D0*                     CYMBY1*CZMAZ1*       GAM( 7, 6)
     & +2.0D0*       CXMBX1*              CZMAZ1*       GAM( 7, 8)
     & +                                  CZMAZ2*       GAM( 7,10)
     & +2.0D0*                            CZMAZ1*CZMBZ1*GAM( 7,10)
     & +      CXMAX1*                                   GAM( 8, 4)
     & +2.0D0*CXMAX1*                            CZMBZ1*GAM( 8, 7)
     & +      CXMAX1*CXMBX1*                            GAM( 8, 9)
     & +      CXMAX1*       CYMAY1*                     GAM( 8,10)
     & +      CXMAX1*              CYMBY1*              GAM( 8,10)
     & +      CXMAX1*                                   GAM( 9, 3)
     & +2.0D0*CXMAX1*              CYMBY1*              GAM( 9, 6)
     & +      CXMAX1*CXMBX1*                            GAM( 9, 8)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9,10)
     & +      CXMAX1*                            CZMBZ1*GAM( 9,10)
     & +                                                GAM(10, 1)
     & +             CXMBX1*                            GAM(10, 2)
     & +                    CYMAY1*                     GAM(10, 3)
     & +                           CYMBY1*              GAM(10, 3)
     & +                                  CZMAZ1*       GAM(10, 4)
       COEF( 6)=COEF( 6)
     & +                                         CZMBZ1*GAM(10, 4)
     & +             CXMBX2*                            GAM(10, 5)
     & +2.0D0*              CYMAY1*CYMBY1*              GAM(10, 6)
     & +                           CYMBY2*              GAM(10, 6)
     & +2.0D0*                            CZMAZ1*CZMBZ1*GAM(10, 7)
     & +                                         CZMBZ2*GAM(10, 7)
     & +             CXMBX1*CYMAY1*                     GAM(10, 8)
     & +             CXMBX1*       CYMBY1*              GAM(10, 8)
     & +             CXMBX1*              CZMAZ1*       GAM(10, 9)
     & +             CXMBX1*                     CZMBZ1*GAM(10, 9)
     & +                    CYMAY1*       CZMAZ1*       GAM(10,10)
     & +                    CYMAY1*              CZMBZ1*GAM(10,10)
     & +                           CYMBY1*CZMAZ1*       GAM(10,10)
     & +                           CYMBY1*       CZMBZ1*GAM(10,10)
  
C      POLY( 7)=1.0D0       *YYMCY2        
       COEF( 7)=0.0D0
     & +                                                GAM( 1, 6)
     & +      CXMAX1*                                   GAM( 2, 6)
     & +                                                GAM( 3, 3)
     & +                    CYMAY1*                     GAM( 3, 6)
     & +2.0D0*                     CYMBY1*              GAM( 3, 6)
     & +             CXMBX1*                            GAM( 3, 8)
     & +                                         CZMBZ1*GAM( 3,10)
     & +                                  CZMAZ1*       GAM( 4, 6)
     & +      CXMAX2*                                   GAM( 5, 6)
     & +                                                GAM( 6, 1)
     & +             CXMBX1*                            GAM( 6, 2)
     & +2.0D0*              CYMAY1*                     GAM( 6, 3)
     & +                           CYMBY1*              GAM( 6, 3)
     & +                                         CZMBZ1*GAM( 6, 4)
     & +             CXMBX2*                            GAM( 6, 5)
     & +                    CYMAY2*                     GAM( 6, 6)
     & +4.0D0*              CYMAY1*CYMBY1*              GAM( 6, 6)
     & +                           CYMBY2*              GAM( 6, 6)
     & +                                         CZMBZ2*GAM( 6, 7)
       COEF( 7)=COEF( 7)
     & +2.0D0*       CXMBX1*CYMAY1*                     GAM( 6, 8)
     & +             CXMBX1*       CYMBY1*              GAM( 6, 8)
     & +             CXMBX1*                     CZMBZ1*GAM( 6, 9)
     & +2.0D0*              CYMAY1*              CZMBZ1*GAM( 6,10)
     & +                           CYMBY1*       CZMBZ1*GAM( 6,10)
     & +                                  CZMAZ2*       GAM( 7, 6)
     & +      CXMAX1*                                   GAM( 8, 3)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 6)
     & +2.0D0*CXMAX1*              CYMBY1*              GAM( 8, 6)
     & +      CXMAX1*CXMBX1*                            GAM( 8, 8)
     & +      CXMAX1*                            CZMBZ1*GAM( 8,10)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 6)
     & +                                  CZMAZ1*       GAM(10, 3)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 6)
     & +2.0D0*                     CYMBY1*CZMAZ1*       GAM(10, 6)
     & +             CXMBX1*              CZMAZ1*       GAM(10, 8)
     & +                                  CZMAZ1*CZMBZ1*GAM(10,10)
  
C      POLY( 8)=      XXMCX1       *ZZMCZ1 
       COEF( 8)=0.0D0
     & +                                                GAM( 1, 9)
     & +                                                GAM( 2, 4)
     & +2.0D0*                                   CZMBZ1*GAM( 2, 7)
     & +      CXMAX1*                                   GAM( 2, 9)
     & +             CXMBX1*                            GAM( 2, 9)
     & +                           CYMBY1*              GAM( 2,10)
     & +                    CYMAY1*                     GAM( 3, 9)
     & +                                                GAM( 4, 2)
     & +2.0D0*       CXMBX1*                            GAM( 4, 5)
     & +                           CYMBY1*              GAM( 4, 8)
     & +                                  CZMAZ1*       GAM( 4, 9)
     & +                                         CZMBZ1*GAM( 4, 9)
     & +2.0D0*CXMAX1*                                   GAM( 5, 4)
     & +4.0D0*CXMAX1*                            CZMBZ1*GAM( 5, 7)
     & +      CXMAX2*                                   GAM( 5, 9)
     & +2.0D0*CXMAX1*CXMBX1*                            GAM( 5, 9)
     & +2.0D0*CXMAX1*              CYMBY1*              GAM( 5,10)
     & +                    CYMAY2*                     GAM( 6, 9)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 2)
       COEF( 8)=COEF( 8)
     & +4.0D0*       CXMBX1*              CZMAZ1*       GAM( 7, 5)
     & +2.0D0*                     CYMBY1*CZMAZ1*       GAM( 7, 8)
     & +                                  CZMAZ2*       GAM( 7, 9)
     & +2.0D0*                            CZMAZ1*CZMBZ1*GAM( 7, 9)
     & +                    CYMAY1*                     GAM( 8, 4)
     & +2.0D0*              CYMAY1*              CZMBZ1*GAM( 8, 7)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 9)
     & +             CXMBX1*CYMAY1*                     GAM( 8, 9)
     & +                    CYMAY1*CYMBY1*              GAM( 8,10)
     & +                                                GAM( 9, 1)
     & +      CXMAX1*                                   GAM( 9, 2)
     & +             CXMBX1*                            GAM( 9, 2)
     & +                           CYMBY1*              GAM( 9, 3)
     & +                                  CZMAZ1*       GAM( 9, 4)
     & +                                         CZMBZ1*GAM( 9, 4)
     & +2.0D0*CXMAX1*CXMBX1*                            GAM( 9, 5)
     & +             CXMBX2*                            GAM( 9, 5)
     & +                           CYMBY2*              GAM( 9, 6)
     & +2.0D0*                            CZMAZ1*CZMBZ1*GAM( 9, 7)
       COEF( 8)=COEF( 8)
     & +                                         CZMBZ2*GAM( 9, 7)
     & +      CXMAX1*              CYMBY1*              GAM( 9, 8)
     & +             CXMBX1*       CYMBY1*              GAM( 9, 8)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 9)
     & +      CXMAX1*                            CZMBZ1*GAM( 9, 9)
     & +             CXMBX1*              CZMAZ1*       GAM( 9, 9)
     & +             CXMBX1*                     CZMBZ1*GAM( 9, 9)
     & +                           CYMBY1*CZMAZ1*       GAM( 9,10)
     & +                           CYMBY1*       CZMBZ1*GAM( 9,10)
     & +                    CYMAY1*                     GAM(10, 2)
     & +2.0D0*       CXMBX1*CYMAY1*                     GAM(10, 5)
     & +                    CYMAY1*CYMBY1*              GAM(10, 8)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 9)
     & +                    CYMAY1*              CZMBZ1*GAM(10, 9)
  
C      POLY( 9)=      XXMCX1*YYMCY1        
       COEF( 9)=0.0D0
     & +                                                GAM( 1, 8)
     & +                                                GAM( 2, 3)
     & +2.0D0*                     CYMBY1*              GAM( 2, 6)
     & +      CXMAX1*                                   GAM( 2, 8)
     & +             CXMBX1*                            GAM( 2, 8)
     & +                                         CZMBZ1*GAM( 2,10)
     & +                                                GAM( 3, 2)
     & +2.0D0*       CXMBX1*                            GAM( 3, 5)
     & +                    CYMAY1*                     GAM( 3, 8)
     & +                           CYMBY1*              GAM( 3, 8)
     & +                                         CZMBZ1*GAM( 3, 9)
     & +                                  CZMAZ1*       GAM( 4, 8)
     & +2.0D0*CXMAX1*                                   GAM( 5, 3)
     & +4.0D0*CXMAX1*              CYMBY1*              GAM( 5, 6)
     & +      CXMAX2*                                   GAM( 5, 8)
     & +2.0D0*CXMAX1*CXMBX1*                            GAM( 5, 8)
     & +2.0D0*CXMAX1*                            CZMBZ1*GAM( 5,10)
     & +2.0D0*              CYMAY1*                     GAM( 6, 2)
     & +4.0D0*       CXMBX1*CYMAY1*                     GAM( 6, 5)
       COEF( 9)=COEF( 9)
     & +                    CYMAY2*                     GAM( 6, 8)
     & +2.0D0*              CYMAY1*CYMBY1*              GAM( 6, 8)
     & +2.0D0*              CYMAY1*              CZMBZ1*GAM( 6, 9)
     & +                                  CZMAZ2*       GAM( 7, 8)
     & +                                                GAM( 8, 1)
     & +      CXMAX1*                                   GAM( 8, 2)
     & +             CXMBX1*                            GAM( 8, 2)
     & +                    CYMAY1*                     GAM( 8, 3)
     & +                           CYMBY1*              GAM( 8, 3)
     & +                                         CZMBZ1*GAM( 8, 4)
     & +2.0D0*CXMAX1*CXMBX1*                            GAM( 8, 5)
     & +             CXMBX2*                            GAM( 8, 5)
     & +2.0D0*              CYMAY1*CYMBY1*              GAM( 8, 6)
     & +                           CYMBY2*              GAM( 8, 6)
     & +                                         CZMBZ2*GAM( 8, 7)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 8)
     & +      CXMAX1*              CYMBY1*              GAM( 8, 8)
     & +             CXMBX1*CYMAY1*                     GAM( 8, 8)
     & +             CXMBX1*       CYMBY1*              GAM( 8, 8)
       COEF( 9)=COEF( 9)
     & +      CXMAX1*                            CZMBZ1*GAM( 8, 9)
     & +             CXMBX1*                     CZMBZ1*GAM( 8, 9)
     & +                    CYMAY1*              CZMBZ1*GAM( 8,10)
     & +                           CYMBY1*       CZMBZ1*GAM( 8,10)
     & +                                  CZMAZ1*       GAM( 9, 3)
     & +2.0D0*                     CYMBY1*CZMAZ1*       GAM( 9, 6)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 8)
     & +             CXMBX1*              CZMAZ1*       GAM( 9, 8)
     & +                                  CZMAZ1*CZMBZ1*GAM( 9,10)
     & +                                  CZMAZ1*       GAM(10, 2)
     & +2.0D0*       CXMBX1*              CZMAZ1*       GAM(10, 5)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 8)
     & +                           CYMBY1*CZMAZ1*       GAM(10, 8)
     & +                                  CZMAZ1*CZMBZ1*GAM(10, 9)
  
C      POLY(10)=      XXMCX2               
       COEF(10)=0.0D0
     & +                                                GAM( 1, 5)
     & +                                                GAM( 2, 2)
     & +      CXMAX1*                                   GAM( 2, 5)
     & +2.0D0*       CXMBX1*                            GAM( 2, 5)
     & +                           CYMBY1*              GAM( 2, 8)
     & +                                         CZMBZ1*GAM( 2, 9)
     & +                    CYMAY1*                     GAM( 3, 5)
     & +                                  CZMAZ1*       GAM( 4, 5)
     & +                                                GAM( 5, 1)
     & +2.0D0*CXMAX1*                                   GAM( 5, 2)
     & +             CXMBX1*                            GAM( 5, 2)
     & +                           CYMBY1*              GAM( 5, 3)
     & +                                         CZMBZ1*GAM( 5, 4)
     & +      CXMAX2*                                   GAM( 5, 5)
     & +4.0D0*CXMAX1*CXMBX1*                            GAM( 5, 5)
     & +             CXMBX2*                            GAM( 5, 5)
     & +                           CYMBY2*              GAM( 5, 6)
     & +                                         CZMBZ2*GAM( 5, 7)
     & +2.0D0*CXMAX1*              CYMBY1*              GAM( 5, 8)
       COEF(10)=COEF(10)
     & +             CXMBX1*       CYMBY1*              GAM( 5, 8)
     & +2.0D0*CXMAX1*                            CZMBZ1*GAM( 5, 9)
     & +             CXMBX1*                     CZMBZ1*GAM( 5, 9)
     & +                           CYMBY1*       CZMBZ1*GAM( 5,10)
     & +                    CYMAY2*                     GAM( 6, 5)
     & +                                  CZMAZ2*       GAM( 7, 5)
     & +                    CYMAY1*                     GAM( 8, 2)
     & +      CXMAX1*       CYMAY1*                     GAM( 8, 5)
     & +2.0D0*       CXMBX1*CYMAY1*                     GAM( 8, 5)
     & +                    CYMAY1*CYMBY1*              GAM( 8, 8)
     & +                    CYMAY1*              CZMBZ1*GAM( 8, 9)
     & +                                  CZMAZ1*       GAM( 9, 2)
     & +      CXMAX1*                     CZMAZ1*       GAM( 9, 5)
     & +2.0D0*       CXMBX1*              CZMAZ1*       GAM( 9, 5)
     & +                           CYMBY1*CZMAZ1*       GAM( 9, 8)
     & +                                  CZMAZ1*CZMBZ1*GAM( 9, 9)
     & +                    CYMAY1*       CZMAZ1*       GAM(10, 5)
  
C      POLY(11)=1.0D0              *ZZMCZ3 
       COEF(11)=0.0D0
     & +                                                GAM( 4, 7)
     & +                                                GAM( 7, 4)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 7)
     & +2.0D0*                                   CZMBZ1*GAM( 7, 7)
     & +             CXMBX1*                            GAM( 7, 9)
     & +                           CYMBY1*              GAM( 7,10)
     & +      CXMAX1*                                   GAM( 9, 7)
     & +                    CYMAY1*                     GAM(10, 7)
  
C      POLY(12)=1.0D0       *YYMCY1*ZZMCZ2 
       COEF(12)=0.0D0
     & +                                                GAM( 3, 7)
     & +                                                GAM( 4,10)
     & +2.0D0*              CYMAY1*                     GAM( 6, 7)
     & +                                                GAM( 7, 3)
     & +2.0D0*                     CYMBY1*              GAM( 7, 6)
     & +             CXMBX1*                            GAM( 7, 8)
     & +2.0D0*                            CZMAZ1*       GAM( 7,10)
     & +                                         CZMBZ1*GAM( 7,10)
     & +      CXMAX1*                                   GAM( 8, 7)
     & +      CXMAX1*                                   GAM( 9,10)
     & +                                                GAM(10, 4)
     & +                                  CZMAZ1*       GAM(10, 7)
     & +2.0D0*                                   CZMBZ1*GAM(10, 7)
     & +             CXMBX1*                            GAM(10, 9)
     & +                    CYMAY1*                     GAM(10,10)
     & +                           CYMBY1*              GAM(10,10)
  
C      POLY(13)=1.0D0       *YYMCY2*ZZMCZ1 
       COEF(13)=0.0D0
     & +                                                GAM( 3,10)
     & +                                                GAM( 4, 6)
     & +                                                GAM( 6, 4)
     & +2.0D0*                                   CZMBZ1*GAM( 6, 7)
     & +             CXMBX1*                            GAM( 6, 9)
     & +2.0D0*              CYMAY1*                     GAM( 6,10)
     & +                           CYMBY1*              GAM( 6,10)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 6)
     & +      CXMAX1*                                   GAM( 8,10)
     & +      CXMAX1*                                   GAM( 9, 6)
     & +                                                GAM(10, 3)
     & +                    CYMAY1*                     GAM(10, 6)
     & +2.0D0*                     CYMBY1*              GAM(10, 6)
     & +             CXMBX1*                            GAM(10, 8)
     & +                                  CZMAZ1*       GAM(10,10)
     & +                                         CZMBZ1*GAM(10,10)
  
C      POLY(14)=1.0D0       *YYMCY3        
       COEF(14)=0.0D0
     & +                                                GAM( 3, 6)
     & +                                                GAM( 6, 3)
     & +2.0D0*              CYMAY1*                     GAM( 6, 6)
     & +2.0D0*                     CYMBY1*              GAM( 6, 6)
     & +             CXMBX1*                            GAM( 6, 8)
     & +                                         CZMBZ1*GAM( 6,10)
     & +      CXMAX1*                                   GAM( 8, 6)
     & +                                  CZMAZ1*       GAM(10, 6)
  
C      POLY(15)=      XXMCX1       *ZZMCZ2 
       COEF(15)=0.0D0
     & +                                                GAM( 2, 7)
     & +                                                GAM( 4, 9)
     & +2.0D0*CXMAX1*                                   GAM( 5, 7)
     & +                                                GAM( 7, 2)
     & +2.0D0*       CXMBX1*                            GAM( 7, 5)
     & +                           CYMBY1*              GAM( 7, 8)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 9)
     & +                                         CZMBZ1*GAM( 7, 9)
     & +                    CYMAY1*                     GAM( 8, 7)
     & +                                                GAM( 9, 4)
     & +                                  CZMAZ1*       GAM( 9, 7)
     & +2.0D0*                                   CZMBZ1*GAM( 9, 7)
     & +      CXMAX1*                                   GAM( 9, 9)
     & +             CXMBX1*                            GAM( 9, 9)
     & +                           CYMBY1*              GAM( 9,10)
     & +                    CYMAY1*                     GAM(10, 9)
  
C      POLY(16)=      XXMCX1*YYMCY1*ZZMCZ1 
       COEF(16)=0.0D0
     & +                                                GAM( 2,10)
     & +                                                GAM( 3, 9)
     & +                                                GAM( 4, 8)
     & +2.0D0*CXMAX1*                                   GAM( 5,10)
     & +2.0D0*              CYMAY1*                     GAM( 6, 9)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 8)
     & +                                                GAM( 8, 4)
     & +2.0D0*                                   CZMBZ1*GAM( 8, 7)
     & +      CXMAX1*                                   GAM( 8, 9)
     & +             CXMBX1*                            GAM( 8, 9)
     & +                    CYMAY1*                     GAM( 8,10)
     & +                           CYMBY1*              GAM( 8,10)
     & +                                                GAM( 9, 3)
     & +2.0D0*                     CYMBY1*              GAM( 9, 6)
     & +      CXMAX1*                                   GAM( 9, 8)
     & +             CXMBX1*                            GAM( 9, 8)
     & +                                  CZMAZ1*       GAM( 9,10)
     & +                                         CZMBZ1*GAM( 9,10)
     & +                                                GAM(10, 2)
       COEF(16)=COEF(16)
     & +2.0D0*       CXMBX1*                            GAM(10, 5)
     & +                    CYMAY1*                     GAM(10, 8)
     & +                           CYMBY1*              GAM(10, 8)
     & +                                  CZMAZ1*       GAM(10, 9)
     & +                                         CZMBZ1*GAM(10, 9)
  
C      POLY(17)=      XXMCX1*YYMCY2        
       COEF(17)=0.0D0
     & +                                                GAM( 2, 6)
     & +                                                GAM( 3, 8)
     & +2.0D0*CXMAX1*                                   GAM( 5, 6)
     & +                                                GAM( 6, 2)
     & +2.0D0*       CXMBX1*                            GAM( 6, 5)
     & +2.0D0*              CYMAY1*                     GAM( 6, 8)
     & +                           CYMBY1*              GAM( 6, 8)
     & +                                         CZMBZ1*GAM( 6, 9)
     & +                                                GAM( 8, 3)
     & +                    CYMAY1*                     GAM( 8, 6)
     & +2.0D0*                     CYMBY1*              GAM( 8, 6)
     & +      CXMAX1*                                   GAM( 8, 8)
     & +             CXMBX1*                            GAM( 8, 8)
     & +                                         CZMBZ1*GAM( 8,10)
     & +                                  CZMAZ1*       GAM( 9, 6)
     & +                                  CZMAZ1*       GAM(10, 8)
  
C      POLY(18)=      XXMCX2       *ZZMCZ1 
       COEF(18)=0.0D0
     & +                                                GAM( 2, 9)
     & +                                                GAM( 4, 5)
     & +                                                GAM( 5, 4)
     & +2.0D0*                                   CZMBZ1*GAM( 5, 7)
     & +2.0D0*CXMAX1*                                   GAM( 5, 9)
     & +             CXMBX1*                            GAM( 5, 9)
     & +                           CYMBY1*              GAM( 5,10)
     & +2.0D0*                            CZMAZ1*       GAM( 7, 5)
     & +                    CYMAY1*                     GAM( 8, 9)
     & +                                                GAM( 9, 2)
     & +      CXMAX1*                                   GAM( 9, 5)
     & +2.0D0*       CXMBX1*                            GAM( 9, 5)
     & +                           CYMBY1*              GAM( 9, 8)
     & +                                  CZMAZ1*       GAM( 9, 9)
     & +                                         CZMBZ1*GAM( 9, 9)
     & +                    CYMAY1*                     GAM(10, 5)
  
C      POLY(19)=      XXMCX2*YYMCY1        
       COEF(19)=0.0D0
     & +                                                GAM( 2, 8)
     & +                                                GAM( 3, 5)
     & +                                                GAM( 5, 3)
     & +2.0D0*                     CYMBY1*              GAM( 5, 6)
     & +2.0D0*CXMAX1*                                   GAM( 5, 8)
     & +             CXMBX1*                            GAM( 5, 8)
     & +                                         CZMBZ1*GAM( 5,10)
     & +2.0D0*              CYMAY1*                     GAM( 6, 5)
     & +                                                GAM( 8, 2)
     & +      CXMAX1*                                   GAM( 8, 5)
     & +2.0D0*       CXMBX1*                            GAM( 8, 5)
     & +                    CYMAY1*                     GAM( 8, 8)
     & +                           CYMBY1*              GAM( 8, 8)
     & +                                         CZMBZ1*GAM( 8, 9)
     & +                                  CZMAZ1*       GAM( 9, 8)
     & +                                  CZMAZ1*       GAM(10, 5)
  
C      POLY(20)=      XXMCX3               
       COEF(20)=0.0D0
     & +                                                GAM( 2, 5)
     & +                                                GAM( 5, 2)
     & +2.0D0*CXMAX1*                                   GAM( 5, 5)
     & +2.0D0*       CXMBX1*                            GAM( 5, 5)
     & +                           CYMBY1*              GAM( 5, 8)
     & +                                         CZMBZ1*GAM( 5, 9)
     & +                    CYMAY1*                     GAM( 8, 5)
     & +                                  CZMAZ1*       GAM( 9, 5)
  
C      POLY(21)=1.0D0              *ZZMCZ4 
       COEF(21)=0.0D0
     & +                                                GAM( 7, 7)
  
C      POLY(22)=1.0D0       *YYMCY1*ZZMCZ3 
       COEF(22)=0.0D0
     & +                                                GAM( 7,10)
     & +                                                GAM(10, 7)
  
C      POLY(23)=1.0D0       *YYMCY2*ZZMCZ2 
       COEF(23)=0.0D0
     & +                                                GAM( 6, 7)
     & +                                                GAM( 7, 6)
     & +                                                GAM(10,10)
  
C      POLY(24)=1.0D0       *YYMCY3*ZZMCZ1 
       COEF(24)=0.0D0
     & +                                                GAM( 6,10)
     & +                                                GAM(10, 6)
  
C      POLY(25)=1.0D0       *YYMCY4        
       COEF(25)=0.0D0
     & +                                                GAM( 6, 6)
  
C      POLY(26)=      XXMCX1       *ZZMCZ3 
       COEF(26)=0.0D0
     & +                                                GAM( 7, 9)
     & +                                                GAM( 9, 7)
  
C      POLY(27)=      XXMCX1*YYMCY1*ZZMCZ2 
       COEF(27)=0.0D0
     & +                                                GAM( 7, 8)
     & +                                                GAM( 8, 7)
     & +                                                GAM( 9,10)
     & +                                                GAM(10, 9)
  
C      POLY(28)=      XXMCX1*YYMCY2*ZZMCZ1 
       COEF(28)=0.0D0
     & +                                                GAM( 6, 9)
     & +                                                GAM( 8,10)
     & +                                                GAM( 9, 6)
     & +                                                GAM(10, 8)
  
C      POLY(29)=      XXMCX1*YYMCY3        
       COEF(29)=0.0D0
     & +                                                GAM( 6, 8)
     & +                                                GAM( 8, 6)
  
C      POLY(30)=      XXMCX2       *ZZMCZ2 
       COEF(30)=0.0D0
     & +                                                GAM( 5, 7)
     & +                                                GAM( 7, 5)
     & +                                                GAM( 9, 9)
  
C      POLY(31)=      XXMCX2*YYMCY1*ZZMCZ1 
       COEF(31)=0.0D0
     & +                                                GAM( 5,10)
     & +                                                GAM( 8, 9)
     & +                                                GAM( 9, 8)
     & +                                                GAM(10, 5)
  
C      POLY(32)=      XXMCX2*YYMCY2        
       COEF(32)=0.0D0
     & +                                                GAM( 5, 6)
     & +                                                GAM( 6, 5)
     & +                                                GAM( 8, 8)
  
C      POLY(33)=      XXMCX3       *ZZMCZ1 
       COEF(33)=0.0D0
     & +                                                GAM( 5, 9)
     & +                                                GAM( 9, 5)
  
C      POLY(34)=      XXMCX3*YYMCY1        
       COEF(34)=0.0D0
     & +                                                GAM( 5, 8)
     & +                                                GAM( 8, 5)
  
C      POLY(35)=      XXMCX4               
       COEF(35)=0.0D0
     & +                                                GAM( 5, 5)
           ARG=EXP(-ARG)
           DO I=1,35
           COEF(I)=COEF(I)*ARG
           END DO
C POLY(1 ) 1.0D0
C POLY(2 ) Z
C POLY(3 ) Y
C POLY(4 ) X
C POLY(5 ) ZZ
C POLY(6 ) YZ
C POLY(7 ) YY
C POLY(8 ) XZ
C POLY(9 ) XY
C POLY(10) XX
C POLY(11) ZZZ
C POLY(12) YZZ
C POLY(13) YYZ
C POLY(14) YYY
C POLY(15) XZZ
C POLY(16) XYZ
C POLY(17) XYY  
C POLY(18) XXZ 
C POLY(19) XXY 
C POLY(20) XXX 
C POLY(21) ZZZZ
C POLY(22) YZZZ
C POLY(23) YYZZ
C POLY(24) YYYZ
C POLY(25) YYYY
C POLY(26) XZZZ
C POLY(27) XYZZ
C POLY(28) XYYZ
C POLY(29) XYYY
C POLY(30) XXZZ
C POLY(31) XXYZ
C POLY(32) XXYY
C POLY(33) XXXZ
C POLY(34) XXXY
C POLY(35) XXXX
C DENSITY= COEF(I)*POLY(I)*EXP(-DLT*(R-C)**2) SUMMED OVER I
       IF(GGA)THEN 
             DO I=2,10
             DO J=1,20
             DOEF(J,I)=0.0D0
             END DO
             END DO
       DOEF( 1, 2)=       COEF( 4)!x  X                 
       DOEF( 2, 2)=       COEF( 8)!x  XZ        Z       
       DOEF( 3, 2)=       COEF( 9)!x  XY        Y       
       DOEF( 4, 2)= 2.0D0*COEF(10)!x  XX        X       
       DOEF( 5, 2)=       COEF(15)!x  XZZ       ZZ      
       DOEF( 6, 2)=       COEF(16)!x  XYZ       YZ      
       DOEF( 7, 2)=       COEF(17)!x  XYY       YY      
       DOEF( 8, 2)= 2.0D0*COEF(18)!x  XXZ       XZ      
       DOEF( 9, 2)= 2.0D0*COEF(19)!x  XXY       XY      
       DOEF(10, 2)= 3.0D0*COEF(20)!x  XXX       XX      
       DOEF(11, 2)=       COEF(26)!x  XZZZ      ZZZ     
       DOEF(12, 2)=       COEF(27)!x  XYZZ      YZZ     
       DOEF(13, 2)=       COEF(28)!x  XYYZ      YYZ     
       DOEF(14, 2)=       COEF(29)!x  XYYY      YYY     
       DOEF(15, 2)= 2.0D0*COEF(30)!x  XXZZ      XZZ     
       DOEF(16, 2)= 2.0D0*COEF(31)!x  XXYZ      XYZ     
       DOEF(17, 2)= 2.0D0*COEF(32)!x  XXYY      XYY     
       DOEF(18, 2)= 3.0D0*COEF(33)!x  XXXZ      XXZ     
       DOEF(19, 2)= 3.0D0*COEF(34)!x  XXXY      XXY     
       DOEF(20, 2)= 4.0D0*COEF(35)!x  XXXX      XXX     
       DOEF( 1, 3)=       COEF( 3)!y  Y                 
       DOEF( 2, 3)=       COEF( 6)!y  YZ        Z       
       DOEF( 3, 3)= 2.0D0*COEF( 7)!y  YY        Y       
       DOEF( 4, 3)=       COEF( 9)!y  XY        X       
       DOEF( 5, 3)=       COEF(12)!y  YZZ       ZZ      
       DOEF( 6, 3)= 2.0D0*COEF(13)!y  YYZ       YZ      
       DOEF( 7, 3)= 3.0D0*COEF(14)!y  YYY       YY      
       DOEF( 8, 3)=       COEF(16)!y  XYZ       XZ      
       DOEF( 9, 3)= 2.0D0*COEF(17)!y  XYY       XY      
       DOEF(10, 3)=       COEF(19)!y  XXY       XX      
       DOEF(11, 3)=       COEF(22)!y  YZZZ      ZZZ     
       DOEF(12, 3)= 2.0D0*COEF(23)!y  YYZZ      YZZ     
       DOEF(13, 3)= 3.0D0*COEF(24)!y  YYYZ      YYZ     
       DOEF(14, 3)= 4.0D0*COEF(25)!y  YYYY      YYY     
       DOEF(15, 3)=       COEF(27)!y  XYZZ      XZZ     
       DOEF(16, 3)= 2.0D0*COEF(28)!y  XYYZ      XYZ     
       DOEF(17, 3)= 3.0D0*COEF(29)!y  XYYY      XYY     
       DOEF(18, 3)=       COEF(31)!y  XXYZ      XXZ     
       DOEF(19, 3)= 2.0D0*COEF(32)!y  XXYY      XXY     
       DOEF(20, 3)=       COEF(34)!y  XXXY      XXX     
       DOEF( 1, 4)=       COEF( 2)!z  Z                 
       DOEF( 2, 4)= 2.0D0*COEF( 5)!z  ZZ        Z       
       DOEF( 3, 4)=       COEF( 6)!z  YZ        Y       
       DOEF( 4, 4)=       COEF( 8)!z  XZ        X       
       DOEF( 5, 4)= 3.0D0*COEF(11)!z  ZZZ       ZZ      
       DOEF( 6, 4)= 2.0D0*COEF(12)!z  YZZ       YZ      
       DOEF( 7, 4)=       COEF(13)!z  YYZ       YY      
       DOEF( 8, 4)= 2.0D0*COEF(15)!z  XZZ       XZ      
       DOEF( 9, 4)=       COEF(16)!z  XYZ       XY      
       DOEF(10, 4)=       COEF(18)!z  XXZ       XX      
       DOEF(11, 4)= 4.0D0*COEF(21)!z  ZZZZ      ZZZ     
       DOEF(12, 4)= 3.0D0*COEF(22)!z  YZZZ      YZZ     
       DOEF(13, 4)= 2.0D0*COEF(23)!z  YYZZ      YYZ     
       DOEF(14, 4)=       COEF(24)!z  YYYZ      YYY     
       DOEF(15, 4)= 3.0D0*COEF(26)!z  XZZZ      XZZ     
       DOEF(16, 4)= 2.0D0*COEF(27)!z  XYZZ      XYZ     
       DOEF(17, 4)=       COEF(28)!z  XYYZ      XYY     
       DOEF(18, 4)= 2.0D0*COEF(30)!z  XXZZ      XXZ     
       DOEF(19, 4)=       COEF(31)!z  XXYZ      XXY     
       DOEF(20, 4)=       COEF(33)!z  XXXZ      XXX     
       DOEF( 1, 5)= 2.0D0*COEF(10)!xx  XX                
       DOEF( 2, 5)= 2.0D0*COEF(18)!xx  XXZ       Z       
       DOEF( 3, 5)= 2.0D0*COEF(19)!xx  XXY       Y       
       DOEF( 4, 5)= 6.0D0*COEF(20)!xx  XXX       X       
       DOEF( 5, 5)= 2.0D0*COEF(30)!xx  XXZZ      ZZ      
       DOEF( 6, 5)= 2.0D0*COEF(31)!xx  XXYZ      YZ      
       DOEF( 7, 5)= 2.0D0*COEF(32)!xx  XXYY      YY      
       DOEF( 8, 5)= 6.0D0*COEF(33)!xx  XXXZ      XZ      
       DOEF( 9, 5)= 6.0D0*COEF(34)!xx  XXXY      XY      
       DOEF(10, 5)=12.0D0*COEF(35)!xx  XXXX      XX      
       DOEF( 1, 6)= 2.0D0*COEF( 7)!yy  YY                
       DOEF( 2, 6)= 2.0D0*COEF(13)!yy  YYZ       Z       
       DOEF( 3, 6)= 6.0D0*COEF(14)!yy  YYY       Y       
       DOEF( 4, 6)= 2.0D0*COEF(17)!yy  XYY       X       
       DOEF( 5, 6)= 2.0D0*COEF(23)!yy  YYZZ      ZZ      
       DOEF( 6, 6)= 6.0D0*COEF(24)!yy  YYYZ      YZ      
       DOEF( 7, 6)=12.0D0*COEF(25)!yy  YYYY      YY      
       DOEF( 8, 6)= 2.0D0*COEF(28)!yy  XYYZ      XZ      
       DOEF( 9, 6)= 6.0D0*COEF(29)!yy  XYYY      XY      
       DOEF(10, 6)= 2.0D0*COEF(32)!yy  XXYY      XX      
       DOEF( 1, 7)= 2.0D0*COEF( 5)!zz  ZZ                
       DOEF( 2, 7)= 6.0D0*COEF(11)!zz  ZZZ       Z       
       DOEF( 3, 7)= 2.0D0*COEF(12)!zz  YZZ       Y       
       DOEF( 4, 7)= 2.0D0*COEF(15)!zz  XZZ       X       
       DOEF( 5, 7)=12.0D0*COEF(21)!zz  ZZZZ      ZZ      
       DOEF( 6, 7)= 6.0D0*COEF(22)!zz  YZZZ      YZ      
       DOEF( 7, 7)= 2.0D0*COEF(23)!zz  YYZZ      YY      
       DOEF( 8, 7)= 6.0D0*COEF(26)!zz  XZZZ      XZ      
       DOEF( 9, 7)= 2.0D0*COEF(27)!zz  XYZZ      XY      
       DOEF(10, 7)= 2.0D0*COEF(30)!zz  XXZZ      XX      
       DOEF( 1, 8)=       COEF( 9)!xy  XY                
       DOEF( 2, 8)=       COEF(16)!xy  XYZ       Z       
       DOEF( 3, 8)= 2.0D0*COEF(17)!xy  XYY       Y       
       DOEF( 4, 8)= 2.0D0*COEF(19)!xy  XXY       X       
       DOEF( 5, 8)=       COEF(27)!xy  XYZZ      ZZ      
       DOEF( 6, 8)= 2.0D0*COEF(28)!xy  XYYZ      YZ      
       DOEF( 7, 8)= 3.0D0*COEF(29)!xy  XYYY      YY      
       DOEF( 8, 8)= 2.0D0*COEF(31)!xy  XXYZ      XZ      
       DOEF( 9, 8)= 4.0D0*COEF(32)!xy  XXYY      XY      
       DOEF(10, 8)= 3.0D0*COEF(34)!xy  XXXY      XX      
       DOEF( 1, 9)= 1.0D0*COEF( 8)!xz  XZ                
       DOEF( 2, 9)= 2.0D0*COEF(15)!xz  XZZ       Z       
       DOEF( 3, 9)= 1.0D0*COEF(16)!xz  XYZ       Y       
       DOEF( 4, 9)= 2.0D0*COEF(18)!xz  XXZ       X       
       DOEF( 5, 9)= 3.0D0*COEF(26)!xz  XZZZ      ZZ      
       DOEF( 6, 9)= 2.0D0*COEF(27)!xz  XYZZ      YZ      
       DOEF( 7, 9)= 1.0D0*COEF(28)!xz  XYYZ      YY      
       DOEF( 8, 9)= 4.0D0*COEF(30)!xz  XXZZ      XZ      
       DOEF( 9, 9)= 2.0D0*COEF(31)!xz  XXYZ      XY      
       DOEF(10, 9)= 3.0D0*COEF(33)!xz  XXXZ      XX      
       DOEF( 1,10)=       COEF( 6)!yz  YZ                
       DOEF( 2,10)= 2.0D0*COEF(12)!yz  YZZ       Z       
       DOEF( 3,10)= 2.0D0*COEF(13)!yz  YYZ       Y       
       DOEF( 4,10)=       COEF(16)!yz  XYZ       X       
       DOEF( 5,10)= 3.0D0*COEF(22)!yz  YZZZ      ZZ      
       DOEF( 6,10)= 4.0D0*COEF(23)!yz  YYZZ      YZ      
       DOEF( 7,10)= 3.0D0*COEF(24)!yz  YYYZ      YY      
       DOEF( 8,10)= 2.0D0*COEF(27)!yz  XYZZ      XZ      
       DOEF( 9,10)= 2.0D0*COEF(28)!yz  XYYZ      XY      
       DOEF(10,10)=       COEF(31)!yz  XXYZ      XX      
       END IF
       RETURN
       END
