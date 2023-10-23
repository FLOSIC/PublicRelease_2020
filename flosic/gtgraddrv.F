C UTEP Electronic Structure Lab (2020)
C YY. Jan 2016
C This subroutine computes d/dR_v d/dx_i phi where R_v can be
C Rx, Ry, or Rz and x_i can be x, y, z for a given point of
C X, Y, Z. X = x - Rx, Y = y - Ry, and Z = z - Rz respectively.
C
C This program is modified version of GTGRAD subroutine which computes
C dPhi/dR_v at a point X, Y, and Z where X = x - Rx and so on.
C
C INPUTS: X,Y,Z: x - Rx, y - Ry, and z - Rz
C         I_ATOM
C OUTPUT: GRAD: d/dRx, d/dRy, and d/dRz of basis functions  
C         GGRAD: d/dx_i d/dRv of basis functions
C

C PPANG(1,:): d/dx d/dRx ANG = -d/dx d/dx ANG = -d/dRx d/dRx ANG
C PPANG(2,:): d/dy d/dRy ANG = -d/dy d/dy ANG = -d/dRy d/dRy ANG
C PPANG(3,:): d/dz d/dRz ANG = -d/dz d/dz ANG = -d/dRz d/dRz ANG
C PPANG(4,:): d/dx d/dRy ANG = -d/d
C PPANG(5,:): d/dx d/dRz ANG
C PPANG(6,:): d/dy d/dRz ANG 
C duplicate/double check
C PPANG(7,:): d/dy d/dRx
C PPANG(8,:): d/dz d/dRx
C PPANG(9,:): d/dz d/dRy

C GGRAD(1,N_CON,MU,LI): d/dx d/dRx phi
C GGRAD(2,N_CON,MU,LI): d/dy d/dRy phi
C GGRAD(3,N_CON,MU,LI): d/dz d/dRz phi
C GGRAD(4,N_CON,MU,LI): d/dx d/dRy phi
C GGRAD(5,N_CON,MU,LI): d/dx d/dRz phi
C GGRAD(6,N_CON,MU,LI): d/dy d/dRz phi


C *********************************************************************
C
C GTGRAD CALCULATES THE NEGATIVE GRADIENT OF THE ATOMIC BASIS 
C FUNCTIONS OF ONE ATOM FOR ONE POINT
C
       SUBROUTINE GTGRADDRV(X,Y,Z,I_ATOM,GRAD,GGRAD)
C ORIGINAL VERSION BY KOBLAR A. JACKSON (1988)
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
       INCLUDE 'PARAMA2'
       SAVE
       INTEGER :: IBEG, IEND, LMAX1, I_ATOM, I, J, LI, IANG, 
     &            NB, IX, I_BARE, MU 
       REAL*8  :: GEPS, ANG, X, Y, Z, RSQR, PANG, PANGM, PPANG, 
     &            GRAD, GGRAD, ALF, REXPON, ALF2, ALF2M, GAUS, 
     &            EXP, DANG, DANGIT

       DIMENSION GRAD(3,MAX_CON,10,3),    !GRAD(d/Rv,N_CON,MU,LI)
     &           GGRAD(6,MAX_CON,10,3),   !GGRAD(d/dx_i,d/dR_v,N_CON,Mu,LI)
     &           ANG(10),        ! ANG(#fn): Angular function/Polynomials: (x-Rx)^l + (y-Ry)^m + (z-Rz)^n
     &           PANG(3,10),     ! PANG(RxRyRz, #fn): derivatives of ANG with respect to Rx, Ry, Rz
     &           PANGM(3,10),    ! derivatives of ANG with respect to x, y, z
     &           PPANG(6,10),    ! PPANG(iv, #fn): drivatives of PANG with respect to x, y, z. d^2/dxdRx
     &           DANG(3,10),
     &           DANGIT(6,10),
     &           ALF2(3),        ! Derivatives of alpha*(X^2+Y^2+Z^2) with respect to Rx, Ry, Rz: 2*alpha*(x-Rx).
     &           ALF2M(3)        ! Derivatives of alpha*(X^2+Y^2+Z^2) with respect to x, y, z: -2*alpha*(x-Rx)
       DIMENSION IBEG(3), IEND(3)
C
       DATA IBEG /1,2,5/
       DATA IEND /1,4,10/
C
       LMAX1 = LSYMMAX(I_ATOM)+1
       GEPS = EXP(-CUTEXP)
C
C GET ANGULAR FUNCTION ARRAY
C
       ANG(1) = 1.0D0      !1.0
       ANG(2) = X          !x-Rx
       ANG(3) = Y          !y-Ry
       ANG(4) = Z          !z-Rz
       ANG(5) = X*X        !(x-Rx)**2
       ANG(6) = Y*Y        !(y-Ry)**2
       ANG(7) = Z*Z        !(z-Rz)**2
       ANG(8) = X*Y        !(x-Rx)(y-Ry)
       ANG(9) = X*Z        !(x-Rx)(z-Rz)
       ANG(10)= Y*Z        !(y-Ry)(z-Rz)
       RSQR=ANG(5)+ANG(6)+ANG(7)
C
C GET DERIVATIVE ARRAY FOR POLYNOMIAL FUNCTIONS:  -D/DX(ANG)
C
       DO 3 I=1,3
        DO 2 J=1,10
         PANG(I,J)=0.0D0
         PANGM(I,J)=0.0d0
         PPANG(2*I-1,J)=0.0d0
         PPANG(2*I,J)=0.0d0
    2   CONTINUE
    3  CONTINUE
C
C PANG's are d/dRv of the polynomial ANG. !d/dx d/dRv in comments
C (PPANG(d/dx,d/dRv,MU)) where 1,2,and 3 are derivatives wrt x,y,and z.
C X-DERIVATIVE
C
       PANG(1,2) = -1.0D0
       PANG(1,5) = -2*X
       PANG(1,8) = -Y
       PANG(1,9) = -Z

       PANGM(1,2) = 1.0D0
       PANGM(1,5) = 2*X
       PANGM(1,8) = Y
       PANGM(1,9) = Z
C
C  Y DERIVATIVE
C
       PANG(2,3) = -1.0D0
       PANG(2,6) = -2*Y
       PANG(2,8) = -X
       PANG(2,10)= -Z

       PANGM(2,3) = 1.0D0
       PANGM(2,6) = 2*Y
       PANGM(2,8) = X
       PANGM(2,10)= Z
C
C  Z DERIVATIVE
C
       PANG(3,4) = -1.0D0
       PANG(3,7) = -2*Z
       PANG(3,9) = -X
       PANG(3,10)= -Y

       PANGM(3,4) = 1.0D0
       PANGM(3,7) = 2*Z
       PANGM(3,9) = X
       PANGM(3,10)= Y
C
C PPANG's are d/dx_i d/dRv of the polynomial ANG.
C d/dx terms:
       PPANG(1,5) = -2.0d0
       PPANG(4,8) = -1.0d0
       PPANG(5,9) = -1.0d0
C d/dy terms:
       PPANG(2,6) = -2.0d0
       PPANG(4,8) = -1.0d0
       PPANG(6,10)= -1.0d0
C d/dz terms:
       PPANG(3,7) = -2.0d0
       PPANG(5,9) = -1.0d0
       PPANG(6,10)= -1.0d0
C
C INITIALIZE GRAD
C
       DO 20 LI=1,LMAX1
        DO 15 IANG=1,IEND(LI)
         DO 10 NB=1,N_CON(LI,I_ATOM)
          DO 5 IX=1,3
           GRAD(IX,NB,IANG,LI)=0.0D0
           GGRAD(2*IX-1,NB,IANG,LI)=0.0d0
           GGRAD(2*IX,NB,IANG,LI)=0.0d0
    5     CONTINUE
   10    CONTINUE
   15   CONTINUE
   20  CONTINUE
C
C  GET GRADIENT TERM FOR DIFFERENT ANGULAR FUNCTIONS
C
       DO 70 I_BARE=1,N_BARE(I_ATOM)
        ALF=BFALP(I_BARE,I_ATOM)
        REXPON=ALF*RSQR
        IF (REXPON.GT.CUTEXP) GOTO 70
        ALF2(1)=2.0d0*ALF*X
        ALF2(2)=2.0d0*ALF*Y
        ALF2(3)=2.0d0*ALF*Z
        ALF2M(1)=-2.0d0*ALF*X
        ALF2M(2)=-2.0d0*ALF*Y
        ALF2M(3)=-2.0d0*ALF*Z
        GAUS=EXP(-REXPON)
C
C X,Y,Z DERIVATIVE
C
! DANG: d/dR[ Ang GAUS ]
        DO 30 I=1,3
         DO 25 J=1,IEND(LMAX1)
          DANG(I,J) = ALF2(I)*ANG(J)+PANG(I,J)
!Old version
!          DANGIT(I,J)=PPANG(I,J)-PANG(I,J)*2.0d0*ALF2(I)
!     &               +ANG(J)*2.0d0*ALF
!     &               -ANG(J)*ALF2(I)*ALF2(I)
!New version: equivalent of old version
          DANGIT(I,J) =PPANG(I,J)+PANGM(I,J)*ALF2(I)
     &                +PANG(I,J)*ALF2M(I)
     &                +ANG(J)*(2.0D0)*ALF
     &                +ANG(J)*ALF2(I)*ALF2M(I)
   25    CONTINUE
   30   CONTINUE

        do J=1,IEND(LMAX1)
! Old
!         DANGIT(4,J)=PPANG(4,J)-PANG(1,J)*ALF2(2)
!     &              -PANG(2,J)*ALF2(1)
!     &              -ANG(J)*ALF2(1)*ALF2(2)
!         DANGIT(5,J)=PPANG(5,J)-PANG(1,J)*ALF2(3)
!     &              -PANG(3,J)*ALF2(1)
!     &              -ANG(J)*ALF2(1)*ALF2(3)
!         DANGIT(6,J)=PPANG(6,J)-PANG(2,J)*ALF2(3)
!     &              -PANG(3,J)*ALF2(2)
!     &              -ANG(J)*ALF2(2)*ALF2(3)

! New: equivalent of old
         DANGIT(4,J)=PPANG(4,J)+PANGM(1,J)*ALF2(2)
     &              +PANG(2,J)*ALF2M(1)
     &              -ANG(J)*(0.0D0)*ALF
     &              +ANG(J)*ALF2(1)*ALF2M(2)
         DANGIT(5,J)=PPANG(5,J)+PANGM(1,J)*ALF2(3)
     &              +PANG(3,J)*ALF2M(1)
     &              -ANG(J)*(0.0D0)*ALF
     &              +ANG(J)*ALF2(1)*ALF2M(3)
         DANGIT(6,J)=PPANG(6,J)+PANGM(2,J)*ALF2(3)
     &              +PANG(3,J)*ALF2M(2)
     &              -ANG(J)*(0.0D0)*ALF
     &              +ANG(J)*ALF2(2)*ALF2M(3)
        end do
C
C CALCULATING DERIVATIVE
C
        DO 60 LI=1,LMAX1
         DO 55 NB=1,N_CON(LI,I_ATOM)
          IF (ABS(BFCON(I_BARE,NB,LI,I_ATOM)*GAUS).GT.GEPS) THEN
           DO 50 MU=IBEG(LI),IEND(LI)
            DO 45 IX=1,3
             GRAD(IX,NB,MU,LI)=GRAD(IX,NB,MU,LI) 
     &       +BFCON(I_BARE,NB,LI,I_ATOM)*GAUS*DANG(IX,MU)
 
             GGRAD(2*IX-1,NB,MU,LI)=GGRAD(2*IX-1,NB,MU,LI)
     &       +BFCON(I_BARE,NB,LI,I_ATOM)*GAUS*DANGIT(2*IX-1,MU)
             GGRAD(2*IX,NB,MU,LI)=GGRAD(2*IX,NB,MU,LI)
     &       +BFCON(I_BARE,NB,LI,I_ATOM)*GAUS*DANGIT(2*IX,MU)
   45       CONTINUE
   50      CONTINUE
          END IF
   55    CONTINUE
   60   CONTINUE
   70  CONTINUE
       RETURN
       END
