C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE THCNOV(A1,A2,A3,B,C,OV)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL NOD
      DIMENSION AE(3),BE(3),OV(ND,ND),B(3),C(3)
      SAVE
      DATA NOD/.FALSE./
      DATA PI/3.141592653589793D0/
      DO 5 K=1,ND
      DO 5 J=1,ND
    5 OV(J,K)=0.0D0
      AT=A1+A2+A3
      EX=(A2*B(1)+A3*C(1))/AT
      EY=(A2*B(2)+A3*C(2))/AT
      EZ=(A2*B(3)+A3*C(3))/AT
      E2=EX**2+EY**2+EZ**2
      B2=B(1)**2+B(2)**2+B(3)**2
      C2=C(1)**2+C(2)**2+C(3)**2
      CCC=AT*E2-A2*B2-A3*C2
      COEF=SQRT(PI/AT)**3*EXP(CCC)
      AE(1)=EX
      AE(2)=EY
      AE(3)=EZ
      BE(1)=(EX-B(1))
      BE(2)=(EY-B(2))
      BE(3)=(EZ-B(3))
C THIS SS
      OV(1,1)=COEF
      DO 105 I=1,3
C THIS IS SPX
      OV(1,I+1)=COEF*BE(I)
C THIS IS PXS
      OV(I+1,1)=COEF*AE(I)
      IF(NOD) GO TO 1000
C THIS IS SDXX
      OV(1,I+4)=COEF*(BE(I)**2+0.5D0/AT)
C THIS IS DXXS
      OV(I+4,1)=COEF*(AE(I)**2+0.5D0/AT)
 1000 CONTINUE
C THIS IS PXPX
      OV(I+1,I+1)=COEF*(AE(I)*BE(I)+0.5D0/AT)
      IF(NOD) GO TO 1001
C THIS IS DXXPX
      OV(I+4,I+1)=COEF*(BE(I)*AE(I)**2+0.5D0*(2.0D0*AE(I)+BE(I))/AT)
C THIS IS PXDXX
      OV(I+1,I+4)=COEF*(AE(I)*BE(I)**2+0.5D0*(2.0D0*BE(I)+AE(I))/AT)
C THIS IS DXXDXX
      OV(I+4,I+4)=COEF*(AE(I)**2*BE(I)**2+0.5D0*(AE(I)**2+BE(I)**2) 
     &           /AT+2*BE(I)*AE(I)/AT+0.75D0/AT**2)
 1001 CONTINUE
  105 CONTINUE
      IF(NOD) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      OV(1,K)=COEF*BE(I)*BE(J)
C THIS IS DXYS
      OV(K,1)=COEF*AE(I)*AE(J)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      OV(I+1,J+1)=COEF*AE(I)*BE(J)
      IF(NOD) GO TO 1003
C THIS IS PXDYY
      OV(I+1,J+4)=COEF*AE(I)*(BE(J)**2+0.5D0/AT)
C THIS IS DXXPY
      OV(I+4,J+1)=COEF*BE(J)*(AE(I)**2+0.5D0/AT)
C THIS IS DXXDYY
      OV(I+4,J+4)=COEF*(AE(I)**2+0.5D0/AT)*(BE(J)**2+0.5D0/AT)
 1003 CONTINUE
  115 CONTINUE
      IF(NOD) GO TO 1004
      DO 130 I=1,3
      DO 130 J=1,2
      DO 130 K=2,3
      IF (J.EQ.K) GO TO 130
      IF (I.NE.J.AND.I.NE.K) GO TO 120
      IF (J.EQ.I) GO TO 118
      L=J
      GO TO 119
  118 L=K
  119 CONTINUE
C THIS IS PXDXY
      OV(I+1,J+K+5)=COEF*BE(L)*(AE(I)*BE(I)+0.5D0/AT)
C THIS IS DXXDXY
      OV(I+4,J+K+5)=COEF*BE(L)*(AE(I)**2*BE(I)
     &             +0.5D0*(2*AE(I)+BE(I))/AT)
C THIS IS DXYPX
      OV(J+K+5,I+1)=COEF*AE(L)*(AE(I)*BE(I)+0.5D0/AT)
C THIS IS DXYDXX
      OV(J+K+5,I+4)=COEF*AE(L)*(AE(I)*BE(I)**2
     &             +0.5D0*(2*BE(I)+AE(I))/AT)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      OV(I+1,J+K+5)=COEF*BE(J)*BE(K)*AE(I)
C THIS IS DXXDYZ
      OV(I+4,J+K+5)=COEF*BE(J)*BE(K)*(AE(I)**2+0.5D0/AT)
C THIS IS DYZPX
      OV(J+K+5,I+1)=COEF*AE(J)*AE(K)*BE(I)
C THIS IS DYZDXX
      OV(J+K+5,I+4)=COEF*AE(J)*AE(K)*(BE(I)**2+0.5D0/AT)
  130 CONTINUE
      DO 146 I=1,2
      DO 146 J=2,3
      IF (I.EQ.J) GO TO 146
      DO 145 K=1,2
      DO 145 L=2,3
      IF (K.EQ.L) GO TO 145
      IF (I.EQ.K.AND.J.EQ.L) GO TO 140
      IF (I.EQ.K) GO TO 131
      IF (I.EQ.L) GO TO 132
      IF (J.EQ.K) GO TO 133
      IF  (J.EQ.L) GO TO 134
  131 M1=I
      M2=J
      M3=L
      GO TO 135
  132 M1=I
      M2=J
      M3=K
      GO TO 135
  133 M1=J
      M2=I
      M3=L
      GO TO 135
  134 M1=J
      M2=I
      M3=K
      GO TO 135
  135 CONTINUE
C THIS IS DXYDYZ
      OV(I+J+5,K+L+5)=COEF*AE(M2)*BE(M3)*(AE(M1)*BE(M1)+0.5D0/AT)
      GO TO 145
  140 CONTINUE
C THIS IS DXYDXY
      OV(I+J+5,K+L+5)=COEF*(AE(I)*BE(I)+0.5D0/AT)
     &               *(AE(J)*BE(J)+0.5D0/AT)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      RETURN
      END
