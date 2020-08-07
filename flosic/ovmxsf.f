C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE OVMXSF(A1,A2,R,T,WO)
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL DYES
      DIMENSION WO(ND,ND),R(3),T(3),X(3)
      SAVE
      DATA DYES/.FALSE./
      X(1)=T(1)-R(1)
      X(2)=T(2)-R(2)
      X(3)=T(3)-R(3)
      A=A1*A2/(A1+A2)
      AS=A*A
      PI=3.14159265358979324D0
      B=PI/(A1+A2)
      B=SQRT(B)
      B=B**3
      DO 100 I=1,ND
      DO 100 J=1,ND
  100 WO(I,J)=0.0D0
      XS=X(1)**2+X(2)**2+X(3)**2
      IF( XS.LT.1.0D-14 ) GO TO 150
      C=A*XS
      C=EXP(-C)
C THIS SS
      WO(1,1)=B*C
      DO 105 I=1,3
C THIS IS SPX
      WO(1,I+1)=-A*B*C*X(I)/A2
C THIS IS PXS
      WO(I+1,1)=A*B*C*X(I)/A1
      IF(DYES) GO TO 1000
C THIS IS SDXX
      WO(1,I+4)=(0.5D0*B*C/A2)*((1.0D0-A/A2)+(2*A*A*X(I)*X(I)/A2))
C THIS IS DXXS
      WO(I+4,1)=(0.5D0*B*C/A1)*((1.0D0-A/A1)+(2*A*A*X(I)*X(I)/A1))
 1000 CONTINUE
C THIS IS PXPX
      WO(I+1,I+1)=A*B*C*(0.5D0-A*X(I)*X(I))/(A1*A2)
      IF(DYES) GO TO 1001
C THIS IS DXXPX
      WO(I+4,I+1)=-A*B*C*X(I)*((1.0D0-3*A/A1)+(2*A*A*X(I)*X(I)/A1))
     &           /(2*A1*A2)
C THIS IS PXDXX
      WO(I+1,I+4)= A*B*C*X(I)*((1.0D0-3*A/A2)+(2*A*A*X(I)*X(I)/A2))
     &           /(2*A1*A2)
C THIS IS DXXDXX
      WOW=3.0D0+X(I)**2*(2*(A1+A2)-12*A)+4*AS*X(I)**4
      WO(I+4,I+4)=WOW*B*C*AS/(4*A1**2*A2**2)
 1001 CONTINUE
  105 CONTINUE
      IF(DYES) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      WO(1,K)=A*A*B*C*X(I)*X(J)/(A2*A2)
C THIS IS DXYS
      WO(K,1)=A*A*B*C*X(I)*X(J)/(A1*A1)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      WO(I+1,J+1)=-A*A*B*C*X(I)*X(J)/(A1*A2)
      IF(DYES) GO TO 1003
C THIS IS PXDYY
      WO(I+1,J+4)=A*B*C*X(I)*((1.0D0-A/A2)+(2*A*A*X(J)*X(J)/A2))
     &           /(2*A1*A2)
C THIS IS DXXPY
      WO(I+4,J+1)=-A*B*C*X(J)*((1.0D0-A/A1)+(2*A*A*X(I)*X(I)/A1))
     &           /(2*A1*A2)
C THIS IS DXXDYY
      WOW=1.0D0+2*(A1-A)*X(J)**2+2*(A2-A)*X(I)**2+4*AS*X(I)**2*X(J)**2
      WO(I+4,J+4)=WOW*B*C*AS/(4*A1**2*A2**2)
 1003 CONTINUE
  115 CONTINUE
      IF(DYES) GO TO 1004
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
      WO(I+1,J+K+5)=-A*A*B*C*(0.5D0-A*X(I)*X(I))*X(L)/(A1*A2*A2)
C THIS IS DXXDXY
      WO(I+4,J+K+5)=A*A*B*C*X(J)*X(K)*((1.0D0-3*A/A1)
     &             +2*A*A*X(I)*X(I)/A1)/(2*A1*A2*A2)
C THIS IS DXYPX
      WO(J+K+5,I+1)=A*A*B*C*(0.5D0-A*X(I)*X(I))*X(L)/(A1*A1*A2)
C THIS IS DXYDXX
      WO(J+K+5,I+4)=A*A*B*C*X(J)*X(K)*((1.0D0-3*A/A2)
     &             +2*A*A*X(I)*X(I)/A2)/(2*A1*A1*A2)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      WO(I+1,J+K+5)=(A**3)*B*C*X(I)*X(J)*X(K)/(A1*A2*A2)
C THIS IS DXXDYZ
      WO(I+4,J+K+5)=A*A*B*C*X(J)*X(K)*((1.0D0-A/A1)
     &             +2*A*A*X(I)*X(I)/A1)/(2*A1*A2*A2)
C THIS IS DXYPZ
      WO(J+K+5,I+1)=-(A**3)*B*C*X(I)*X(J)*X(K)/(A1*A1*A2)
C THIS IS DXYDZZ
      WO(J+K+5,I+4)=A*A*B*C*X(J)*X(K)*((1.0D0-A/A2)
     &             +2*A*A*X(I)*X(I)/A2)/(2*A1*A1*A2)
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
C THIS IS DXYDYZ
  135 WO(I+J+5,K+L+5)=-(A**3)*B*C*X(M2)*X(M3)
     &               *(1.0D0-2*A*X(M1)*X(M1))/(2*A1*A1*A2*A2)
      GO TO 145
C THIS IS DXYDXY
  140 WO(I+J+5,K+L+5)=A*A*B*C*(1.0D0-2*A*X(I)*X(I)-2*A*X(J)*X(J)
     &               +4*A*A*X(I)*X(I)*X(J)*X(J))/(4*A1*A1*A2*A2)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      GO TO 200
  150 CONTINUE
C THIS SS
      WO(1,1)=B
      DO 155 I=1,3
      IF(DYES) GO TO 1005
C THIS IS SDXX
      WO(1,I+4)=B*(1.0D0-A/A2)/(2*A2)
C THIS IS DXXS
      WO(I+4,1)=B*(1.0D0-A/A1)/(2*A1)
 1005 CONTINUE
C THIS IS PXPX
      WO(I+1,I+1)=0.5D0*A*B/(A1*A2)
      IF(DYES) GO TO 1006
C THIS IS DXXDXX
      WO(I+4,I+4)=3*A**2*B/(4*A1**2*A2**2)
 1006 CONTINUE
  155 CONTINUE
      IF(DYES) GO TO 1007
      DO 165 I=1,3
      DO 165 J=1,3
      IF (I.EQ.J) GO TO 165
C THIS IS DXXDYY
      WO(I+4,J+4)=A*A*B/(4*A1**2*A2**2)
C THIS IS DXYDXY
      WO(I+J+5,I+J+5)=WO(I+4,J+4)
  165 CONTINUE
 1007 CONTINUE
  200 CONTINUE
      RETURN
      END
