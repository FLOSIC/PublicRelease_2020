C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE KNMXSF(A1,A2,R,T,WK)
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL DYES
      DIMENSION R(3),T(3),X(3),WK(ND,ND)
      SAVE
      DATA DYES/.FALSE./
      X(1)=T(1)-R(1)
      X(2)=T(2)-R(2)
      X(3)=T(3)-R(3)
      A=A1*A2/(A1+A2)
      AS=A*A
      AC=AS*A
      PI=3.14159265358979324D0
      B=PI/(A1+A2)
      B=SQRT(B)
      B=B**3
      DO 100 I=1,ND
      DO 100 J=1,ND
  100 WK(I,J)=0.0D0
      XS=X(1)**2+X(2)**2+X(3)**2
      IF( XS.LT.1.0D-14 ) GO TO 150
      C=A*XS
      C=EXP(-C)
C THIS IS SS
      WK(1,1)=A*B*C*(3.0D0-2*A*XS)
      DO 105 I=1,3
C THIS IS SPX
      WK(1,I+1)=-AS*B*C*(5.0D0-2*A*XS)*X(I)/A2
C THIS IS PXS
      WK(I+1,1)= AS*B*C*(5.0D0-2*A*XS)*X(I)/A1
      IF(DYES) GO TO 1000
C THIS IS SDXX
      WK(1,I+4)=0.5D0*A*B*C*(3*A2-5*A+2*A*(A-A2)*XS+14*AS*X(I)**2
     &         -4*AC*XS*X(I)**2)/(A2**2)
C THIS IS DXXS
      WK(I+4,1)=0.5D0*A*B*C*(3*A1-5*A+2*A*(A-A1)*XS+14*AS*X(I)**2
     &         -4*AC*XS*X(I)**2)/(A1**2)
 1000 CONTINUE
C THIS IS PXPX
      WK(I+1,I+1)=AS*B*C*(5.0D0-2*A*XS-14*A*X(I)**2+4*AS*XS*X(I)**2) 
     &           /(2*A1*A2)
      IF(DYES) GO TO 1001
C THIS IS DXXPX
      WK(I+4,I+1)=-AS*B*C*(5*A1-21*A+18*AS*X(I)**2+2*A*(3*A-A1)*XS
     &           -4*AC*XS*X(I)**2)*X(I)/(2*A1**2*A2)
C THIS IS PXDXX
      WK(I+1,I+4)= AS*B*C*(5*A2-21*A+18*AS*X(I)**2+2*A*(3*A-A2)*XS
     &           -4*AC*XS*X(I)**2)*X(I)/(2*A1*A2**2)
C THIS IS DXXDXX
      WKW=2*A1*A2-21*AS-(14*A*A1*A2-108*AC)*X(I)**2+6*AC*XS
     &   +(4*AS*A1*A2-24*A**4)*XS*X(I)**2-44*A**4*X(I)**4
     &   +8*A**5*XS*X(I)**4
      WK(I+4,I+4)=-A*B*C*WKW/(4*A1**2*A2**2)
 1001 CONTINUE
  105 CONTINUE
      IF(DYES) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      WK(1,K)=AC*B*C*X(I)*X(J)*(7.0D0-2*A*XS)/(A2**2)
C THIS IS DXYS
      WK(K,1)=AC*B*C*X(I)*X(J)*(7.0D0-2*A*XS)/(A1**2)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      WK(I+1,J+1)=-AC*B*C*X(I)*X(J)*(7.0D0-2*A*XS)/(A1*A2)
      IF(DYES) GO TO 1003
C THIS IS PXDYY
      WK(I+1,J+4)= AS*B*C*X(I)*(5*A2-7*A-2*A*(A2-A)*XS+18*AS*X(J)**2
     &           -4*AC*XS*X(J)**2)/(2*A1*A2**2)
C THIS IS DXXPY
      WK(I+4,J+1)=-AS*B*C*X(J)*(5*A1-7*A-2*A*(A1-A)*XS+18*AS*X(I)**2
     &           -4*AC*XS*X(I)**2)/(2*A1**2*A2)
C THIS IS DXXDYY
      WKW=7*AS-2*A1*A2+2*AS*(7*A2-9*A)*X(I)**2+2*AS*(7*A1-9*A)*X(J)**2
     &   -2*AC*XS+4*AC*(A-A2)*XS*X(I)**2+4*AC*(A-A1)*XS*X(J)**2
     &   +44*A**4*X(I)**2*X(J)**2-8*A**5*XS*X(I)**2*X(J)**2
      WK(I+4,J+4)=A*B*C*WKW/(4*A1**2*A2**2)
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
      WK(I+1,J+K+5)=-AC*B*C*(7.0D0-18*A*X(I)**2
     &             -2*A*XS+4*AS*XS*X(I)**2)*X(L)/(2*A1*A2**2)
C THIS IS DXYPX
      WK(J+K+5,I+1)= AC*B*C*(7.0D0-18*A*X(I)**2
     &             -2*A*XS+4*AS*XS*X(I)**2)*X(L)/(2*A1**2*A2)
C THIS IS DXXDXY
      WK(I+4,J+K+5)=AC*B*C*(7*A1-27*A+22*AS*X(I)**2-2*A*(A1-3*A)*XS
     &             -4*AC*XS*X(I)**2)*X(I)*X(L)/(2*A1**2*A2**2)
C THIS IS DXYDXX
      WK(J+K+5,I+4)=AC*B*C*(7*A2-27*A+22*AS*X(I)**2-2*A*(A2-3*A)*XS
     &             -4*AC*XS*X(I)**2)*X(I)*X(L)/(2*A1**2*A2**2)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      WK(I+1,J+K+5)= A**4*B*C*X(I)*X(J)*X(K)*(9.0D0-2*A*XS)/(A1*A2**2) 
C THIS IS DYZPX
      WK(J+K+5,I+1)=-A**4*B*C*X(I)*X(J)*X(K)*(9.0D0-2*A*XS)/(A1**2*A2)
C THIS IS DXXDYZ
      WK(I+4,J+K+5)=AC*B*C*X(J)*X(K)*(7*A1-9*A+22*AS*X(I)**2
     &             +2*A*(A-A1)*XS-4*AC*XS*X(I)**2)/(2*A1**2*A2**2)
C THIS IS DYZDXX
      WK(J+K+5,I+4)=AC*B*C*X(J)*X(K)*(7*A2-9*A+22*AS*X(I)**2
     &             +2*A*(A-A2)*XS-4*AC*XS*X(I)**2)/(2*A1**2*A2**2)
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
  135 CONTINUE
      WK(I+J+5,K+L+5)=-A**4*B*C*(9.0D0-2*A*XS-22*A*X(M1)**2
     &               +4*AS*XS*X(M1)**2)*X(M2)*X(M3)/(2*A1**2*A2**2)
      GO TO 145
C THIS IS DXYDXY
  140 CONTINUE
      WK(I+J+5,K+L+5)= AC*B*C*(7.0D0-2*A*XS+(4*AS*XS-18*A)*(X(I)**2
     &               +X(J)**2)+(44*AS-8*AC*XS)*X(I)**2*X(J)**2)
     &               /(4*A1**2*A2**2)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      GO TO 200
  150 CONTINUE
C THIS SS
      WK(1,1)=3.0D0*A*B
      DO 155 I=1,3
      IF(DYES) GO TO 1005
C THIS IS SDXX
      WK(1,I+4)=0.5D0*A*B*(3.0D0-5*A/A2)/A2
C THIS IS DXXS
      WK(I+4,1)=0.5D0*A*B*(3.0D0-5*A/A1)/A1
 1005 CONTINUE
C THIS IS PXPX
      WK(I+1,I+1)= 5*B*AS/(2*A1*A2)
      IF(DYES) GO TO 1006
C THIS IS DXXDXX
      WK(I+4,I+4)=-B*A*(2.0D0-21*AS/(A1*A2))/(4*A1*A2)
 1006 CONTINUE
  155 CONTINUE
      IF(DYES) GO TO 1007
      DO 165 I=1,3
      DO 165 J=1,3
      IF (I.EQ.J) GO TO 165
C THIS IS DXXDYY
      WK(I+4,J+4)=-B*A*(2.0D0-7*AS/(A1*A2))/(4*A1*A2)
  165 CONTINUE
      DO 175 I=8,10
C THIS IS DXYDXY
  175 WK(I,I)=7*A*AS*B/(4*A1**2*A2**2)
 1007 CONTINUE
  200 CONTINUE
      RETURN
      END
