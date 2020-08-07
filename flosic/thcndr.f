C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE THCNDR(A1,A2,A3,B,C,DR)
      PARAMETER (ND=10)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL NOD
      DIMENSION F(3),AE(3),BE(3),DR(ND,ND),B(3),C(3)
      SAVE
      DATA PI/3.141592653589793D0/
      DATA NOD/.FALSE./
      DO 5 K=1,ND
      DO 5 J=1,ND
    5 DR(J,K)=0.0D0
      AT=A1+A2+A3
      EX=(A2*B(1)+A3*C(1))/AT
      EY=(A2*B(2)+A3*C(2))/AT
      EZ=(A2*B(3)+A3*C(3))/AT
      E2=EX**2+EY**2+EZ**2
      B2=B(1)**2+B(2)**2+B(3)**2
      C2=C(1)**2+C(2)**2+C(3)**2
      CCC=AT*E2-A2*B2-A3*C2
      COEF=2.0D0*PI*EXP(CCC)/AT
      F(1)=(A2*B(1)-(A1+A2)*C(1))/AT
      F(2)=(A2*B(2)-(A1+A2)*C(2))/AT
      F(3)=(A2*B(3)-(A1+A2)*C(3))/AT
      F2=F(1)**2+F(2)**2+F(3)**2
      AE(1)=EX
      AE(2)=EY
      AE(3)=EZ
      BE(1)=EX-B(1)
      BE(2)=EY-B(2)
      BE(3)=EZ-B(3)
      XX=AT*F2
      CALL FMTCAL(XX,S0,S2,S4,S6,S8)
C THIS SS
      DR(1,1)=COEF*S0
      DO 105 I=1,3
C THIS IS SPX
      DR(1,I+1)=COEF*(BE(I)*S0-F(I)*S2)
C THIS IS PXS
      DR(I+1,1)=COEF*(AE(I)*S0-F(I)*S2)
      IF(NOD) GO TO 1000
C THIS IS SDXX
      DR(1,I+4)=COEF*((BE(I)**2+0.5D0/AT)*S0-(2*BE(I)*F(I)
     &         +0.5D0/AT)*S2+F(I)**2*S4)
C THIS IS DXXS
      DR(I+4,1)=COEF*((AE(I)**2+0.5D0/AT)*S0-(2*AE(I)*F(I)
     &         +0.5D0/AT)*S2+F(I)**2*S4)
 1000 CONTINUE
C THIS IS PXPX
      DR(I+1,I+1)=COEF*((AE(I)*BE(I)+0.5D0/AT)*S0-(AE(I)*F(I)
     &           +BE(I)*F(I)+0.5D0/AT)*S2+F(I)**2*S4)
      IF(NOD) GO TO 1001
C THIS IS DXXPX
      DP=(AE(I)**2*BE(I)+0.5D0*(2*AE(I)+BE(I))/AT)*S0
      DP=DP-(AE(I)**2*F(I)+0.5D0*(3*F(I)+2*AE(I)+BE(I))/AT
     &  +2*AE(I)*BE(I)*F(I))*S2
      DP=DP+(F(I)**2*(2*AE(I)+BE(I))+1.5D0*F(I)/AT)*S4-F(I)**3*S6
      DR(I+4,I+1)=COEF*DP
C THIS IS PXDXX
      PD=(AE(I)*BE(I)**2+0.5D0*(AE(I)+2*BE(I))/AT)*S0
      PD=PD-(BE(I)**2*F(I)+0.5D0*(3*F(I)+2*BE(I)+AE(I))/AT
     &  +2*AE(I)*BE(I)*F(I))*S2
      PD=PD+(F(I)**2*(AE(I)+2*BE(I))+1.5D0*F(I)/AT)*S4-F(I)**3*S6
      DR(I+1,I+4)=COEF*PD
C THIS IS DXXDXX
      DD=(AE(I)**2*BE(I)**2+0.5D0*(AE(I)**2+BE(I)**2+4*AE(I)*BE(I))/AT
     &  +0.75D0/AT**2)*S0
      DD=DD-((2*AE(I)*BE(I)+3.0D0/AT)*(AE(I)+BE(I))*F(I)
     &  +0.5D0*(AE(I)**2+BE(I)**2+4*BE(I)*AE(I))/AT+1.5D0/AT**2)*S2
      DD=DD+((AE(I)**2+BE(I)**2+4*AE(I)*BE(I))*F(I)**2
     &  +3*F(I)*(AE(I)+BE(I)+F(I))/AT+0.75D0/AT**2)*S4
      DD=DD-(2*F(I)**3*(AE(I)+BE(I))+3*F(I)**2/AT)*S6+F(I)**4*S8
      DR(I+4,I+4)=COEF*DD
 1001 CONTINUE
  105 CONTINUE
      IF(NOD) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      DR(1,K)=COEF*(BE(I)*BE(J)*S0-(BE(J)*F(I)+BE(I)*F(J))*S2+F(I)
     &       *F(J)*S4)
C THIS IS DXYS
      DR(K,1)=COEF*(AE(I)*AE(J)*S0-(AE(J)*F(I)+AE(I)*F(J))*S2+F(I)
     &       *F(J)*S4)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      DR(I+1,J+1)=COEF*(AE(I)*BE(J)*S0-(AE(I)*F(J)+BE(J)*F(I))*S2
     &           +F(I)*F(J)*S4)
      IF(NOD) GO TO 1003
C THIS IS PXDYY
      PD=(BE(J)**2+0.5D0/AT)*AE(I)*S0-F(J)**2*F(I)*S6
      PD=PD-(2*BE(J)*F(J)*AE(I)+0.5D0*AE(I)/AT+BE(J)**2*F(I)
     &  +0.5D0*F(I)/AT)*S2
      PD=PD+(F(J)**2*AE(I)+2*BE(J)*F(J)*F(I)+0.5D0*F(I)/AT)*S4
      DR(I+1,J+4)=COEF*PD
C THIS IS DXXPY
      DP=(AE(I)**2+0.5D0/AT)*BE(J)*S0-F(I)**2*F(J)*S6
      DP=DP-(2*AE(I)*F(I)*BE(J)+0.5D0*BE(J)/AT+AE(I)**2*F(J)
     &  +0.5D0*F(J)/AT)*S2
      DP=DP+(F(I)**2*BE(J)+2*AE(I)*F(I)*F(J)+0.5D0*F(J)/AT)*S4
      DR(I+4,J+1)=COEF*DP
C THIS IS DXXDYY
      DD=(AE(I)**2+0.5D0/AT)*(BE(J)**2+0.5D0/AT)*S0
      DD=DD-(2*AE(I)*BE(J)*(BE(J)*F(I)+AE(I)*F(J))+0.5D0*(AE(I)**2
     &  +BE(J)**2)/AT+(AE(I)*F(I)+BE(J)*F(J))/AT+0.5D0/AT**2)*S2
      DD=DD+(F(I)**2*BE(J)**2+F(J)**2*AE(I)**2+4*AE(I)*BE(J)*F(I)
     &  *F(J)+(AE(I)*F(I)+BE(J)*F(J))/AT+0.5D0*(F(I)**2+F(J)**2)/AT
     &  +0.25D0/AT**2)*S4
      DD=DD-(2*F(I)*F(J)*(BE(J)*F(I)+AE(I)*F(J))
     &  +0.5D0*(F(I)**2+F(J)**2)/AT)*S6
      DD=DD+F(I)**2*F(J)**2*S8
      DR(I+4,J+4)=COEF*DD
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
      PD=(BE(L)*S0-F(L)*S2)*(AE(I)*BE(I)+0.5D0/AT)
      PD=PD+(F(L)*S4-BE(L)*S2)*(AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)
      PD=PD+(BE(L)*S4-F(L)*S6)*F(I)**2
      DR(I+1,J+K+5)=COEF*PD
C THIS IS DXXDXY
      DD=(BE(L)*S0-F(L)*S2)*(AE(I)**2*BE(I)+0.5D0*(2*AE(I)+BE(I))/AT)
      DD=DD+(F(L)*S4-BE(L)*S2)*(AE(I)**2*F(I)+0.5D0*(3*F(I)
     &  +2*AE(I)+BE(I))/AT+2*AE(I)*BE(I)*F(I))
      DD=DD+(BE(L)*S4-F(L)*S6)*(F(I)**2*(2*AE(I)+BE(I))+1.5D0*F(I)/AT)
      DD=DD+(F(L)*S8-BE(L)*S6)*F(I)**3
      DR(I+4,J+K+5)=COEF*DD
C THIS IS DXYPX
      DP=(AE(L)*S0-F(L)*S2)*(AE(I)*BE(I)+0.5D0/AT)
      DP=DP+(F(L)*S4-AE(L)*S2)*(AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)
      DP=DP+(AE(L)*S4-F(L)*S6)*F(I)**2
      DR(J+K+5,I+1)=COEF*DP
C THIS IS DXYDXX
      DD=(AE(L)*S0-F(L)*S2)*(AE(I)*BE(I)**2+0.5D0*(AE(I)+2*BE(I))/AT)
      DD=DD+(F(L)*S4-AE(L)*S2)*(BE(I)**2*F(I)+0.5D0*(3*F(I)
     &  +2*BE(I)+AE(I))/AT+2*AE(I)*BE(I)*F(I))
      DD=DD+(AE(L)*S4-F(L)*S6)*(F(I)**2*(AE(I)+2*BE(I))+1.5D0*F(I)/AT)
      DD=DD+(F(L)*S8-AE(L)*S6)*F(I)**3
      DR(J+K+5,I+4)=COEF*DD
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      PD=AE(I)*BE(J)*BE(K)*S0-(AE(I)*(BE(J)*F(K)+BE(K)*F(J))
     &  +F(I)*BE(J)*BE(K))*S2
      PD=PD+(AE(I)*F(J)*F(K)+BE(K)*F(I)*F(J)+BE(J)*F(I)*F(K))*S4
     &  -F(I)*F(J)*F(K)*S6
      DR(I+1,J+K+5)=COEF*PD
C THIS IS DXXDYZ
      DD=(BE(K)*S0-F(K)*S2)*((AE(I)**2+0.5D0/AT)*BE(J))
      DD=DD+(F(K)*S4-BE(K)*S2)*(2*AE(I)*F(I)*BE(J)+0.5D0*BE(J)/AT
     &  +AE(I)**2*F(J)+0.5D0*F(J)/AT)
      DD=DD+(BE(K)*S4-F(K)*S6)*(F(I)**2*BE(J)+2*AE(I)*F(I)*F(J)
     &  +0.5D0*F(J)/AT)
      DD=DD+(F(K)*S8-BE(K)*S6)*(F(I)**2*F(J))
      DR(I+4,J+K+5)=COEF*DD
C THIS IS DYZPX
      DP=AE(J)*AE(K)*BE(I)*S0-(F(J)*AE(K)*BE(I)+AE(J)*F(K)*BE(I)
     &  +AE(J)*AE(K)*F(I))*S2
      DP=DP+(F(J)*F(K)*BE(I)+F(J)*AE(K)*F(I)+AE(J)*F(K)*F(I))*S4
      DP=DP-F(J)*F(K)*F(I)*S6
      DR(J+K+5,I+1)=COEF*DP
C THIS IS DYZDXX
      DD=(AE(K)*S0-F(K)*S2)*(BE(I)**2+0.5D0/AT)*AE(J)
      DD=DD+(F(K)*S4-AE(K)*S2)*(2*BE(I)*F(I)*AE(J)+0.5D0*AE(J)/AT
     &  +BE(I)**2*F(J)+0.5D0*F(J)/AT)
      DD=DD+(AE(K)*S4-F(K)*S6)*(F(I)**2*AE(J)+2*BE(I)*F(I)*F(J)
     &  +0.5D0*F(J)/AT)
      DD=DD+(F(K)*S8-AE(K)*S6)*F(I)**2*F(J)
      DR(J+K+5,I+4)=COEF*DD
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
      DD=(BE(M3)*S0-F(M3)*S2)*(AE(M1)*BE(M1)+0.5D0/AT)*AE(M2)
      DD=DD+(F(M3)*S4-BE(M3)*S2)*((AE(M1)*F(M1)+BE(M1)*F(M1)
     &  +0.5D0/AT)*AE(M2)+(AE(M1)*BE(M1)+0.5D0/AT)*F(M2))
      DD=DD+(BE(M3)*S4-F(M3)*S6)*(F(M1)**2*AE(M2)+(AE(M1)*F(M1)
     &  +BE(M1)*F(M1)+0.5D0/AT)*F(M2))
      DD=DD+(F(M3)*S8-BE(M3)*S6)*F(M2)*F(M1)**2
      DR(I+J+5,K+L+5)=COEF*DD
      GO TO 145
  140 CONTINUE
C THIS IS DXYDXY
      DD=(AE(I)*BE(I)+0.5D0/AT)*(AE(J)*BE(J)+0.5D0/AT)*S0
      DD=DD-((AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)*(AE(J)*BE(J)+0.5D0/AT)
     &  +(AE(J)*F(J)+BE(J)*F(J)+0.5D0/AT)*(AE(I)*BE(I)+0.5D0/AT))*S2
      DD=DD+((AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)*(AE(J)*F(J)+BE(J)*F(J)
     &  +0.5D0/AT)+(AE(I)*BE(I)+0.5D0/AT)*F(J)**2+(AE(J)*BE(J)
     &  +0.5D0/AT)*F(I)**2)*S4
      DD=DD-((AE(I)*F(I)+BE(I)*F(I)+0.5D0/AT)*F(J)**2
     &  +(AE(J)*F(J)+BE(J)*F(J)+0.5D0/AT)*F(I)**2)*S6
      DD=DD+F(I)**2*F(J)**2*S8
      DR(I+J+5,K+L+5)=COEF*DD
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      RETURN
      END
