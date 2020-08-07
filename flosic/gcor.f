C UTEP Electronic Structure Lab (2020)
C
C ******************************************************************
C
      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,RS,GG,GGRS)
C
C  CALLED BY SUBROUTINE PW91LC
C
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
      Q0 = -2*A*(1.0D0+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      Q1 = 2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RS)
      Q2 = LOG(1.0D0+1.0D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2*B2+3*B3*RS12+4*B4*RS)
      GGRS = -2*A*A1*Q2-Q0*Q3/(Q1*Q1+Q1)
      RETURN
      END
