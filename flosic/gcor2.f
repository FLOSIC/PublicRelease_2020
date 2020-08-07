C UTEP Electronic Structure Lab (2020)
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,RTRS,GG,GGRS)
c slimmed down version of GCOR used in PW91 routines, to interpolate
c LSD correlation energy, as given by (10) of
c J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c K. Burke, May 11, 1996.
      IMPLICIT REAL*8 (A-H,O-Z)
      SAVE
      Q0 = -2.D0*A*(1.D0+A1*RTRS*RTRS)
      Q1 = 2.D0*A*RTRS*(B1+RTRS*(B2+RTRS*(B3+B4*RTRS)))
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RTRS+2.D0*B2+RTRS*(3.D0*B3+4.D0*B4*RTRS))
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1*(1.D0+Q1))
      RETURN
      END
