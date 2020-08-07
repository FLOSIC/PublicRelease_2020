C UTEP Electronic Structure Lab (2020)
***********************************************************************c
*  Becke exchange for a spin-unpolarized electronic system 
*
*  Gradient-corrected exchange energy based on
*     [A.D. Becke, J.Chem.Phys.96, 2155, 1992].
*  The LSDA energy functional, obtained as E{n+,n-}=(E{2n+}+E{2n-})/2,
*     and the functional derivative formula are given by 
c     [J.P. Perdew , PRB 33, 8800, 1986].
*     [J.P. Perdew , PRB 34, 7406, 1986].
*  see also [G.Ortiz ...,PRB 43, 6376 (1991)] eq. (A2)
*  
*  Hartree a.u.
*
*  Inputs
*  d            density
*  s            abs(grad d)/(2kf*d)
*  u            (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
*           >>  grad(abs(grad d) has mixed derivatives ! <<
*  v            (laplacian d)/(d*(2*kf)**2)
*
*  Outputs
*  ex           exchange energy per electron
*  vx           exchange potential
*
*                  Adapted Raja.  June 1998.
c**********************************************************************
c
      subroutine xbecke(d,s,u,v,ex,dex,vx)
c
      implicit real*8 (a-h,o-z)
      data c / .779555417944150792d1/
      data b / .42d-2/
      data bb/-.451357747124625192d-2/
      data ax/-.738558766382022406d0/
      data thrd,thrd4/.333333333333333333d0,.1333333333333333333d1/
      intrinsic dlog, dsqrt
c
c exchange enhancement factor f
      x  = c*s
      y1 = 1.d0/dsqrt(1.d0+x*x)
      y0 = dlog(x+1.d0/y1)
      y2 = -x*y1*y1*y1
      ddi= 1.d0/(1.d0 + 6.d0*b*x*y0)
      dd1= 6.d0*b*(y0+x*y1)
      g  = 1.d0 - 0.5d0*x*dd1*ddi
      fs = -2.d0*bb*c*c*ddi
      g1 = -3.d0*b*(y0+x*(3.d0*y1+x*y2-dd1*dd1*ddi/(6.d0*b)))
      fss= fs*c*(g1 - g*dd1)*ddi
      fs = fs*g
      f  = 1.d0 - bb*x*x*ddi
c
c LDA only
      fac= ax*d**thrd
C Local contribution
      ex=fac
C Non Local contribution
      dex=(f-1.d0)*fac
c
c energy
C      ex = fac*f
c
c potential
      vx = fac*(thrd4*f-(u-thrd4*s*s*s)*fss-v*fs)
c
      return
      end
*e
*

