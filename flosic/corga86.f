C UTEP Electronic Structure Lab (2020)
************************************************************************
*
*  gradient-correction to correlation energy from
*  
*  [J.P.Perdew, PRB 33, 8822 (1986) and PRB 34, 7406 (1986)]
*
*  to correlation part of the Becke-Perdew gradient-corrected 
*  xc-functional
*
*  Input
*  d1,d2 : up/down spindensity
*  dp12..: grad(d1)*grad(d2) {* == vector product}, MUST NOT == 0
*  uu    : (grad d)*grad(abs(grad d)) , d = d1 + d2
*  vv    : laplacian d
*
*  Output
*  dec  : correction to correlation energy per electron      
*  dvcup :         - "" -            potential maj. spin
*  dvcdn :         - "" -            potential min. spin
*
*  Hartree a.u.
*
************************************************************************
*
      subroutine corga86(d1,d2,dp11,dp22,dp12,uu,vv,dec,dvcup,dvcdn)
c
      implicit real*8 (a-h,o-z)
c
      data a1,a2,a3,a4,a5,a6,a7/ 2.568d-3
     &,                           1.443307452d-2
     &,                           2.843543831d-6
     &,                           5.411317331d0
     &,                           1.816419933d-1
     &,                           1.763993811d-2
     &,                           8.12908d-4/
      data t13,t23,t43,t53,t76/	.33333333333333333d0
     &,				.66666666666666667d0
     &,				.13333333333333333d1
     &,				.16666666666666667d1
     &,				.11666666666666667d1/
      data crt2/.1587401052d1/
c
      d    = d1+d2
      dm13 = 1.d0/d**t13
      d43  = d**t43
c
c gradient expansion coefficient
      c1 = a1 + dm13*(a2+dm13*a3)
      c2 = 1.d0 + dm13*(a4+dm13*(a5+dm13*a6))
      c  = 1.667d-3 + c1/c2
*                write(6,*)  'Hello from cor86'
c
      dpnorm = sqrt(dp11+dp22+2.d0*dp12)
      dpnorm2= dpnorm*dpnorm
c
      fi = a7*dpnorm/(c*d**t76)
c
c spin interpolation
      zet= (d1-d2)/d
      dd = sqrt(.5d0*((1.d0+zet)**t53 + (1.d0-zet)**t53))
c
c dC(n)/dn
      cp = -t13/(d43*c2*c2)
     &    *(c2*(a2+2.d0*a3*dm13)-c1*(a4+dm13*(2.d0*a5+dm13*3.d0*a6)))
c
c spin-independent terms
      www=( (fi-1.d0)*(cp/c-t43/d)+(fi-2.d0)*fi*(cp/c+t76/d) )*dpnorm2
      uuu=(3.d0-fi)*fi*uu/dpnorm
      vvv=(fi-2.d0)*vv
c
c spin dependent term     
      zzz=crt2*.5d0*t53/(dd*d43)**2*(d1**t23-d2**t23)
      zz1= zzz*((1.d0-fi)*d2*dpnorm2-(2.d0-fi)*d*(dp12+dp22))
      zz2=-zzz*((1.d0-fi)*d1*dpnorm2-(2.d0-fi)*d*(dp12+dp11))
c
      dec=c/(dd*exp(fi)*d43)
c
      dvcup=dec*(www+uuu+vvv+zz1)
      dvcdn=dec*(www+uuu+vvv+zz2)
      dec=dec*dpnorm2/d
c
      return
      end
