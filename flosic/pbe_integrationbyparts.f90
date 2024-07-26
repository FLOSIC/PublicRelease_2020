! Source code is from:
! https://github.com/QEF/q-e/blob/master/Modules/more_functionals.f90
! Edited and adjusted for our code.
!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!______________________________________________________________________
      subroutine exchpbe(rho,agrad,ex,dexdrho,dexdg)
      !subroutine exchpbe(rho,agrad,exloc,exnlc,dexdrho,dexdg)
!     _________________________________________________________________
!
! Perdew-Burke-Ernzerhof gga, Exchange term:
! Calculates the exchange energy density and the two functional derivative
! that will be used to calculate the potential
!
!     USE kinds, ONLY: DP
      implicit none
! input
! input rho:     charge density
! input agrad:   abs(grad rho)
      real(8) rho, agrad
! ouput
! output ex: Ex[rho,grad_rho] = \int ex dr
! output dexdrho: d ex / d rho
! output dexdg:   d ex / d grad_rho(i) = dexdg*grad_rho(i)/abs(grad_rho)
      real(8) ex, dexdrho, dexdg
! local
      real(8) thrd, thrd4, pi32td, ax, al, um, uk, ul
      parameter(thrd=.33333333333333333333d0,thrd4=4.d0/3.d0)
      parameter(pi32td=3.09366772628014d0) ! pi32td=(3.d0*pi*pi)**0.333d0
      parameter(al=0.161620459673995d0)    ! al=1.0/(2.0*(pi32)**0.333d0)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
!
      real(8) rhothrd, exunif, dexunif, kf, s, s2, p0, fxpbe, fs
!----------------------------------------------------------------------
! construct LDA exchange energy density
!
      rhothrd = rho**thrd
      dexunif = ax*rhothrd
      exunif  = rho*dexunif
!----------------------------------------------------------------------
! construct PBE enhancement factor
!
      kf = pi32td*rhothrd
      s = agrad/(2.d0*kf*rho)
      s2 = s*s
      p0 = 1.d0 + ul*s2
      fxpbe = 1.d0 + uk - uk/p0
      !ex = exunif*fxpbe
      ex = dexunif*fxpbe !Edited to make it Ex per density
      !exloc = dexunif
      !exnlc = dexunif*(fxpbe-1.0d0)
!----------------------------------------------------------------------
! now calculates the potential terms
!
!  fs=(1/s)*d fxPBE/ ds
!
      fs=2.d0*uk*ul/(p0*p0)
      dexdrho = dexunif*thrd4*(fxpbe-s2*fs)
      dexdg = ax*al*s*fs
!
      return
      end subroutine exchpbe

!----------------------------------------------------------------------
      subroutine ecorpbe(rho,agrad,zet,ectot,decup,decdn,decdg,nspin)
      !subroutine ecorpbe(rho,agrad,zet,ec,h,decup,decdn,decdg,nspin)
!     -----------------------------------------------------------------
!
!  Adapted from the Official PBE correlation code. K. Burke, May 14, 1996.
!
!   input: rho   = rho_up + rho_down; total  charge density
!   input: agrad = abs( grad(rho) )
!   input: zet   = (rho_up-rho_down)/rho
!   input: nspin
!  output: ectot = ec*rho       ---correlation energy density---
!  output: decup = d ( ec*rho ) / d (rho_up)
!  output: decdn = d ( ec*rho ) / d (rho_down)
!  output: decdg = (d ( ec*rho ) / d (grad(rho)_i)) * agrad / grad_i
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!     USE kinds, ONLY: DP
!     USE constants, ONLY: pi
      implicit none
      real(8) rho, agrad, zet, ectot, decup, decdn, decdg
      integer nspin
      real(8) pi32, alpha, thrd, thrdm, thrd2, sixthm, thrd4,  &
          gam, fzz, gamma, bet, delt, eta, pi
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      gam= 2^(4/3)-2
!      fzz=f''(0)= 8/(9*gam)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      parameter(pi=3.14159265358979323846d0)
      parameter(pi32=29.608813203268075856503472999628d0)
      parameter(alpha=1.91915829267751300662482032624669d0)
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(gam=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*gam))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
      real(8) g, fk, rs, sk, twoksg, t
      real(8) rtrs, eu, eurs, ep, eprs, alfm, alfrsm, z4, f, ec
      real(8) ecrs, fz, eczet, comm, vcup, vcdn, g3, pon, b, b2, t2, t4
      real(8) q4, q5, h, g4, t6, rsthrd, gz, fac
      real(8) bg, bec, q8, q9, hb, hrs, hz, ht, pref
!----------------------------------------------------------------------
      if (nspin.eq.1) then
         g=1.d0
      else
         g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)*0.5d0
      endif
      fk=(pi32*rho)**thrd
      rs=alpha/fk
      sk=sqrt(4.d0*fk/pi)
      twoksg=2.d0*sk*g
      t=agrad/(twoksg*rho)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! eu=unpolarized LSD correlation energy
! eurs=deu/drs
! ep=fully polarized LSD correlation energy
! eprs=dep/drs
! alfm=-spin stiffness, [c](3).
! alfrsm=-dalpha/drs
! f=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      call gcor2(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,      &
     &    0.49294d0,rtrs,eu,eurs)
      if (nspin.eq.2) then
         call gcor2(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0, &
     &       0.62517d0,rtrs,ep,eprs)
         call gcor2(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,  &
     &       0.49671d0,rtrs,alfm,alfrsm)
         z4 = zet**4
         f=((1.d0+zet)**thrd4+(1.d0-zet)**thrd4-2.d0)/gam
         ec = eu*(1.d0-f*z4)+ep*f*z4-alfm*f*(1.d0-z4)/fzz
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ecrs = dec/drs [c](A2)
! eczet=dec/dzeta [c](A3)
! fz = df/dzeta [c](A4)
         ecrs = eurs*(1.d0-f*z4)+eprs*f*z4-alfrsm*f*(1.d0-z4)/fzz
         fz = thrd4*((1.d0+zet)**thrd-(1.d0-zet)**thrd)/gam
         eczet = 4.d0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu       &
     &           -(1.d0-z4)*alfm/fzz)
         comm = ec -rs*ecrs/3.d0-zet*eczet
         vcup = comm + eczet
         vcdn = comm - eczet
      else
         ecrs = eurs
         ec = eu
         vcup = ec -rs*ecrs/3.d0
      endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! PBE correlation energy
! g=phi(zeta), given after [a](3)
! delt=bet/gamma
! b=a of [a](8)
!      g=((1.d0+zet)**thrd2+(1.d0-zet)**thrd2)/2.d0
      g3 = g**3
      pon=-ec/(g3*gamma)
      b = delt/(dexp(pon)-1.d0)
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      q4 = 1.d0+b*t2
      q5 = 1.d0+b*t2+b2*t4
      h = g3*(bet/delt)*dlog(1.d0+delt*Q4*t2/Q5)
      !ectot = rho*(ec + h)
      ectot = ec + h ! Edited to make it Ec per density
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! energy done. Now the potential, using appendix e of [b].
      t6 = t4*t2
      rsthrd = rs/3.d0
      fac = delt/b+1.d0
      bec = b2*fac/(bet*g3)
      q8 = q5*q5+delt*q4*q5*t2
      q9 = 1.d0+2.d0*b*t2
      hb = -bet*g3*b*t6*(2.d0+b*t2)/q8
      hrs = -rsthrd*hb*bec*ecrs
      ht = 2.d0*bet*g3*q9/q8
      comm = h+hrs-7.d0*t2*ht/6.d0
      if (nspin.eq.2) then
         g4 = g3*g
         bg = -3.d0*b2*ec*fac/(bet*g4)
         gz=(((1.d0+zet)**2+eta)**sixthm-                               &
     &   ((1.d0-zet)**2+eta)**sixthm)/3.d0
         hz = 3.d0*gz*h/g + hb*(bg*gz+bec*eczet)
         pref = hz-gz*t2*ht/g
         decup = vcup + comm + pref*(  1.d0 - zet)
         decdn = vcdn + comm + pref*( -1.d0 - zet)
      else
         decup = vcup + comm
      endif
      decdg = t*ht/twoksg
!
      return
      end subroutine ecorpbe
!______________________________________________________________________

