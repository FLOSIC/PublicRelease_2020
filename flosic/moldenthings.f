C UTEP Electronic Structure Lab (2020)
!4/13/2017
!Fixed indexing bug
!L=1,2,3: S,P,D types in the power function.

!This subroutine is referenced and adapted from MODEN2AIM 
!sourcecode (https://github.com/zorkzou/Molden2AIM)

! Ordering seems compatible with NRLMOL
c-----------------------------------------------------------------------
c---  calculate the normalization factor for GTO(l,m,n)
c---  it = 1,...,35
c---  Ordering type:
c---  S,P,D: MOLDEN, Gaussian, GAMESS, WFN, ...
c---  F: WFN (for MOLDEN and Gaussian, 14~19 are different)
c---  G: WFN, MOLDEN (not for Gaussian)
c-----------------------------------------------------------------------

      function fnorm_lmn(a,it)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      pi=acos(-1.d0)
      pi3=pi**3.d0

      select case(it)
        case(1)            ! 1S
          n1=3
          n2=3
          nf=1
        case(2:4)          ! 3P: x, y, z
          n1=7
          n2=5
          nf=1
        case(5:10)         ! 6D
          n1=11
          n2=7
          select case(it)
            case(5:7)        ! xx, yy, zz
              nf=9
            case(8:10)       ! xy, xz, yz
              nf=1
          end select
        case(11:20)        ! 10F
          n1=15
          n2=9
          select case(it)
            case(11:13)      ! xxx, yyy, zzz
              nf=225
            case(14:19)      ! xxy, xxz, yyz, xyy, xzz, yzz
              nf=9
            case(20)         ! xyz
              nf=1
          end select
        case(21:35)        ! 15G
          n1=19
          n2=11
          select case(it)
            case(21:23)      ! xxxx yyyy zzzz
              nf=11025
            case(24:29)      ! xxxy xxxz yyyx yyyz zzzx zzzy
              nf=225
            case(30:32)      ! xxyy xxzz yyzz
              nf=81
            case(33:35)      ! xxyz yyxz zzxy
              nf=9
          end select
cc<<<     Gaussian (subroutine pattml should also be modified)
c          select case(it)
c            case(21,25,35)               ! xxxx yyyy zzzz
c              nf=11025
c            case(22,24,26,29,33,34)      ! xxxy xxxz yyyx yyyz zzzx zzzy
c              nf=225
c            case(23,30,32)               ! xxyy xxzz yyzz
c              nf=81
c            case(27,28,31)               ! xxyz yyxz zzxy
c              nf=9
c          end select
cc>>>
      end select

c--- Normal^4 = 2^n1 * a^n2 / (pi^3 * nf)
c---   n1=3+4*(l+m+n)
c---   n2=3+2*(l+m+n)
c---   nf=[(2l-1)!!(2m-1)!!(2n-1)!!]^2
      f = (2.d0**dble(n1)) * (a**dble(n2)) / (pi3 * dble(nf))
      fnorm_lmn=sqrt(sqrt(f))

      return
      end


!Taken from MOLDEN2AIM sourcecode
!Renormalize a contracted basis function
!       subroutine renorm(al,a,ci,ngauss)
       subroutine renorm(it,a,ci,ngauss)
       implicit double precision (a-h,o-z)
       parameter (tol=1.0d-10,maxpgc=10000)
       dimension a(*),ci(*),c(maxpgc)
!       character*1 al

       pi=acos(-1.d0)
       pi3=pi**3.d0

       print *, "Renormalizing basis function"

c--- ngauss.le.maxpgc has been checked in subroutine bknorm
 
c--- unnormalize primitives
c--- Normal^4 = 2^n1 * a^n2 / (pi^3 * nf)
c---   n1=3+4*L; n2=3+2*L, nf=[(2L-1)!!]^2
!       call power(al,n1,n2,nf)
       call power(it,n1,n2,nf)
       fc = (2.d0**dble(n1)) / (pi3 * dble(nf))
       do i = 1,ngauss
         f = fc * (a(i)**dble(n2))
         f = sqrt(sqrt(f))
         c(i) = ci(i)*f
       end do

       fsum = 0.d0
       do i = 1,ngauss
         do j = 1,i
           a2 = (a(i)+a(j))/2.d0
           f = fc * (a2**dble(n2))
           f = sqrt(f)
           f = c(i)*c(j)/f
           if (i .ne. j) f = f*2.d0
           fsum = fsum+f
         end do
       end do
       if (fsum .gt. tol) fsum = 1.d0/sqrt(fsum)
       do i = 1,ngauss
         print *, ci(i),"x",fsum,"=",ci(i)*fsum
         ci(i) = ci(i) * fsum    !This should go to ISYMGEN/INPUT
       end do
 
       return
       end


!Taken from MOLDEN2AIM sourcecode
!get n1=3+4*L, n2+3+2*L, nf (2L-1)!!^2
!       subroutine power(al,n1,n2,nf)
       subroutine power(it,n1,n2,nf)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!       character*1 al
 
       n1=0
       n2=0
       nf=0
!       select case(al)
       select case(it)
!         case('S')
         case(1)
           n1=3
           n2=3
           nf=1
!         case('P')
!         case(2:4)
         case(2)
           n1=7
           n2=5
           nf=1
!         case('D')
!         case(5:10)
         case(3)
           n1=11
           n2=7
           nf=9
!         case('F')
         case(11:20)
           n1=15
           n2=9
           nf=225
!         case('G')
         case(21:35)
           n1=19
           n2=11
           nf=11025
       end select
 
       return
       end
