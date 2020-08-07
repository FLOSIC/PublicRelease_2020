C UTEP Electronic Structure Lab (2020)
C> @file diagv.f
C> Diagonalization routine used for SIC newwave
       SUBROUTINE DIAGV(matlen,N,M,FOHAM,PHIRES,IPREONLY)
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:41 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MATLEN, N, M, IPREONLY, I, ITER, J, K, L, MAXITER
       REAL*8 :: SYMBOL , FOHAM, PHIRES, CONVG, FOEIG, FOHAMAVG, FOOVR,
     & FOSC, HAMTEMP, HII, HIJ, HJJ, PHIINEW, PHIJNEW, T, TEMP
!       INCLUDE 'commons.inc'
       DIMENSION FOHAM(matlen,matlen), PHIRES(matlen,matlen)

       DIMENSION FOHAMAVG(M,M), FOOVR(M,M), FOEIG(M), FOSC(M)
       DIMENSION PHIINEW(M), PHIJNEW(M)

       PARAMETER (MAXITER=100)
       PARAMETER (CONVG=1.0d-6)
       SAVE

C      Pre-diagonalize it to generate initial guesses
C      1.Generate average hamiltonian
       FOOVR=0.0d0
       do i=1,M
        FOOVR(i,i)=1.0d0
        do j=i,M
         FOHAMAVG(i,j)=0.5d0*(FOHAM(i,j)+FOHAM(j,i))
         FOHAMAVG(j,i)=FOHAMAVG(i,j)
        end do
       end do

       call DIAGGE(M,M,FOHAMAVG,FOOVR,FOEIG,FOSC,1)

C phi_new_i = sum_j phires(i,j) phi_old_j
       PHIRES=0.0d0
       do i=1,M
        do j=1,M
         PHIRES(i,j)=FOHAMAVG(j,i)
        end do
       end do

C      Iteratively solve
       if(IPREONLY==1) return
       do iter=1,maxiter
        do i=1,N
         do j=N+1,M
          hii=0.0d0
          hij=0.0d0
          hjj=0.0d0
          do k=1,M
           do l=1,M
            hii=hii+PHIRES(i,k)*PHIRES(i,l)*FOHAM(k,l)
            hij=hij+PHIRES(i,k)*PHIRES(j,l)*FOHAM(k,l)
            hjj=hjj+PHIRES(j,k)*PHIRES(j,l)*FOHAM(k,l)
           end do
          end do

          t=0.5d0*atan(2.0d0*hij/(hii-hjj))

          PHIINEW=0.0d0
          PHIJNEW=0.0d0
          do k=1,M
           PHIINEW(k)=PHIINEW(k)+cos(t)*PHIRES(i,k)+sin(t)*PHIRES(j,k)
           PHIJNEW(k)=PHIJNEW(k)-sin(t)*PHIRES(i,k)+cos(t)*PHIRES(j,k)
          end do
          PHIRES(i,1:M)=PHIINEW(1:M) !YY Array sizes are different NDH vs. NBAS
          PHIRES(j,1:M)=PHIJNEW(1:M)
         end do
        end do

        temp=0.0d0
        do i=1,N
         do j=N+1,M
          HAMTEMP=0.0d0
          do k=1,M
           do l=1,M
            HAMTEMP=HAMTEMP+PHIRES(i,k)*PHIRES(j,l)*FOHAM(k,l)
           end do
          end do
          if(abs(HAMTEMP)>temp) temp=abs(HAMTEMP)
         end do
        end do

        if(temp<CONVG) exit
       end do

       if(temp>CONVG) then
        print*,'WARNING: DIAGV not converged! ',temp
       end if

       END SUBROUTINE DIAGV
