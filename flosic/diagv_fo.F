C UTEP Electronic Structure Lab (2020)
C> @file diagv_fo.f
C *********************************************************************
C
       SUBROUTINE DIAGV_FO(IPREONLY)
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:41 MDT 2017

!      INCLUDE  'PARAMAS' 
       use debug1
       use global_inputs,only : M=>INBAS
       use for_diag1,only : FOHAM=>HAM
       use diagv1,only : N=>NORB,PHIRES
       INCLUDE  'PARAMA2'  
       INTEGER :: IPREONLY, ITER, MAXITER, NTID
       REAL*8 :: SYMBOL , CONVG, FOEIG, FOHAMAVG, FOOVR,
     & FOSC, HAMTEMP, HII, HIJ, HJJ, PHIINEW, PHIJNEW, T, TEMP
!       INCLUDE 'commons.inc'

!       INTEGER, intent(in) :: N
!       INTEGER, intent(in) :: M

       INTEGER I,J,K,L,matlen

!       DIMENSION FOHAM(matlen,matlen), PHIRES(matlen,matlen)
!       DIMENSION FOHAM(matlen,matlen), PHIRES(M,M)

       DIMENSION FOHAMAVG(M,M), FOOVR(M,M), FOEIG(M), FOSC(M)
       DIMENSION PHIINEW(M), PHIJNEW(M)

       REAL(8) :: H1(N,M), H2(N,M-N) !dgemm work arrays
       REAL(8) :: alpha,beta,phii,phij,sum1
       REAL(8) :: time1,time2,time3,time4
       CHARACTER(1) :: transT,transn
       REAL(8),external :: ddot
       INTEGER :: N2, ierr

       PARAMETER (MAXITER=600)
       PARAMETER (CONVG=1.0d-6)
       SAVE

       matlen=NDH
C BLAS SETUP
       transT='T'
       transn='N'
       alpha=1.0d0
       beta= 0.0d0
       N2=M-N !N unoccupied PHIRES

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

       call gttime(time1)
#ifdef SCALAPACK
#else
       call DIAGGE(M,M,FOHAMAVG,FOOVR,FOEIG,FOSC,1)
#endif
       call gttime(time2)
       call timout('DIAGGE',time2-time1)

C phi_new_i = sum_j phires(i,j) phi_old_j
       PHIRES=0.0d0
       do i=1,M
        do j=1,M
         PHIRES(j,i)=FOHAMAVG(j,i) !transposing PHIRES later for better access
!         PHIRES(i,j)=FOHAMAVG(j,i)
        end do
       end do
       if(IPREONLY==1) return
       call gttime(time3)
!      Iteratively solve
       do iter=1,maxiter
        do i=1,N
          do j=N+1,M
            call dgemv(transn,M,M,alpha,FOHAM,NDH, 
     &          PHIRES(1,i),1,beta,FOSC,1)
            hii=ddot(M,PHIRES(1,i),1,FOSC,1)
            call dgemv(transn,M,M,alpha,FOHAM,NDH,
     &             PHIRES(1,j),1,beta,FOSC,1)
            hij=ddot(M,PHIRES(1,i),1,FOSC,1)
            hjj=ddot(M,PHIRES(1,j),1,FOSC,1)
! ORIGINAL CODE
!            hii=0.0d0
!            hij=0.0d0
!            hjj=0.0d0
!            do l=1,M
!              do k=1,M
!                hii=hii+PHIRES(k,i)*PHIRES(l,i)*FOHAM(k,l)
!                hij=hij+PHIRES(k,i)*PHIRES(l,j)*FOHAM(k,l)
!                hjj=hjj+PHIRES(k,j)*PHIRES(l,j)*FOHAM(k,l)
!              end do
!            end do
            t=0.5d0*atan(2.0d0*hij/(hii-hjj))

            PHIINEW=0.0d0
            PHIJNEW=0.0d0
            do k=1,M
              PHIINEW(k)=PHIINEW(k)+cos(t)*PHIRES(k,i)+
     &                              sin(t)*PHIRES(k,j)
              PHIJNEW(k)=PHIJNEW(k)-sin(t)*PHIRES(k,i)+
     &                              cos(t)*PHIRES(k,j)
            end do
            PHIRES(1:M,i)=PHIINEW(1:M) !YY Array sizes are different NDH vs. NBAS
            PHIRES(1:M,j)=PHIJNEW(1:M)
          end do
        end do
! CMD transforming loops to matrix multiplications
!  H2 = Phi_occupied^T * FOHAM * Phi_unoccupied
!
!  H1 = Phi_occupied^T * FOHAM
!        call write_phires
        if(N > 0) then !Condition needed for 1e system
        call dgemm(transT,transn,N,M,M,alpha, PHIRES, M,
     &      FOHAM,matlen, beta, H1,N)
!  H2 = H1 * Phi_unoccupied
        call dgemm(transn,transn,N,N2,M,alpha, H1, N,
     &      PHIRES(1,N+1),M, beta, H2,N)
        end if

        temp=maxval(abs(H2))
!        CALL TRACER('ITER',ITER,TEMP)
! ORIGINAL CODE
!        temp=0.0d0
!        do i=1,N !occupied
!         do j=N+1,M !unoccupied
!          HAMTEMP=0.0d0
!          do k=1,M
!           do l=1,M
!            HAMTEMP=HAMTEMP+PHIRES(k,i)*PHIRES(l,j)*FOHAM(k,l)
!           end do
!          end do
!          if(abs(HAMTEMP)>temp) temp=abs(HAMTEMP)
!         end do
!        end do
        if(temp<CONVG) then
         PHIRES=transpose(PHIRES)
         CALL TRACER('ITER=',ITER)
         call gttime(time4)
         call timout('DIAGV DIAGONALIZATION',time4-time3)
         exit
        end if
       end do

       if(temp>CONVG) then
        print*,'WARNING: DIAGV not converged! ',temp
        PHIRES=transpose(PHIRES)
       end if

       END SUBROUTINE DIAGV_FO
C *********************************************************************
C
       SUBROUTINE DIAGV_FO_UNI(matlen,N,M,FOHAM,PHIRES,IPREONLY)
! Converion to implicit none.  Raja Zope Thu Aug 17 14:34:41 MDT 2017
! Diagonalize unified Hamiltonian in the local orbital basis
!

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'
       INTEGER :: IPREONLY, ITER, MAXITER
       REAL*8 :: SYMBOL , FOHAM, PHIRES, CONVG, FOEIG, FOHAMAVG, FOOVR,
     & FOSC, HAMTEMP, HII, HIJ, HJJ, PHIINEW, PHIJNEW, T, TEMP
!       INCLUDE 'commons.inc'

       INTEGER, intent(in) :: N
       INTEGER, intent(in) :: M
       INTEGER, intent(in) :: matlen

       INTEGER I,J,K,L
       LOGICAL PRINTHAM

       DIMENSION FOHAM(matlen,matlen), PHIRES(matlen,matlen)

       DIMENSION FOHAMAVG(M,M), FOOVR(M,M), FOEIG(M), FOSC(M)
       DIMENSION PHIINEW(M), PHIJNEW(M)

       PARAMETER (MAXITER=100)
       PARAMETER (CONVG=1.0d-6)
       SAVE

       PRINTHAM = .FALSE.

C      1.Generate OVERLAP Matrix
       FOOVR=0.0d0
       do i=1,M
        FOOVR(i,i)=1.0d0
        do j = 1,M
           FOHAMAVG(j,i) = FOHAM(j,i)
        end do
       end do

c
c  print unified HAM
c
       if(PRINTHAM) then 
        print *, 'in diagge_uni M N matlen', M, N, matlen
        do i = 1, matlen
          print 26, (FOHAMAVG(j,i), j = 1,matlen)
        end do
       end if
26     format(1x, 10(f8.3,1x))
       call DIAGGE(M,M,FOHAMAVG,FOOVR,FOEIG,FOSC,1)
       if(PRINTHAM) then
        print *, 'UNIHAM evalues'
        print 26, (FOEIG(j),j=1,matlen)
c
        print *, 'AFTER UNIHAM DIAG'
c
        do i = 1, matlen
          print 26, (FOHAMAVG(j,i), j = 1,matlen)
        end do
       end if

c
c  evectors should be placed into hamiltonian matrix
c
       PHIRES=0.0d0
       do i=1,M
        do j=1,M
         PHIRES(i,j)=FOHAMAVG(j,i)
        end do
       end do

       END SUBROUTINE DIAGV_FO_UNI
