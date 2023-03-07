            subroutine gint_new_test
            implicit real*8 (a-h,o-z)
            dimension w7(10,10,7),a3(3),b3(3)
C            print*,'alp,bet,a,b=?'
c           read*,alp,bet,A3,B3
            alp=0.03 
            bet=0.04
            a3(1)=0.0
            a3(2)=0.0
            a3(3)=0.0
            b3(1)=1.0
            b3(2)=1.0
            b3(3)=1.0
            call GINT_NEW(10,2,2,alp,bet,A3,B3,W7)
            end 

            subroutine GINT_NEW(mxL,np,nk,alp,bet,A3,B3,W7) 
            implicit real*8 (a-h,o-z)
            parameter (nd=8)
            character*3 coord
            character*60 line
            dimension oned(-1:nd,-1:nd,0:4,3)
            dimension savd(-1:nd,-1:nd)
            dimension a3(3),b3(3)
            dimension w7(mxL,mxL,7),lp(3,mxL),w0(10,10)
            dimension w4(10,10,4)
            data coord/'xyz'/
         
c 
c
c calculate int dx (x-A)^n exp(-alp*(x-a)^2) (x-B)^m exp(-alp*(x-B)^2  x^k 
c n = 0 to np
c m = 0 to np
c k = 0 to nk
c
c
c           oned       =0.0d0
c           w7         =0.0D0
            do 100 ix=1,3
            savd       =0.0d0
c           print*,'Starting Integrals for:',coord(ix:ix)
              a=a3(ix)
              b=b3(ix)
            arg=alp*bet*(A-B)*(A-B)/(alp+bet)
             f1=1.0D0/(alp+bet) 
             f2=f1*(A*alp+B*bet)
            fact=sqrt(4.0d0*atan(1.0D0)/(alp+bet))  !check this formulae
            c=(alp*A+bet*B)/(alp+bet)
            cmb= (c-b)
            cma= (c-a)
            fct=0.5d0/(alp+bet) 
            oned(0,0,0,ix)=fact*exp(-arg)
            mp=np+nk+2
                   if(mp.gt.nd)then
                   print*,'nd must be>',mp,' in gint_new'
                   call stopit
                   end if
            nway=1 ! change to 2 for debugging
            do iway=1,nway
            do n=0,mp
            do m=n,mp  
            oned(n,m+1,0,ix)=cmb*oned(n,m,0,ix)
     &                 +(n*fct)*oned(n-1,m,0,ix)
     &                 +(m*fct)*oned(n,m-1,0,ix)
            oned(m+1,n,0,ix)=cma*oned(m,n,0,ix)
     &                 +(n*fct)*oned(m,n-1,0,ix)
     &                 +(m*fct)*oned(m-1,n,0,ix)
            end do
c           print*,'n:',n
c           do l=0,mp
c           print 50, (oned(l,m,0,ix),m=0,mp)
c           end do
        
            alt3         =oned(n+2,n,0,ix)+(a-b)*oned(n+1,n,0,ix)
            alt4         =oned(n,n+2,0,ix)+(b-a)*oned(n,n+1,0,ix)
                oned(n+1,n+1,0,ix)=(alt3+alt4)/2.0D0
c           if(iway.eq.1)oned(n+1,n+1,0,ix)=alt3
c           if(iway.eq.2)oned(n+1,n+1,0,ix)=alt4
            end do
            if(nway.eq.2)then
            if(iway.eq.2)print*,'Error Check:'
            do n=0,mp
            print 50,(oned(n,m,0,ix)-savd(n,m),m=0,mp)
            end do
            end if
                      if(nway.eq.2.and.iway.eq.1)then
                      do n=0,mp
                      do m=0,mp
                      savd(n,m)=oned(n,m,0,ix)
                      end do
                      end do
                      end if
            end do
c Inmk = Int (x-a)^n        (x-b)^m x^k  exp[ ]
c      = Int (x-a)^n        (x-b)^m (x-a+a) x^(k-1)
c      = Int (x-a)^(n+1)    (x-b)^m x^(k-1)  + a Int (x-a)^n (x-b)^m x^(k-1)
            do k=1,nk
                do m=0,np+nk-k
                do n=0,np+nk-k
           oned(n,m,k,ix)=(oned(n+1,m,k-1,ix)+a*oned(n,m,k-1,ix)
     &                   + oned(n,m+1,k-1,ix)+b*oned(n,m,k-1,ix))/2.0D0
c
c I think the right strategy above is to use is to  us a if |a|<|b| an to use |b| if |b|<a.
c
c             oned(n,m,k,ix)= f1*(alp*oned(n+1,m,k-1,ix)
c    &                          +bet*oned(n,m+1,k-1,ix))
c    &                      +f2     *oned(n,m,k-1,ix)
                 end do
                end do
             end do
             if(nway.eq.2) then 
             do k=0, nk
      line=
     &'<(#-A#)^n exp(-alp(#-A#)^2)| #^9 |(#-B#)^n exp(-bet(#-B#)^2>' 
            do i=1,60
               if(line(i:i).eq.'#')line(i:i)=coord(ix:ix)
               if(line(i:i).eq.'9')then
                   write(line(i:i),'(I1)')k 
               end if
            end do
            print*, line
 15         format
     &('<(x-A)^n exp(-alp(x-A)^2| x^',i1,'|(x-B)^n exp(-bet(x-B)^2>')
              do  n=0,np
              print 50,(oned(n,m,k,ix),m=0,np)
              end do
            end do
            end if   ! NWAY 
 100        continue
c now assemble for d only case:
               lp(1,1)=0
               lp(2,1)=0
               lp(3,1)=0
c
               lp(1,2)=1
               lp(2,2)=0
               lp(3,2)=0
c
               lp(1,3)=0
               lp(2,3)=1
               lp(3,3)=0
c
               lp(1,4)=0
               lp(2,4)=0
               lp(3,4)=1
c
               lp(1,5)=2
               lp(2,5)=0
               lp(3,5)=0
c
               lp(1,6)=0
               lp(2,6)=2
               lp(3,6)=0
c
               lp(1,7)=0
               lp(2,7)=0
               lp(3,7)=2
c
               lp(1,8)=1
               lp(2,8)=1
               lp(3,8)=0
c
               lp(1,9)=1
               lp(2,9)=0
               lp(3,9)=1
c
               lp(1,10)=0
               lp(2,10)=1
               lp(3,10)=1
c
c 7x2x100=1400 multiplications
c could be reduced significantly
            do il=1,10
            do ir=1,10
c
            w7(il,ir,1)=oned(lp(1,il),lp(1,ir),0,1)
     &                 *oned(lp(2,il),lp(2,ir),0,2)
     &                 *oned(lp(3,il),lp(3,ir),0,3)
c
            w7(il,ir,2)=oned(lp(1,il),lp(1,ir),1,1)
     &                 *oned(lp(2,il),lp(2,ir),0,2)
     &                 *oned(lp(3,il),lp(3,ir),0,3)
c
            w7(il,ir,3)=oned(lp(1,il),lp(1,ir),0,1)
     &                 *oned(lp(2,il),lp(2,ir),1,2)
     &                 *oned(lp(3,il),lp(3,ir),0,3)
c
            w7(il,ir,4)=oned(lp(1,il),lp(1,ir),0,1)
     &                 *oned(lp(2,il),lp(2,ir),0,2)
     &                 *oned(lp(3,il),lp(3,ir),1,3)
c
            w7(il,ir,5)=oned(lp(1,il),lp(1,ir),2,1)
     &                 *oned(lp(2,il),lp(2,ir),0,2)
     &                 *oned(lp(3,il),lp(3,ir),0,3)
c
            w7(il,ir,6)=oned(lp(1,il),lp(1,ir),0,1)
     &                 *oned(lp(2,il),lp(2,ir),2,2)
     &                 *oned(lp(3,il),lp(3,ir),0,3)
c
            w7(il,ir,7)=oned(lp(1,il),lp(1,ir),0,1)
     &                 *oned(lp(2,il),lp(2,ir),0,2)
     &                 *oned(lp(3,il),lp(3,ir),2,3)
            end do
            end do
            if(nway.eq.2) then
            call svmxsf(alp,bet,A3,B3,w0)
            call gintee(alp,bet,A3,B3,w4)
            err=0.0d0
            do i=1,10
            do j=1,10
            err=err+abs(w0(j,i)-w7(j,i,1))
            end do
            print 50,(w7(j,i,1),j=1,10)
            print 50,(w0(j,i  ),j=1,10)
            print*,' '
            end do
            print*,'Error:',err
            do k=1,3
            print*,'Error for <',coord(k:k),'>'
            err=0.0d0
            do i=1,10
            do j=1,10
            err=err+abs(w4(j,i,k)-w7(j,i,k+1))
            end do
            print 50,(w7(j,i,k+1),j=1,10)
            print 50,(w4(j,i,k  ),j=1,10)
            print*,' '
            end do
            print*,'Total Error:',err
            end do
            end if !nway 
 
 50         format(' ',20g15.6)
            return
            end   
         
      SUBROUTINE GINTEE(A1,A2,A,B,W)
      IMPLICIT  REAL*8 (A-H,O-Z)
      DIMENSION A(3),B(3),W(10,10,4),L3(3,4)
      DIMENSION V(3,3,5,3)
      DIMENSION ML(10,3)
      SAVE
      DATA  L3/2,1,1,  1,2,1,  1,1,2,  1,1,1/
      DATA ML/1,2,1,1,3,1,1,2,2,1,
     1        1,1,2,1,1,3,1,2,1,2,
     2        1,1,1,2,1,1,3,1,2,2/
      AT=1.0D0/(A1+A2)
      Q1=0.5D0*AT
      Q2=Q1*Q1
      XS=0.0D0
      DO 10 IFT=1,3
      XS=XS+(A1*A2*((A(IFT)-B(IFT))**2))*AT
      D=(A1*A(IFT)+A2*B(IFT))*AT
      F1=D-A(IFT)
      F2=D-B(IFT)
      F3=D
      F001=F3
      F010=F2
      F011=F010*F3
      F020=F010*F2
      F021=F020*F3
      F100=F1
      F101=F100*F3
      F110=F100*F2
      F111=F110*F3
      F120=F110*F2
      F121=F120*F3
      F200=F100*F1
      F201=F200*F3
      F210=F200*F2
      F211=F210*F3
      F220=F210*F2
      F221=F220*F3
      V(1,1,1,IFT)=1
      V(1,2,1,IFT)=+F010
      V(1,3,1,IFT)=+F020+Q1
      V(2,1,1,IFT)=+F100
      V(2,2,1,IFT)=+F110+Q1
      V(2,3,1,IFT)=+F120+Q1*F100+2*Q1*F010
      V(3,1,1,IFT)=+F200+Q1
      V(3,2,1,IFT)=+F210+2*Q1*F100+Q1*F010
      V(3,3,1,IFT)=+F220+Q1*F200+4*Q1*F110+Q1*F020+3*Q2
      V(1,1,2,IFT)=+F001
      V(1,2,2,IFT)=+F011+Q1
      V(1,3,2,IFT)=+F021+2*Q1*F010+Q1*F001
      V(2,1,2,IFT)=+F101+Q1
      V(2,2,2,IFT)=+F111+Q1*F100+Q1*F010+Q1*F001
      V(2,3,2,IFT)=+F121+2*Q1*F110+Q1*F101+Q1*F020+2*Q1*F011+3*Q2
      V(3,1,2,IFT)=+F201+2*Q1*F100+Q1*F001
      V(3,2,2,IFT)=+F211+Q1*F200+2*Q1*F110+2*Q1*F101+Q1*F011+3*Q2
      V(3,3,2,IFT)=+F221+2*Q1*F210+Q1*F201+2*Q1*F120+4*Q1*F111
     &             +6*Q2*F100+Q1*F021+6*Q2*F010+3*Q2*F001
 10   CONTINUE
      D=EXP(-XS)*(SQRT(3.14159265358979324D0*AT)**3)
      DO 40 K=1,4 
       DO 35 J=1,10
        DO 30 I=1,10
         W(I,J,K)=V(ML(I,1),ML(J,1),L3(1,K),1)*
     &            V(ML(I,2),ML(J,2),L3(2,K),2)*
     &            V(ML(I,3),ML(J,3),L3(3,K),3)
   30   CONTINUE
   35  CONTINUE
   40 CONTINUE
      DO 50 K=1,4
       DO 48 J=1,10
        DO 46 I=1,10
         W(I,J,K)=W(I,J,K)*D
  46    CONTINUE
  48   CONTINUE
  50  CONTINUE
      END


      SUBROUTINE SVMXSF(A1,A2,R,T,WO)
      IMPLICIT  REAL*8 (A-H,O-Z)
      PARAMETER (ND=10)
      LOGICAL DYES
      DIMENSION WO(ND,ND),R(3),T(3),X(3)
      SAVE
      DATA DYES/.FALSE./
      X(1)=T(1)-R(1)
      X(2)=T(2)-R(2)
      X(3)=T(3)-R(3)
      A=A1*A2/(A1+A2)
      AS=A*A
      PI=3.14159265358979324D0
      B=PI/(A1+A2)
      B=SQRT(B)
      B=B**3
      DO 100 I=1,ND
      DO 100 J=1,ND
  100 WO(I,J)=0.0D0
      XS=X(1)**2+X(2)**2+X(3)**2
      IF( XS.LT.1.0D-14 ) GO TO 150
      C=A*XS
      C=EXP(-C)
C THIS SS
      WO(1,1)=B*C
      DO 105 I=1,3
C THIS IS SPX
      WO(1,I+1)=-A*B*C*X(I)/A2
C THIS IS PXS
      WO(I+1,1)=A*B*C*X(I)/A1
      IF(DYES) GO TO 1000
C THIS IS SDXX
      WO(1,I+4)=(0.5D0*B*C/A2)*((1.0D0-A/A2)+(2*A*A*X(I)*X(I)/A2))
C THIS IS DXXS
      WO(I+4,1)=(0.5D0*B*C/A1)*((1.0D0-A/A1)+(2*A*A*X(I)*X(I)/A1))
 1000 CONTINUE
C THIS IS PXPX
      WO(I+1,I+1)=A*B*C*(0.5D0-A*X(I)*X(I))/(A1*A2)
      IF(DYES) GO TO 1001
C THIS IS DXXPX
      WO(I+4,I+1)=-A*B*C*X(I)*((1.0D0-3*A/A1)+(2*A*A*X(I)*X(I)/A1))
     &           /(2*A1*A2)
C THIS IS PXDXX
      WO(I+1,I+4)= A*B*C*X(I)*((1.0D0-3*A/A2)+(2*A*A*X(I)*X(I)/A2))
     &           /(2*A1*A2)
C THIS IS DXXDXX
      WOW=3.0D0+X(I)**2*(2*(A1+A2)-12*A)+4*AS*X(I)**4
      WO(I+4,I+4)=WOW*B*C*AS/(4*A1**2*A2**2)
 1001 CONTINUE
  105 CONTINUE
      IF(DYES) GO TO 1002
      K=7
      DO 110 I=1,2
      DO 110 J=2,3
      IF (I.EQ.J) GO TO 110
      K=K+1
C THIS IS SDXY
      WO(1,K)=A*A*B*C*X(I)*X(J)/(A2*A2)
C THIS IS DXYS
      WO(K,1)=A*A*B*C*X(I)*X(J)/(A1*A1)
  110 CONTINUE
 1002 CONTINUE
      DO 115 I=1,3
      DO 115 J=1,3
      IF (I.EQ.J) GO TO 115
C THIS IS PXPY
      WO(I+1,J+1)=-A*A*B*C*X(I)*X(J)/(A1*A2)
      IF(DYES) GO TO 1003
C THIS IS PXDYY
      WO(I+1,J+4)=A*B*C*X(I)*((1.0D0-A/A2)+(2*A*A*X(J)*X(J)/A2))
     &           /(2*A1*A2)
C THIS IS DXXPY
      WO(I+4,J+1)=-A*B*C*X(J)*((1.0D0-A/A1)+(2*A*A*X(I)*X(I)/A1))
     &           /(2*A1*A2)
C THIS IS DXXDYY
      WOW=1.0D0+2*(A1-A)*X(J)**2+2*(A2-A)*X(I)**2+4*AS*X(I)**2*X(J)**2
      WO(I+4,J+4)=WOW*B*C*AS/(4*A1**2*A2**2)
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
      WO(I+1,J+K+5)=-A*A*B*C*(0.5D0-A*X(I)*X(I))*X(L)/(A1*A2*A2)
C THIS IS DXXDXY
      WO(I+4,J+K+5)=A*A*B*C*X(J)*X(K)*((1.0D0-3*A/A1)
     &             +2*A*A*X(I)*X(I)/A1)/(2*A1*A2*A2)
C THIS IS DXYPX
      WO(J+K+5,I+1)=A*A*B*C*(0.5D0-A*X(I)*X(I))*X(L)/(A1*A1*A2)
C THIS IS DXYDXX
      WO(J+K+5,I+4)=A*A*B*C*X(J)*X(K)*((1.0D0-3*A/A2)
     &             +2*A*A*X(I)*X(I)/A2)/(2*A1*A1*A2)
      GO TO 130
  120 CONTINUE
C THIS IS PXDYZ
      WO(I+1,J+K+5)=(A**3)*B*C*X(I)*X(J)*X(K)/(A1*A2*A2)
C THIS IS DXXDYZ
      WO(I+4,J+K+5)=A*A*B*C*X(J)*X(K)*((1.0D0-A/A1)
     &             +2*A*A*X(I)*X(I)/A1)/(2*A1*A2*A2)
C THIS IS DXYPZ
      WO(J+K+5,I+1)=-(A**3)*B*C*X(I)*X(J)*X(K)/(A1*A1*A2)
C THIS IS DXYDZZ
      WO(J+K+5,I+4)=A*A*B*C*X(J)*X(K)*((1.0D0-A/A2)
     &             +2*A*A*X(I)*X(I)/A2)/(2*A1*A1*A2)
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
  135 WO(I+J+5,K+L+5)=-(A**3)*B*C*X(M2)*X(M3)
     &               *(1.0D0-2*A*X(M1)*X(M1))/(2*A1*A1*A2*A2)
      GO TO 145
C THIS IS DXYDXY
  140 WO(I+J+5,K+L+5)=A*A*B*C*(1.0D0-2*A*X(I)*X(I)-2*A*X(J)*X(J)
     &               +4*A*A*X(I)*X(I)*X(J)*X(J))/(4*A1*A1*A2*A2)
  145 CONTINUE
  146 CONTINUE
 1004 CONTINUE
      GO TO 200
  150 CONTINUE
C THIS SS
      WO(1,1)=B
      DO 155 I=1,3
      IF(DYES) GO TO 1005
C THIS IS SDXX
      WO(1,I+4)=B*(1.0D0-A/A2)/(2*A2)
C THIS IS DXXS
      WO(I+4,1)=B*(1.0D0-A/A1)/(2*A1)
 1005 CONTINUE
C THIS IS PXPX
      WO(I+1,I+1)=0.5D0*A*B/(A1*A2)
      IF(DYES) GO TO 1006
C THIS IS DXXDXX
      WO(I+4,I+4)=3*A**2*B/(4*A1**2*A2**2)
 1006 CONTINUE
  155 CONTINUE
      IF(DYES) GO TO 1007
      DO 165 I=1,3
      DO 165 J=1,3
      IF (I.EQ.J) GO TO 165
C THIS IS DXXDYY
      WO(I+4,J+4)=A*A*B/(4*A1**2*A2**2)
C THIS IS DXYDXY
      WO(I+J+5,I+J+5)=WO(I+4,J+4)
  165 CONTINUE
 1007 CONTINUE
  200 CONTINUE
      RETURN
      END
C
