C UTEP Electronic Structure Lab (2020)
         SUBROUTINE OVLP1D(ALP,BET,A,B,NFMX,OV)
         IMPLICIT REAL*8 (A-H,O-Z)
         PARAMETER (MPMX=17)
         DATA NFMXS/-1/
         DIMENSION POW(0:MPMX,2),TT(0:MPMX)
         DIMENSION BIN(0:MPMX,0:MPMX),FCTRL(0:MPMX)
         DIMENSION OV(0:NFMX,0:NFMX)
         SAVE
c integrand = (x-a)^n exp[-alp*(x-a)**2] (x-b)^m exp[-alp*(x-b)**2]   
c (n,m)=0,...,nfmx
         IF(NFMXS.LT.NFMX)THEN
           PI=4.0D0*ATAN(1.0D0)
           NFMXS=NFMX
           IF(NFMX.GT.MPMX)STOP'OVLP1D'
           FCTRL(0)=1.0D0
               DO I=1,NFMX
               FCTRL(I)=FCTRL(I-1)*I
               END DO
C          write(6,*)(FCTRL(J),J=0,NFMX)
c calculate binomial coeficients:
         DO N=0,NFMX
         DO M=0,NFMX
         BIN(N,M)=0.0D0                          
         END DO
         END DO
         DO N=0,NFMX
         DO M=0,N      
         BIN(N,M)=FCTRL(N)/FCTRL(M)/FCTRL(N-M)
         END DO
         END DO
         END IF
         RECPR=1.0D0/(ALP+BET)
         C=(ALP*A+BET*B)*RECPR         
         ARG=ALP*BET*(A-B)*(A-B)*RECPR    
         ARG=EXP(-ARG)
         TT(0)=SQRT(PI*RECPR)
c           write(6,*)'tt:',tt(0),pi*recpr
         RECPR=RECPR*0.5
           DO M=1,NFMX
           TT(M)=(2*M-1)*TT(M-1)*RECPR
           END DO
C           write(6,*)(TT(I),I=0,M)
           DO M=0,NFMX
           TT(M)=TT(M)*ARG
           END DO
c           write(6,*)(tt(i),i=0,m)
               CMA=C-A
               CMB=C-B
               POW(0,1)=1.0D0
               DO M=1,NFMX
               POW(M,1)=POW(M-1,1)*CMA
               END DO
               POW(0,2)=1.0D0
               DO M=1,NFMX
               POW(M,2)=POW(M-1,2)*CMB
               END DO
         DO M=0,NFMX
         DO N=0,NFMX
         OV(M,N)=0.0D0
             SS=0.0D0
             DO MP=0,M
               IF(2*(MP/2).EQ.MP)THEN
                  NPB=0
               ELSE
                  NPB=1
               END IF
             DO NP=NPB,N,2
               COEF=BIN(M,MP)*BIN(N,NP)*POW(M-MP,1)*POW(N-NP,2)
               SS=SS+COEF*TT((MP+NP)/2) 
             END DO
             END DO
           OV(M,N)=SS
         END DO
         END DO
         RETURN
         END 
