C UTEP Electronic Structure Lab (2020)
C
         SUBROUTINE FONEDM(ALPI,ALPJ,AI,AJ,NFMX,NDMX,SS)
c m.r. pederson 21-june 2000
         IMPLICIT REAL*8 (A-H,O-Z)
         PARAMETER (MXDIM=17)
         DIMENSION SS(0:NFMX,0:NFMX,0:NDMX)
         DIMENSION SP(0:MXDIM,0:MXDIM,0:NDMX)
         DIMENSION OV(0:MXDIM,0:MXDIM)
             IF(NFMX+NDMX.GT.MXDIM)THEN
             write(6,*)'PARAMETER MXDIM IN FONEDM TOO SMALL'
             write(6,*)'MUST BE:',NFMX+NDMX
             STOP
             END IF
c calculate all interesting integrals:
c 
c integrand = (x-ai)^n exp[-alpi*(x-ai)**2] 
c             d^l/dx^l  (x-aj)^m exp[-alpj*(x-aj)**2]   
c (n,m)=0,...,nfmx
c    l =0,...,ndmx
c use recursion relations...
             CALL OVLP1D(ALPI,ALPJ,AI,AJ,MXDIM,SP(0,0,0))
                   BET=ALPJ*2.0D0
             DO ID=1,NDMX
              DO N=0,NFMX
              SP(N,0,ID)=                -BET*SP(N,  1,ID-1)
              DO M=1,NFMX+NDMX-1 
              SP(N,M,ID)=M*SP(N,M-1,ID-1)-BET*SP(N,M+1,ID-1)
              END DO
              END DO
             END DO 
             DO ID=0,NDMX
             DO N =0,NFMX
             DO M =0,NFMX
              SS(N,M,ID)=SP(N,M,ID)
             END DO
             END DO
             END DO
         RETURN
         END
