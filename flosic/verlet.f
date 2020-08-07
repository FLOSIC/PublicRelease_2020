C UTEP Electronic Structure Lab (2020)
C
C *****************************************************************
C
       SUBROUTINE VERLET(NPAR,X,G,GAMMA)
       use common3,only : RMAT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NPAR, I, MPAR, MXPAR
       REAL*8 :: SYMBOL , X, G, GAMMA, DRAG1, DRAG2, DT, TNEW, TOLD,
     & X0, X1, XOLD
       SAVE
       PARAMETER (MXPAR=3*MAX_IDENT)
       LOGICAL EXIST,JEREMIES_WAY
       DIMENSION X(NPAR),G(NPAR),XOLD(MXPAR)
       DATA JEREMIES_WAY/.TRUE./
       DATA DT/0.03/
       INQUIRE(FILE='VERLET',EXIST=EXIST)
       OPEN(60,FILE='VERLET',STATUS='UNKNOWN',FORM='FORMATTED')
       IF(NPAR.GT.MXPAR)THEN
        write(6,*)'VERLET: MXPAR MUST BE AT LEAST: ',NPAR
        GOTO 900
       END IF
       IF(DT*GAMMA .GT. 1.0D0)THEN
        write(6,*)'VERLET: GAMMA*DT IS TOO BIG'
        GOTO 900
       END IF 
       IF(EXIST) THEN
        READ(60,*) MPAR
        IF(MPAR.NE.NPAR) THEN
         PRINT *,'MPAR AND NPAR ARE NOT EQUAL IN VERLET'
         GOTO 900
        END IF
        READ(60,*) (XOLD(I), I=1,MPAR)
        READ(60,*) TOLD
        REWIND(60)
        DRAG1 = 1.0D0 - 0.5D0*GAMMA*DT
        DRAG2 = 1.0D0 + 0.5D0*GAMMA*DT
        DRAG2=1.0D0/DRAG2
        TNEW=0.0D0
        DO 200 I=1,MPAR
         X0=XOLD(I)
         X1 = 2*X(I) - XOLD(I)*DRAG1 - G(I)*DT*DT
         XOLD(I)=X(I)
         X(I)   =X1*DRAG2
         TNEW=TNEW+0.5D0*( (X(I)-X0)/(2*DT) )**2
  200   CONTINUE
       ELSE
        TOLD=-1.0D0
        TNEW= 0.0D0
        MPAR=NPAR
        DO 400 I=1,MPAR
         XOLD(I)=X(I)
C
C INITIAL VELOCITIES (MORE OR LESS)
C
         X(I)=X(I) - (2*DT)*G(I)/(1.0D0+ABS(G(I)))
  400   CONTINUE
       END IF
       WRITE(60,*) MPAR
       WRITE(60,*) (XOLD(I),I=1,MPAR)
       WRITE(60,*)TNEW
       IF((TNEW.LT.TOLD).AND.JEREMIES_WAY)THEN
        CLOSE(60,STATUS='DELETE')
       ELSE
        CLOSE(60)
       END IF
       RETURN
  900  CLOSE(60)
       RETURN
       END
