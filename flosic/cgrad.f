C UTEP Electronic Structure Lab (2020)
c
c ************************************************************
c
c cgrad (Dirk Porezag, November 1998)
c performs a conjugate-gradient relaxation
c uses units 43 and 45 to write files
c
c input:  nopt:   number of degrees of freedom
c         fuval:  function value for current structure
c         xvec:   coordinates for current structure
c         gvec:   gradient for current structure
c         gtol:   convergence margin for gradient
c         ftol:   accuracy of fuval (if two function values differ by 
c                 less than fuval, they are considered identical
c         scrv:   scratch vector of size >= 6*nopt
c output: xvec:   new coordinates 
c         istat:  status: 0: converged
c                         1: interval expansion
c                         2: quadratic interpolation
c                         3: linear interpolation
c                         4: bisection
c
       SUBROUTINE CGRAD(NOPT,FUVAL,XVEC,GVEC,GTOL,FTOL,SCRV,ISTAT)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION XVEC(NOPT),GVEC(NOPT),SCRV(NOPT,6)
       PARAMETER (MAXLINE=30)
       PARAMETER (MXOPTIM=5)
       LOGICAL   IRESET,LINPOS,LINNEG
       CHARACTER ICHRA,ICHRB
       DIMENSION GAMMA(MAXLINE+1),FUNCT(MAXLINE),DERIV(MAXLINE)
       SAVE
       DATA EPS  /1.0D-6/
c
c setup. gltol is the convergence margin for the line minimization
c
       ISTAT= 0
       GLTOL= 0.5D0*ABS(GTOL)
       IF (NOPT .LT. 1) GOTO 900
c
c check if gradients have converged
c
       GMAX= 0.0D0
       DO IOPT= 1,NOPT
        GMAX= MAX(GMAX,ABS(GVEC(IOPT)))
       END DO
       IF (GMAX .LE. ABS(GTOL)) GOTO 900
c
c check if file cgrad is okay. if not, start new optimization
c if first calculation, start also new optimization
c
       IRESET= .TRUE.
       ISTEP= 0
       ILINE= 0
       ILNOPT= 0
       DELTA= 0.0D0
       NOPTRD= 0
       OPEN(43,FILE='CGRAD',FORM='unformatted',STATUS='unknown')
       REWIND(43)
       READ(43,END=10) ISTEP,ILINE,ILNOPT,NOPTRD,DELTA
       IF (NOPTRD .EQ. NOPT) IRESET=.FALSE.
   10  CONTINUE
c
c read in small, dmag from cgrlog
c
       SMALL= 0.1D0
       DMAG=  2.0D0
       OPEN(45,FILE='CGRLOG',FORM='formatted',STATUS='unknown') 
       REWIND(45)
       READ(45,*,END=20) SMALL,DMAG
   20  CONTINUE 
       IF (IRESET) THEN
        REWIND(45)
        WRITE(45,*) SMALL,DMAG
       END IF 
       CLOSE(45)
c
c set up values for conjugate gradient minimization
c if iline is larger than zero, we are in a line minimization
c
       REWIND(43)
       IF (IRESET) THEN
        ISTEP= 0
        ILINE= 0
        ILNOPT= 0
        DELTA= SMALL
       ELSE
        IF (ILINE .GE. 1) GOTO 100
        READ(43) ISTEP,ILINE,ILNOPT,NOPTRD,DELTA
        DO I= 1,5
         READ(43)(SCRV(IOPT,I),IOPT= 1,NOPT)
        END DO
       END IF
c
c start of conjugate-gradient method
c
   40  CONTINUE
       ISTEP= ISTEP+1
       ILINE= 0
       ILNOPT= 0
       IF (ISTEP .GT. NOPT) ISTEP= 1
c
c first step: start with negative gradient as search direction
c
       IF (ISTEP .EQ. 1) THEN
        DO IOPT= 1,NOPT
         SCRV(IOPT,5)= -GVEC(IOPT)
        END DO  
       ELSE
c
c construction of new hcgr= scrv(*,5) using polak-ribiere formula
c hcgr is the conjugate direction of travel
c ucgr= scrv(*,6) is the unit vector corresponding to hcgr
c
        GAM= 0.0D0
        GDIV= 0.0D0
        DO IOPT= 1,NOPT
         GAM= GAM+(GVEC(IOPT)-SCRV(IOPT,3))*GVEC(IOPT)
         GDIV= GDIV+SCRV(IOPT,3)**2
        END DO
        GAM= GAM/GDIV 
        DO IOPT= 1,NOPT
         SCRV(IOPT,5)= -GVEC(IOPT)+GAM*SCRV(IOPT,5)
        END DO
       END IF
       HNRM= 0.0D0
       DO IOPT= 1,NOPT
        HNRM= HNRM+SCRV(IOPT,5)**2
       END DO
       HNRM= 1.0D0/SQRT(HNRM)
       DO IOPT= 1,NOPT
        SCRV(IOPT,6)= SCRV(IOPT,5)*HNRM
       END DO  
c
c assign initial best values for line minimization
c write data to file cgrad
c
       BEST= FUVAL
       DO IOPT= 1,NOPT
        SCRV(IOPT,1)= XVEC(IOPT)
        SCRV(IOPT,2)= XVEC(IOPT)
        SCRV(IOPT,3)= GVEC(IOPT)
        SCRV(IOPT,4)= GVEC(IOPT)
       END DO  
       GAMMA(1)= 0.0D0
       FUNCT(1)= 0.0D0
       DERIV(1)= 0.0D0
       REWIND(43)
       WRITE(43) ISTEP,ILINE,ILNOPT,NOPT,DELTA
       DO I= 1,6
        WRITE(43)(SCRV(IOPT,I),IOPT= 1,NOPT)
       END DO
       WRITE(43)(GAMMA(I),I=1,ILINE+1)
       WRITE(43)(FUNCT(I),I=1,ILINE)
       WRITE(43)(DERIV(I),I=1,ILINE)
       WRITE(43) BEST
c
c begin/continue line minimiztion in direction of u
c first, read in data
c
  100  CONTINUE
       IWARN= 0
       IMODE= 0
       ICHRA= ' '
       ICHRB= ' '
       REWIND(43)
       READ(43) ISTEP,ILINE,ILNOPT,NOPTRD,DELTA
       DO I= 1,6
        READ(43)(SCRV(IOPT,I),IOPT= 1,NOPT)
       END DO
       READ(43)(GAMMA(I),I=1,ILINE+1)
       READ(43)(FUNCT(I),I=1,ILINE)
       READ(43)(DERIV(I),I=1,ILINE)
       READ(43) BEST
       ILINE= ILINE+1
c
c get funct, deriv, and umax for current point
c
       FUNCT(ILINE)= FUVAL
       DERIV(ILINE)= 0.0D0
       UMAX= 0.0D0
       DO IOPT= 1,NOPT
        DERIV(ILINE)= DERIV(ILINE)+SCRV(IOPT,6)*GVEC(IOPT)
        UMAX= MAX(UMAX,ABS(SCRV(IOPT,6)))
       END DO  
c
c update best values 
c
       IF (FUVAL .LE. BEST) THEN
        BEST= FUVAL
        DO IOPT= 1,NOPT
         SCRV(IOPT,2)= XVEC(IOPT)
         SCRV(IOPT,4)= GVEC(IOPT)
        END DO 
       END IF
c
c check for convergence of line minimization
c check if converged forces have lead to the best function value
c if not, print warning and take best value available
c
       IF (ABS(DERIV(ILINE)) .LT. GLTOL) THEN
        IF (FUVAL-BEST .GT. ABS(FTOL)) THEN
         IWARN= 1
         CALL LOGCGR(IWARN,IMODE,ICHRA,ICHRB,ILINE,ISTEP)
         DO IOPT= 1,NOPT
          XVEC(IOPT)= SCRV(IOPT,2)
          GVEC(IOPT)= SCRV(IOPT,4)
         END DO
         FUVAL= BEST
        END IF 
        GOTO 40
       END IF
c
c calculate next value for gamma
c first step: go finite distance into the search direction
c
       IF (ILINE.EQ.1) THEN
        IMODE= 1
        GAMMA(1)= 0.0D0
        IF (DERIV(1) .GT. 0.0D0) THEN
         GAMMA(2)= -DELTA/UMAX
        ELSE
         GAMMA(2)=  DELTA/UMAX
        END IF
       ELSE 
c
c check whether brackets already available
c 
        LINPOS= .FALSE.
        LINNEG= .FALSE.
        DO JLINE= 1,ILINE
         IF (DERIV(JLINE) .GT. 0.0D0) THEN
          LINPOS= .TRUE.
         ELSE
          LINNEG= .TRUE.
         END IF
        END DO
c
c if no brackets available, magnify interval
c
        IF (.NOT. (LINPOS .AND. LINNEG)) THEN
         IMODE= 1
         GAMMA(ILINE+1)= DMAG*GAMMA(ILINE)+GAMMA(2)
        ELSE
c
c do better update, for (iline .eq. 2) try linear interpolation 
c
         ILNOPT= ILNOPT+1
         IF (ILINE .EQ. 2) THEN
          GAMMA(3)= 0.5D0*GAMMA(2)
          IF (ABS(DERIV(1)-DERIV(2)) .GT. EPS) THEN
           IMODE= 3
           GAMMA(3)= GAMMA(1)*DERIV(2)-GAMMA(2)*DERIV(1)
           GAMMA(3)= GAMMA(3)/(DERIV(2)-DERIV(1))
          ELSE
           IMODE= 4
           GAMMA(3)= 0.5D0*GAMMA(2)
          END IF
         ELSE 
c
c we have brackets and at least three points
c we need to find the gamma value for which deriv is zero
c first try quadratic interpolation, then linear interpolation,
c and if this fails too fall back to bisection        
c sort gamma, funct and deriv (best value -> index 0)
c          
          DO JLINE= 1,ILINE
           DO KLINE= JLINE+1,ILINE
            IF (FUNCT(KLINE) .LT. FUNCT(JLINE)) THEN
             CALL SWAP(GAMMA(KLINE),GAMMA(JLINE))
             CALL SWAP(FUNCT(KLINE),FUNCT(JLINE))
             CALL SWAP(DERIV(KLINE),DERIV(JLINE))
            END IF
           END DO
          END DO
c
c get indices of the three interesting points
c if deriv(1) and deriv(2) have different signs, 
c the third point has index 3
c          
          IF (DERIV(1)*DERIV(2) .GT. 0.0D0) THEN
           DO JLINE= 3,ILINE
            IF (DERIV(JLINE)*DERIV(1) .LT. 0.0D0) GOTO 105
           END DO
  105      CONTINUE
           CALL SWAP(GAMMA(3),GAMMA(JLINE))
           CALL SWAP(FUNCT(3),FUNCT(JLINE))
           CALL SWAP(DERIV(3),DERIV(JLINE))
          END IF
c
c sort the first three values so that gamma(1), gamma(2)
c and gamma(3) are in increasing order 
c
          DO JLINE= 1,2
           DO KLINE= JLINE+1,3
            IF (GAMMA(KLINE) .LT. GAMMA(JLINE)) THEN
             CALL SWAP(GAMMA(KLINE),GAMMA(JLINE))
             CALL SWAP(FUNCT(KLINE),FUNCT(JLINE))
             CALL SWAP(DERIV(KLINE),DERIV(JLINE))
            END IF
           END DO
          END DO
c
c get bracketing triple a,b,c
c
          AG= GAMMA(1)
          BG= GAMMA(2)
          CG= GAMMA(3)
          AF= FUNCT(1)
          CF= FUNCT(3)
          A1= DERIV(1)
          B1= DERIV(2)
          C1= DERIV(3)
c
c check for strange behavior of function points
c if (ag .gt. 0) and (cg .lt. 0), expand interval in direction
c of smallest function value
c
          IF ((A1 .GT. 0.0D0) .AND. (C1 .LT. 0.0D0)) THEN
           IWARN= 2
           IMODE= 1
           ICHRA= 'a' 
           ICHRB= 'c' 
           IF (AF .LT. CF) THEN
            GNEW= AG-DELTA/UMAX
           ELSE
            GNEW= CG+DELTA/UMAX
           END IF
           GOTO 200
          END IF
c
c if either one of the left or right brackets has the wrong sign,
c there is a maximum in between -> try linear interpolation 
c between the two points having different signs
          
          IF ((A1 .GT. 0.0D0) .AND. (C1 .GT. 0.0D0)) THEN
           IWARN= 2
           ICHRA= 'a'
           ICHRB= 'b'
           AG= BG
           BG= CG
           A1= B1
           B1= C1
           GOTO 120
          END IF
          IF ((A1 .LT. 0.0D0) .AND. (C1. LT. 0.0D0)) THEN
           IWARN= 2
           ICHRA= 'b'
           ICHRB= 'c'
           GOTO 120
          END IF
c
c try quadratic interpolation
c
          G2= BG-AG
          G3= CG-AG
          D2= B1-A1
          D3= C1-A1
          DN= (G3-G2)*G3*G2
          IF (ABS(DN) .LE. 1D-20) GOTO 120
          DN= 1.0D0/DN
          APOL= (D2*G3*G3-D3*G2*G2)*DN
          BPOL= (D3*G2-D2*G3)*DN
          IF (ABS(BPOL).LT.1D-10) GOTO 120
          BPOLR= 1.0D0/BPOL
          DHLP= -0.5D0*APOL*BPOLR
          DRT= DHLP*DHLP-A1*BPOLR
          IF (DRT.LT.1D-20) GOTO 120
          DRT= SQRT(DRT)
          GNEW= DHLP+DRT
          IF ((2*BPOL*GNEW+APOL) .LT. 0.0D0) GNEW= DHLP-DRT
          GNEW= GNEW+AG
          IF ((GNEW .LT. AG) .OR. (GNEW .GT. CG)) GOTO 120 
          IMODE= 2
          GOTO 200
c
c try linear interpolation
c
  120     IF (B1 .LT. 0.0D0) THEN
           AG= BG
           BG= CG
           A1= B1
           B1= C1
          END IF
          IF ((B1-A1) .LE. 1D-10) GOTO 140
          GNEW= (B1*AG-A1*BG)/(B1-A1)
          IF ((GNEW.LT.AG) .OR. (GNEW.GT.CG)) GOTO 140 
          IMODE= 3
          GOTO 200     
c
c bisection
c
  140     IMODE= 4
          GNEW= 0.5D0*(AG+BG)   
  200     GAMMA(ILINE+1)= GNEW
         END IF
        END IF
       END IF
c
c write results in logfile
c
       CALL LOGCGR(IWARN,IMODE,ICHRA,ICHRB,ILINE,ISTEP)
c
c if maxline or mxoptim exceeded, return best value for xvec and 
c start new line minimization
c
       IF ((ILINE .GE. MAXLINE) .OR. (ILNOPT .GE. MXOPTIM)) THEN
        DO IOPT= 1,NOPT
         XVEC(IOPT)= SCRV(IOPT,2)
         GVEC(IOPT)= SCRV(IOPT,4)
        END DO
        FUVAL=BEST
        GOTO 40
       END IF
c
c update xvec
c
       DO IOPT= 1,NOPT
        XVEC(IOPT)= SCRV(IOPT,1)+GAMMA(ILINE+1)*SCRV(IOPT,6)
       END DO
c
c restore history for line minimization
c
       REWIND(43)
       WRITE(43) ISTEP,ILINE,ILNOPT,NOPTRD,DELTA
       DO I= 1,6
        WRITE(43)(SCRV(IOPT,I),IOPT= 1,NOPT)
       END DO
       WRITE(43)(GAMMA(I),I=1,ILINE+1)
       WRITE(43)(FUNCT(I),I=1,ILINE)
       WRITE(43)(DERIV(I),I=1,ILINE)
       WRITE(43) BEST
       CLOSE(43)
       ISTAT= IMODE
c
c if (istat .eq. 0) remove file cgrad
c
  900  IF (ISTAT .EQ. 0) THEN
        OPEN(43,FILE='CGRAD',FORM='unformatted',STATUS='unknown')
        CLOSE(43,STATUS='delete')
       END IF
       RETURN
       END
