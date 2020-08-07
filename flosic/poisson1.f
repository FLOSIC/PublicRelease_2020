C UTEP Electronic Structure Lab (2020)
C
C DIRECTIONS FOR USE:
C
C IF YOU ARE SUSPICIOUS OF A BUG:
C (1)Check if Answer is converged as a function of EPSILON in POISSON1
C (2)Check if Answer is converged as a function of THIRTYSOMETHING
C (3)Compare to exact answer that can be obtained by setting
C
C     EPSILON=0.0 THIRTYSOMETHING=1.0D30
C
C IF YOU CONFIRM THE EXISTENCE OF A BUG, PLEASE FIND IT :-)
C AND NOTIFY MARK PEDERSON. THANKS.
C
       SUBROUTINE POISSON1(NPAIR,ND,MD,ALPHAV,A,BETAV,B,RHO)
C ORIGINAL VERSION BY MARK R PEDERSON 4-July 1988
C$     use omp_lib
       use pot_dens,only : COULOMB,RHOG
       use mesh1,only : rmsh,nmsh
       use debug1
       use common2,only : IGGA, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI
       use common7,only : MODDEN
       use common8,only : REP
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:57 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: LM, LMX, MAX_ANG, NPAIR, ND, MD, MXPR, MPX, I, I1,
     & IANG, IC, IERR, IFMT, II, IP, IPOT, IPTS, IX, J, J1, K1,
     & KPTS_MAX, KPTS_MIN, L, L1, LL, LM_MAX, LMAX, LPTS, LPTS_BEG, M,
     & M1, M2, M3, MAXFMT, MAXN, MI, MM, MOM, MP, MPM, MPP1, MPP2,
     & MPP3, MPTS_CLS, MPTS_FAR, MX, MXX, MY, MYY, MZ, MZZ, N1, N2, N3,
     & NANG, NDMGGA, NGGA, NI, NP, NP_TOT, NPOLY, NPOT, NPP1, NPP2,
     & NPP3, NPTP1, NPTP2, NPTP3, NPTS, NPX, NTIMES, NX, NXX, NY, NYY,
     & NZ, NZZ
       REAL*8 :: SYMBOL , ALPHAV, A, BETAV, B, RHO, A0, AA, ACHRG,
     & ADDP, ADDR, ALPHA, ANG, ARG, ASYMP, ASYMPTOT, ATIME, B0, BB,
     & BETA, C, CCC, CGRAVITY, CGRTEMP, CHGD, COEFICIENT, COEFTMP,
     & COST1, COST2, COSTPROJ, D1X2, D1X4, D2X4, DAVG, DELT, DELTA,
     & DELTA_MIN, DELTAV, DIFF, DIFF_MAX, DIST, DLT, DMAX, DOMEGA,
     & EPSILON, EPSMULTI, ERP, ERROR, FACT, FACTOR, FACTOR1, FLDNR,
     & FMB, GT1, GT2, GT3, GT4, HNORM, PI, PI2, PI4, PIRC, POLY, POTLR,
     & PXY, Q, RATGGA, RCDELT, RCUT, RCUTSQ, RECIP, RECIPR, RHOC, RHOD,
     & RHOP, RI1, RMAX1, RMAX2, RMAX3, RMAX4, RMOMENT, RSQR, RTEST,
     & SKIPTEST, SOLANG, SSS, SSSS, TB, TDIFF, TF, TFMTTM,
     & THIRTYSOMETHING, TIME1, TIME2, TIME3, TIME4, TIMGGA, TIMGGN,
     & TMAX, TMGG1, TMGG2, TOLER, TPOLYS, TPOT, TPOT2, TRUNIT, TSETUP,
     & TWASTED, VLONG, W4, WASTE1, WASTE2, WT1, WT2, XX
       SAVE
       CHARACTER*6 FLAG
       LOGICAL GGA
c
c lm is the maximum value of l
c NOTE: lm must be >= 5 for poly to be dimensioned big enough
c
       PARAMETER (LM=06)
       PARAMETER (LMX=3*LM)
       PARAMETER (MXPR=MXPOISS)
       PARAMETER (MAX_ANG=((LMX+1)*(LMX+2))/2)
       PARAMETER (MPX=MAX_ANG)
       LOGICAL SKIP((LM+1)**2)
       LOGICAL FIRST
       LOGICAL CALC(35),SSCL(5)
C
C USE MIXPOT AS TEMPORARY STORAGE 
C
        REAL*8,ALLOCATABLE :: POTNL(:)
        INTEGER,ALLOCATABLE :: IPTR(:),JPTR(:)
C       COMMON/MIXPOT/POTNL(MAX_PTS),IPTR(MAX_PTS),JPTR(MAX_PTS)
C       COMMON/TMP1/ACOULOMB(MAX_PTS),ARHOG(MAX_PTS,KRHOG,MXSPN)
       COMMON/PTIME/TSETUP,TPOLYS,TFMTTM,TRUNIT,FLDNR,ASYMP,ACHRG
     &  ,ATIME,TWASTED,COSTPROJ 
C      
       DIMENSION RHOC(35,MXPR,MXSPN),RHOD(20,10,MXPR,MXSPN)
       DIMENSION AA(5,5,3,3,3),BB(5,5,5,5)
       DIMENSION CGRAVITY(4),W4(10,10,4)
       DIMENSION ALPHAV(MXPR),BETAV(MXPR),DELT(MXPR)
       DIMENSION A(3,MXPR),B(3,MXPR),C(3,MXPR),A0(3),B0(3),CCC(3)
       DIMENSION RHO(10,10,MXPR)
       DIMENSION RMOMENT((LM+1)**2)
       DIMENSION COEFICIENT(60,MXPR),DLT(MXPR)
       DIMENSION RHOP(3,3,3,3,3,3),ERP(MPX,5),ADDP(MPX),ADDR(MPX,20)
       DIMENSION CHGD(MPX)
       DIMENSION SSS(MAX_ANG,5),XX(MPX),POLY(MPX,(LM+1)**2)
       DIMENSION ANG(3,MAX_ANG),ASYMPTOT(MAX_ANG),Q(3,3*MAX_ANG)
       DIMENSION RECIP(3*MAX_ANG,0:LM),RECIPR(0:LM),RTEST(0:LM)
       DIMENSION DOMEGA(MAX_ANG),POTLR(MPX)
       DIMENSION TF(5),TB(5),FMB(5)
       DATA PI  /3.141592653589793D0/
       DATA PI2 /6.283185307179586D0/
       DATA PI4/12.566370614359172D0/
       DATA DELTA_MIN/2.0D0/
       DATA TOLER/1.0D-5/
C
C FOR USUAL CALCULATIONS (ANALYTICAL AND NUMERICAL SOLUTION)
C EPSILON MUST BE ZERO TO BE ON THE SAFE SIDE
C
c      DATA THIRTYSOMETHING/45.0D0/
       DATA THIRTYSOMETHING/30.0D0/
       DATA EPSMULTI/1.0D-6/
c      DATA THIRTYSOMETHING/20.0D0/
       DATA EPSILON/0.0D0/
C
C FORCE EVERYTHING TO BE DONE ANALYTICALLY 
C
c      DATA THIRTYSOMETHING/1.0D30/
c      DATA EPSILON/0.0D0/
C
       DATA FIRST/.TRUE./

       IF (LM .LT. 5) THEN
        write(6,*)'POISSON: LM MUST BE >= 5'
        CALL STOPIT
       END IF
       GGA=.FALSE.
       TIMGGA=0.0D0
       TIMGGN=0.0D0
       IF(IGGA(1).GT.0.OR.IGGA(2).GT.0)GGA=.TRUE.
C
C ALLLOCATE LOCAL ARRAYS
C
      IF(.NOT.ALLOCATED(POTNL)) ALLOCATE(POTNL(MAX_PTS),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'POISSON1:ERROR ALLOCATING POTNL'
      IF(.NOT.ALLOCATED(IPTR)) ALLOCATE(IPTR(MAX_PTS),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'POISSON1:ERROR ALLOCATING IPTR'
      IF(.NOT.ALLOCATED(JPTR)) ALLOCATE(JPTR(MAX_PTS),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'POISSON1:ERROR ALLOCATING JPTR'
C
C POREZAG 9/94 
C ALL ELEMENTS OF AA THAT ARE CONSTANTS NEED TO BE SET UP
C ONLY ONCE
C
C$ OMP PARALLEL
C$ OMP DO
       DO 760 LL=1,3
        DO 755 MM=1,3
         DO 750 NX=1,3
          DO 745 NY=1,5
           DO 740 NZ=1,5
            AA(NZ,NY,NX,MM,LL)=0.0D0
  740      CONTINUE
  745     CONTINUE
  750    CONTINUE
  755   CONTINUE
        AA(1,1,1,1,LL)=1.0D0
        AA(2,2,1,2,LL)=1.0D0
        AA(3,3,1,3,LL)=1.0D0
        AA(2,2,2,1,LL)=1.0D0
        AA(3,3,2,2,LL)=1.0D0
        AA(4,4,2,3,LL)=1.0D0
        AA(3,3,3,1,LL)=1.0D0
        AA(4,4,3,2,LL)=1.0D0
        AA(5,5,3,3,LL)=1.0D0
  760  CONTINUE
C$ OMP ENDDO
C$ OMP END PARALLEL
C
C IF DEBUG=.TRUE. , CHECK ASYMPTOTIC FORM OF POTENTIAL WITH ANALYTICAL
C FORM
C BOUNDS ON THE INCOMPLETE GAMMA FUNCTION:
C
       TF(1)=0.5D0*SQRT(PI)
       TF(2)=TF(1)*0.5D0
       TF(3)=TF(2)*1.5D0
       TF(4)=TF(3)*2.5D0
       TF(5)=TF(4)*3.5D0
       TB(1)=(TF(1)*(2*0+1))**(1.0D0/0.5D0)
       TB(2)=(TF(2)*(2*1+1))**(1.0D0/1.5D0)
       TB(3)=(TF(3)*(2*2+1))**(1.0D0/2.5D0)
       TB(4)=(TF(4)*(2*3+1))**(1.0D0/3.5D0)
       TB(5)=(TF(5)*(2*4+1))**(1.0D0/4.5D0)
C
C FOR T>TB(I) FM(T) < TF(I)/T**(M+0.5)
C FOR T<TB(I) FM(T) < 1/(2M+1)
C
       TDIFF=0.0D0
       DIFF_MAX=0.0D0
       IF(FIRST)THEN
        FIRST=.FALSE.
        IF(DEBUG)THEN
         OPEN(73,FILE='DEBUG',FORM='FORMATTED',STATUS='UNKNOWN') 
         REWIND(73)
        END IF
        LMAX=LMX
        CALL ANGMSH(MAX_ANG,LMAX,NANG,ANG,DOMEGA)
        IF(NANG.GT.MAX_ANG)THEN
         write(6,*)'POISSON: MAX_ANG MUST BE AT LEAST: ',NANG
         CALL STOPIT
        END IF
c 
c coordinates in ang are already normalized for harmonics routine
c
        CALL HARMONICS(MPX,NANG,LM,ANG,POLY,NPOLY)
        ERROR=0.0D0
c 
c no longer need to calculate RNORM but still check for ERROR 
c
        DO I=1,(LM+1)**2
         DO J=I,(LM+1)**2
          SSSS=0.0D0
          DO IANG=1,NANG
           SSSS=SSSS+POLY(IANG,I)*POLY(IANG,J)*DOMEGA(IANG)
          END DO
          IF(I.NE.J)THEN
           ERROR=MAX(ERROR,ABS(SSSS))
           IF(DEBUG.AND.(ABS(SSSS).GT.TOLER))write(6,*)I,J,SSSS
          ELSE
           IF (DEBUG) write(6,*)I,J,SSSS
c          RNORM(I)=1.0D0/SSSS
          END IF
         END DO
        END DO
        IF(DEBUG)THEN
          SOLANG=0.0D0
          DO IANG=1,NANG
            SOLANG=SOLANG+DOMEGA(IANG)
          END DO
          write(6,*)'SOLANG= ',SOLANG
        END IF
        IF (ERROR.GT.TOLER) THEN
         write(6,*)'POISSON: SPHERICAL HARMONICS TEST FAILED'
         CALL STOPIT
        END IF
       END IF
C
       IF(NMSH.LE.0)THEN
         write(6,*)'POISSON: NMSH IS <= ZERO'
         CALL STOPIT
       END IF
       IF(NMSH.GT.MAX_PTS)THEN
         write(6,*)'POISSON: MAX_PTS MUST BE AT LEAST:',NMSH
         CALL STOPIT
       END IF
       IF(NPAIR.GT.MXPR)THEN
         write(6,*)'POISSON: MXPR MUST BE AT LEAST: ',NPAIR
         CALL STOPIT
       END IF
       CALL GTTIME(TIME1)
       NP=2
       IF(ND.LT.5)NP=1
       IF(ND.LT.2)NP=0
       MP=2
       IF(MD.LT.5)MP=1
       IF(MD.LT.2)MP=0
C
C MAXIMUM DEGREE OF POLYNOMIAL:
C
       NP_TOT=NP+MP
       NPTP1=NP_TOT+1
       NPTP2=NP_TOT+2
       NPTP3=NP_TOT+3
       DELTA_MIN=1.0D30
       CGRAVITY(1)=0.0D0
       CGRAVITY(2)=0.0D0
       CGRAVITY(3)=0.0D0
       CGRAVITY(4)=0.0D0
       DO 200 IP=1,NPAIR,NSPN
        ALPHA=ALPHAV(IP)
        BETA =BETAV (IP)
        IF(MODDEN.EQ.1)THEN
         CALL GTTIME(TMGG1)
         DO ISPN=1,NSPN
           CALL GTDNCF(GGA,ALPHA,A(1,IP),BETA,B(1,IP)
     &                 ,RHO(1,1,IP+ISPN-1),RHOC(1,IP,ISPN)
     &                 ,RHOD(1,1,IP,ISPN),DLT(IP),CCC)
         END DO
         CALL GTTIME(TMGG2)
         TIMGGA=TIMGGA+TMGG2-TMGG1
        END IF
        DO I=1,10
         DO J=1,10
          RHO(J,I,IP)=RHO(J,I,IP)+RHO(J,I,IP+NSPN-1)
         END DO
        END DO
        IF(NPAIR.NE.1)THEN
         CALL GTTIME(GT1)
         CALL GINTED(ALPHA,BETA,A(1,IP),B(1,IP),W4)
         DO IC=1,4
          CGRTEMP=0.0D0
          DO J=1,10
           DO I=1,10
            CGRTEMP=CGRTEMP+W4(I,J,IC)*RHO(I,J,IP)
           END DO
          END DO
          CGRAVITY(IC)=CGRAVITY(IC)+CGRTEMP
         END DO
         CALL GTTIME(GT2)
        END IF
        IF((ALPHA .LT. 0.0D0).OR.(BETA .LT. 0.0D0))THEN
         write(6,*)'POISSON: ALPHA OR BETA < 0'
         CALL STOPIT
        END IF
        DELTA=ALPHA+BETA
        IF(DELTA.LT.DELTA_MIN)DELTA_MIN=DELTA
        RCDELT=1.0D0/DELTA
        ARG=ALPHA*BETA*RCDELT*((A(1,IP)-B(1,IP))**2
     &     +(A(2,IP)-B(2,IP))**2+(A(3,IP)-B(3,IP))**2)
        FACTOR=PI2*EXP(-ARG)*RCDELT
C
C DVP: FIX G77 PROBLEM (PROGRAM WOULD STOP DUE TO TINY INACCURACY)
C
C       IF (EXARG .GT. 1.0D0)THEN
C        write(6,*)'POISSON: FACTOR TOO LARGE: ',FACTOR,ARG,
C        CALL STOPIT
C       END IF
C
        DO 15 I=1,3
         C(I,IP)=(ALPHA*A(I,IP)+BETA*B(I,IP))*RCDELT
         A0(I)=A(I,IP)-C(I,IP)
         B0(I)=B(I,IP)-C(I,IP)
   15   CONTINUE
C
C MOVE RHO TO RHOP:
C

C$ OMP PARALLEL
C$ OMP DO
        DO I1=1,3
         DO J1=1,3
          DO K1=1,3
           DO L1=1,3
            DO M1=1,3
             DO N1=1,3
              RHOP(N1,M1,L1,K1,J1,I1)=0.0D0
             END DO
            END DO
           END DO
          END DO
         END DO
        END DO
C$ OMP ENDDO
C$ OMP END PARALLEL

        DO NXX=0,NP
         DO NYY=0,NP-NXX
          DO NZZ=0,NP-NXX-NYY
           NI=1+NXX*(NXX+4*NYY)+4*NZZ*(NXX+NYY)+3*NZZ
           NI=NI+(NYY*(NYY+3))/2
           N1=NXX+1
           N2=NYY+1
           N3=NZZ+1
           DO MXX=0,MP
            DO MYY=0,MP-MXX
             DO MZZ=0,MP-MXX-MYY
              MI=1+MXX*(MXX+4*MYY)+4*MZZ*(MXX+MYY)+3*MZZ
              MI=MI+(MYY*(MYY+3))/2
              M1=MXX+1
              M2=MYY+1
              M3=MZZ+1
              RHOP(M3,N3,M2,N2,M1,N1)=RHO(NI,MI,IP)
             END DO
            END DO
           END DO
          END DO
         END DO
        END DO

C$ OMP PARALLEL
C$ OMP DO
        DO MM=1,5
         DO NX=1,5
          DO NY=1,5
           DO NZ=1,5
            BB(NZ,NY,NX,MM)=0.0D0
           END DO
          END DO
         END DO
        END DO
C$ OMP ENDDO
C$ OMP END PARALLEL

        DO 10 I=1,3
         CALL POLYX(A0(I),B0(I),DELTA,AA(1,1,1,1,I))
   10   CONTINUE
        NPP1=NP+1
        NPP2=NP+2
        NPP3=NP+3
        MPP1=MP+1
        MPP2=MP+2
        MPP3=MP+3
        CALL GTTIME(GT3)
        DO 85 M1=1,MPP1
        DO 85 N1=1,NPP1
        DO 85 NX=1,NPTP1
         DO 80 MX=NX,NPTP1
          IF(AA(NX,MX,N1,M1,1).EQ.0.0D0)GO TO 79
          DO 75 M2=1,MPP2-M1
          DO 75 N2=1,NPP2-N1
          DO 75 NY=1,NPTP2-NX
           DO 70 MY=NY,NPTP2-MX
            IF(AA(NY,MY,N2,M2,2).EQ.0.0D0)GO TO 69
            PXY=AA(NX,MX,N1,M1,1)*AA(NY,MY,N2,M2,2)
            DO 65 M3=1,MPP3-M1-M2
             DO 60 N3=1,NPP3-N1-N2
              FACTOR1=RHOP(M3,N3,M2,N2,M1,N1)*PXY
              IF(FACTOR1.EQ.0.0D0)GO TO 59
              DO 55 NZ=1,NPTP3-NX-NY
               MM=MX+MY-2
               DO 50 MZ=NZ,NPTP3-MX-MY
                BB(MM+MZ,NX,NY,NZ)=BB(MM+MZ,NX,NY,NZ)
     &                            +FACTOR1*AA(NZ,MZ,N3,M3,3)
   50          CONTINUE
   55         CONTINUE
   59         CONTINUE
   60        CONTINUE
   65       CONTINUE
   69       CONTINUE
   70      CONTINUE
   75     CONTINUE
   79     CONTINUE
   80    CONTINUE
   85   CONTINUE
        CALL GTTIME(GT4)
        TWASTED=TWASTED+GT4-GT3
        DO 86 I=1,46
         COEFICIENT(I,IP)=0.0D0
   86   CONTINUE
C
C POLYNOMIALS OF DEGREE:  0
C
        COEFICIENT(  1,IP)= BB( 1, 1, 1, 1)
        COEFICIENT(  2,IP)= BB( 2, 1, 1, 1)
        COEFICIENT(  3,IP)= BB( 3, 1, 1, 1)
        MAXN=3
        IF(NP_TOT.GE.1)THEN
C
C POLYNOMIALS OF DEGREE 1
C
         COEFICIENT(  4,IP)= BB( 2, 1, 1, 2)
         COEFICIENT(  5,IP)= BB( 2, 1, 2, 1)
         COEFICIENT(  6,IP)= BB( 2, 2, 1, 1)
         COEFICIENT(  7,IP)= BB( 3, 1, 1, 2)
         COEFICIENT(  8,IP)= BB( 3, 1, 2, 1)
         COEFICIENT(  9,IP)= BB( 3, 2, 1, 1)
         MAXN=9
        END IF
        IF(NP_TOT.GE.2)THEN
C
C POLYNOMIALS OF DEGREE:  2
C
         COEFICIENT( 10,IP)= BB( 3, 1, 1, 3)
         COEFICIENT( 11,IP)= BB( 3, 1, 2, 2)
         COEFICIENT( 12,IP)= BB( 3, 1, 3, 1)
         COEFICIENT( 13,IP)= BB( 3, 2, 1, 2)
         COEFICIENT( 14,IP)= BB( 3, 2, 2, 1)
         COEFICIENT( 15,IP)= BB( 3, 3, 1, 1)
         COEFICIENT( 16,IP)= BB( 4, 1, 1, 3)
         COEFICIENT( 17,IP)= BB( 4, 1, 2, 2)
         COEFICIENT( 18,IP)= BB( 4, 1, 3, 1)
         COEFICIENT( 19,IP)= BB( 4, 2, 1, 2)
         COEFICIENT( 20,IP)= BB( 4, 2, 2, 1)
         COEFICIENT( 21,IP)= BB( 4, 3, 1, 1)
         MAXN=21
        END IF
        IF(NP_TOT.GE.3)THEN
C
C POLYNOMIALS OF DEGREE:  3
C
         COEFICIENT( 22,IP)= BB( 4, 1, 1, 4)
         COEFICIENT( 23,IP)= BB( 4, 1, 2, 3)
         COEFICIENT( 24,IP)= BB( 4, 1, 3, 2)
         COEFICIENT( 25,IP)= BB( 4, 1, 4, 1)
         COEFICIENT( 26,IP)= BB( 4, 2, 1, 3)
         COEFICIENT( 27,IP)= BB( 4, 2, 2, 2)
         COEFICIENT( 28,IP)= BB( 4, 2, 3, 1)
         COEFICIENT( 29,IP)= BB( 4, 3, 1, 2)
         COEFICIENT( 30,IP)= BB( 4, 3, 2, 1)
         COEFICIENT( 31,IP)= BB( 4, 4, 1, 1)
         MAXN=31
        END IF
        IF(NP_TOT.GE.4)THEN
C
C POLYNOMIALS OF DEGREE:  4
C
         COEFICIENT( 32,IP)= BB( 5, 1, 1, 5)
         COEFICIENT( 33,IP)= BB( 5, 1, 2, 4)
         COEFICIENT( 34,IP)= BB( 5, 1, 3, 3)
         COEFICIENT( 35,IP)= BB( 5, 1, 4, 2)
         COEFICIENT( 36,IP)= BB( 5, 1, 5, 1)
         COEFICIENT( 37,IP)= BB( 5, 2, 1, 4)
         COEFICIENT( 38,IP)= BB( 5, 2, 2, 3)
         COEFICIENT( 39,IP)= BB( 5, 2, 3, 2)
         COEFICIENT( 40,IP)= BB( 5, 2, 4, 1)
         COEFICIENT( 41,IP)= BB( 5, 3, 1, 3)
         COEFICIENT( 42,IP)= BB( 5, 3, 2, 2)
         COEFICIENT( 43,IP)= BB( 5, 3, 3, 1)
         COEFICIENT( 44,IP)= BB( 5, 4, 1, 2)
         COEFICIENT( 45,IP)= BB( 5, 4, 2, 1)
         COEFICIENT( 46,IP)= BB( 5, 5, 1, 1)
         MAXN=46
        END IF
        DO 40 I=1, MAXN
         COEFICIENT(I,IP)=FACTOR*COEFICIENT(I,IP)
   40   CONTINUE
        DELT(IP)=DELTA
        C(1,IP)=-C(1,IP)
        C(2,IP)=-C(2,IP)
        C(3,IP)=-C(3,IP)
  200  CONTINUE
C
       IF(CGRAVITY(4).NE.0.0D0.AND.NPAIR.GT.1)THEN
        CGRAVITY(1)=CGRAVITY(1)/CGRAVITY(4)
        CGRAVITY(2)=CGRAVITY(2)/CGRAVITY(4)
        CGRAVITY(3)=CGRAVITY(3)/CGRAVITY(4)
       ELSE
        DO IX=1,3
         CGRAVITY(IX)=0.0D0
        END DO
        DO IP=1,NPAIR,NSPN
         DO IX=1,3
          CGRAVITY(IX)=CGRAVITY(IX)-C(IX,IP)
         END DO
        END DO
        DO IX=1,3
         CGRAVITY(IX)=(NSPN*CGRAVITY(IX))/NPAIR
        END DO
       END IF
C
       DMAX=0.0D0 
       DO IP=1,NPAIR,NSPN
        DIST=SQRT((CGRAVITY(1)+C(1,IP))**2+
     &            (CGRAVITY(2)+C(2,IP))**2+
     &            (CGRAVITY(3)+C(3,IP))**2)
        IF(DIST.GT.DMAX)DMAX=DIST
       END DO
       RCUT=SQRT(THIRTYSOMETHING/DELTA_MIN)+DMAX
C
C MODIFIED WAY TO DETERMINE RCUT INTRODUCED BY DVP 05/98
C
       PIRC=1.0D0/(4*ATAN(1.0D0))
       FACT=SQRT(PIRC/(PIRC+LM+1))
       DIST=2*DMAX/(EPSMULTI/FACT)**(1.0D0/(LM+1))
       RCUT=MAX(RCUT,DIST)
C
       DO IANG=1,NANG
        ASYMPTOT(IANG)=0.0D0
        Q(1,IANG)=CGRAVITY(1)+ANG(1,IANG)*RCUT
        Q(2,IANG)=CGRAVITY(2)+ANG(2,IANG)*RCUT
        Q(3,IANG)=CGRAVITY(3)+ANG(3,IANG)*RCUT
       END DO
       CGRAVITY(1)=-CGRAVITY(1)
       CGRAVITY(2)=-CGRAVITY(2)
       CGRAVITY(3)=-CGRAVITY(3)
       RCUTSQ=RCUT*RCUT
C
C FIND NUMBER OF POINTS THAT ARE IN THE ASYMPTOTIC REGION:
C
       MPTS_CLS=0
       MPTS_FAR=0
       CALL GTTIME(WT1)
       DO IPTS=1,NMSH
        RSQR= (RMSH(1,IPTS)+CGRAVITY(1))**2+
     &        (RMSH(2,IPTS)+CGRAVITY(2))**2+
     &        (RMSH(3,IPTS)+CGRAVITY(3))**2
        IF(RSQR.GT.RCUTSQ)THEN
         MPTS_FAR=MPTS_FAR+1
         IPTR(MPTS_FAR)=IPTS
        ELSE
         MPTS_CLS=MPTS_CLS+1
         JPTR(MPTS_CLS)=IPTS
        END IF
       END DO
       CALL GTTIME(WT2)
       NPOT=2
C
C      PROJECTED COSTS:
C
       COST1=1.36D0*(MPTS_CLS+NANG)*NPAIR + 0.99D0*MPTS_FAR
       COST2=1.36D0*(MPTS_FAR+MPTS_CLS)*NPAIR
       IF(COST2.LT.COST1)THEN
        COST1=COST2
        NPOT=1
        MPTS_FAR=0
        MPTS_CLS=NMSH
        DO IPTS=1,NMSH
         JPTR(IPTS)=IPTS
        END DO
        FLDNR=FLDNR-0.001D0*NANG*NPAIR
       END IF
       ASYMP=ASYMP+0.001D0*MPTS_FAR
       FLDNR=FLDNR+0.001D0*(MPTS_CLS+NANG)*NPAIR
       COSTPROJ=COSTPROJ+COST1*1.0D-6
       CALL GTTIME(TIME2)
       TSETUP=TSETUP+TIME2-TIME1
C
C END OF SETTING UP:
C
       DO 600 IPOT=1,NPOT
        CALL GTTIME(TPOT)
        IF(IPOT.EQ.1)THEN
         IF(DEBUG)THEN
          NPTS=NMSH
         ELSE
          NPTS=MPTS_CLS
         END IF
        ELSE
         NPTS=NANG  
        END IF
        IF(IPOT.EQ.1)THEN
         DO IPTS=1,NPTS
          POTNL(IPTS)=0.0D0
C         CHGDN(IPTS)=0.0D0
         END DO
        END IF
C
        LPTS_BEG=0
 250    CONTINUE
        NPX=MIN(MPX,NPTS-LPTS_BEG)
        DO 550 IP=1,NPAIR,NSPN
         CALL GTTIME(TIME1)
         IF(IPOT.EQ.1)THEN
          IF(DEBUG)THEN
           DO 398 IPTS=1,NPX
            POLY(IPTS,4)=RMSH(1,LPTS_BEG+IPTS)+C(1,IP)
            POLY(IPTS,3)=RMSH(2,LPTS_BEG+IPTS)+C(2,IP)
            POLY(IPTS,2)=RMSH(3,LPTS_BEG+IPTS)+C(3,IP)
            POLY(IPTS,1)=1.0D0
 398       CONTINUE
          ELSE
           DO 399 IPTS=1,NPX
            POLY(IPTS,4)=RMSH(1,JPTR(LPTS_BEG+IPTS))+C(1,IP)
            POLY(IPTS,3)=RMSH(2,JPTR(LPTS_BEG+IPTS))+C(2,IP)
            POLY(IPTS,2)=RMSH(3,JPTR(LPTS_BEG+IPTS))+C(3,IP)
            POLY(IPTS,1)=1.0D0
 399       CONTINUE
          END IF
         ELSE
          DO 400 IPTS=1,NPX
           POLY(IPTS,4)=Q(1,LPTS_BEG+IPTS)+C(1,IP)
           POLY(IPTS,3)=Q(2,LPTS_BEG+IPTS)+C(2,IP)
           POLY(IPTS,2)=Q(3,LPTS_BEG+IPTS)+C(3,IP)
           POLY(IPTS,1)=1.0D0
 400      CONTINUE
         END IF
         DO 406 IPTS=1,NPX
          POLY(IPTS,5)=POLY(IPTS,2)**2
          POLY(IPTS,7)=POLY(IPTS,3)**2
          POLY(IPTS,10)=POLY(IPTS,4)**2
 406     CONTINUE
         TMAX=0.0D0 
         DO 432 IPTS=1,NPX
          XX(IPTS)=DELT(IP)*(POLY(IPTS,5)+POLY(IPTS,7)+POLY(IPTS,10))
          TMAX=MAX(TMAX,XX(IPTS))
 432     CONTINUE
        DO IPTS=1,NPX
          CHGD(IPTS)=EXP(-XX(IPTS))
        END DO
         RMAX1=TMAX/DELT(IP)
         RMAX2=RMAX1*RMAX1
         RMAX3=RMAX1*RMAX2
         RMAX4=RMAX2*RMAX2
         FMB(1)=SQRT(TMAX)
         FMB(2)=FMB(1)*TMAX
         FMB(3)=FMB(2)*TMAX
         FMB(4)=FMB(3)*TMAX
         FMB(5)=FMB(4)*TMAX
         DO M=0,4
          IF(TMAX.GT.TB(M+1))THEN
           FMB(M+1)=TF(M+1)/FMB(M+1)
          ELSE
           FMB(M+1)=1.0D0/(2.0D0*M+1.0D0)
          END IF
         END DO
         DO I=1,10
          CALC(I)=.TRUE.
         END DO
C
C OLD VERSION BELOW HAS BEEN REPLACED BY CALC(I)=.TRUE.
C
C        DO I=11,35
C         CALC(I)=.FALSE.
C        END DO
         DO I=11,35
          CALC(I)=.TRUE.
         END DO
         SSCL(1)=.TRUE.
         SSCL(2)=.TRUE.
         DO I=3,5
          SSCL(I)=.FALSE.
         END DO
         DO I=2,4 
          IF(ABS(COEFICIENT(I+5,IP)*RMAX1*FMB(3)).GT.EPSILON)THEN
           SSCL(3)=.TRUE.
          END IF
         END DO
         DO I=5,10
          IF(ABS(COEFICIENT(I+5,IP)*RMAX2*FMB(3)).GT.EPSILON)THEN
           SSCL(3)=.TRUE.
          END IF
         END DO   
         DO I=5,10
          IF(ABS(COEFICIENT(I+11,IP)*RMAX2*FMB(4)).GT.EPSILON)THEN
           SSCL(4)=.TRUE.
          END IF
         END DO
         DO I=11,20
          IF(ABS(COEFICIENT(I+11,IP)*RMAX3*FMB(4)).GT.EPSILON)THEN
           CALC(I)=.TRUE.
           SSCL(4)=.TRUE.
          END IF
         END DO   
         DO I=21,35
          IF(ABS(COEFICIENT(I+11,IP)*RMAX4*FMB(5)).GT.EPSILON)THEN
           CALC(I)=.TRUE.
           SSCL(5)=.TRUE.
          END IF
         END DO
         MAXFMT=3
         IF(SSCL(4))MAXFMT=4
         IF(SSCL(5))MAXFMT=5
         IF(NP_TOT.GT.1)THEN
          DO 402 IPTS=1,NPX
           POLY(IPTS,6)=POLY(IPTS,3)*POLY(IPTS,2)
           POLY(IPTS,8)=POLY(IPTS,4)*POLY(IPTS,2)
           POLY(IPTS,9)=POLY(IPTS,4)*POLY(IPTS,3)
 402      CONTINUE
         END IF
         IF(NP_TOT.GT.2)THEN
          IF(CALC(11))THEN
           DO 407 IPTS=1,NPX
            POLY(IPTS,11)=POLY(IPTS,2)*POLY(IPTS,5)
 407       CONTINUE
          END IF
          IF(CALC(12))THEN
           DO 408 IPTS=1,NPX
            POLY(IPTS,12)=POLY(IPTS,3)*POLY(IPTS,5)
 408       CONTINUE
          END IF
          IF(CALC(13))THEN
           DO 409 IPTS=1,NPX
            POLY(IPTS,13)=POLY(IPTS,2)*POLY(IPTS,7)
 409       CONTINUE
          END IF
          IF(CALC(14))THEN
           DO 410 IPTS=1,NPX
            POLY(IPTS,14)=POLY(IPTS,3)*POLY(IPTS,7)
 410       CONTINUE
          END IF
          IF(CALC(15))THEN
           DO 411 IPTS=1,NPX
            POLY(IPTS,15)=POLY(IPTS,4)*POLY(IPTS,5)
 411       CONTINUE
          END IF
          IF(CALC(16))THEN
           DO 412 IPTS=1,NPX
            POLY(IPTS,16)=POLY(IPTS,4)*POLY(IPTS,6)
 412       CONTINUE
          END IF
          IF(CALC(17))THEN
           DO 413 IPTS=1,NPX
            POLY(IPTS,17)=POLY(IPTS,4)*POLY(IPTS,7)
 413       CONTINUE
          END IF
          IF(CALC(18))THEN
           DO 414 IPTS=1,NPX
            POLY(IPTS,18)=POLY(IPTS,2)*POLY(IPTS,10)
 414       CONTINUE
          END IF
          IF(CALC(19))THEN
           DO 415 IPTS=1,NPX
            POLY(IPTS,19)=POLY(IPTS,3)*POLY(IPTS,10)
 415       CONTINUE
          END IF
          IF(CALC(20))THEN
           DO 416 IPTS=1,NPX
            POLY(IPTS,20)=POLY(IPTS,4)*POLY(IPTS,10)
 416       CONTINUE
          END IF
         END IF
         IF(NP_TOT.GT.3)THEN
          IF(CALC(21))THEN
           DO 417 IPTS=1,NPX
            POLY(IPTS,21)=POLY(IPTS,5)**2
 417       CONTINUE
          END IF
          IF(CALC(22))THEN
           DO 420 IPTS=1,NPX
            POLY(IPTS,22)=POLY(IPTS,5)*POLY(IPTS,6)
 420       CONTINUE
          END IF
          IF(CALC(23))THEN
           DO 418 IPTS=1,NPX
            POLY(IPTS,23)=POLY(IPTS,7)*POLY(IPTS,5)
 418       CONTINUE
          END IF
          IF(CALC(24))THEN
           DO 421 IPTS=1,NPX
            POLY(IPTS,24)=POLY(IPTS,7)*POLY(IPTS,6)
 421       CONTINUE
          END IF
          IF(CALC(25))THEN
           DO 422 IPTS=1,NPX
            POLY(IPTS,25)=POLY(IPTS,7)**2
 422       CONTINUE
          END IF
          IF(CALC(26))THEN
           DO 423 IPTS=1,NPX
            POLY(IPTS,26)=POLY(IPTS,5)*POLY(IPTS,8)
 423       CONTINUE
          END IF
          IF(CALC(27))THEN
           DO 419 IPTS=1,NPX
            POLY(IPTS,27)=POLY(IPTS,9)*POLY(IPTS,5)
 419       CONTINUE
          END IF
          IF(CALC(28))THEN
           DO 424 IPTS=1,NPX
            POLY(IPTS,28)=POLY(IPTS,9)*POLY(IPTS,6)
 424       CONTINUE
          END IF
          IF(CALC(29))THEN
           DO 425 IPTS=1,NPX
            POLY(IPTS,29)=POLY(IPTS,7)*POLY(IPTS,9)
 425       CONTINUE
          END IF
          IF(CALC(30))THEN
           DO 426 IPTS=1,NPX
            POLY(IPTS,30)=POLY(IPTS,10)*POLY(IPTS,5)
 426       CONTINUE
          END IF
          IF(CALC(31))THEN
           DO 427 IPTS=1,NPX
            POLY(IPTS,31)=POLY(IPTS,10)*POLY(IPTS,6)
 427       CONTINUE
          END IF
          IF(CALC(32))THEN
           DO 428 IPTS=1,NPX
            POLY(IPTS,32)=POLY(IPTS,10)*POLY(IPTS,7)
 428       CONTINUE
          END IF
          IF(CALC(33))THEN
           DO 429 IPTS=1,NPX
            POLY(IPTS,33)=POLY(IPTS,8)*POLY(IPTS,10)
 429       CONTINUE
          END IF
          IF(CALC(34))THEN
           DO 430 IPTS=1,NPX
            POLY(IPTS,34)=POLY(IPTS,9)*POLY(IPTS,10)
 430       CONTINUE
          END IF
          IF(CALC(35))THEN
           DO 431 IPTS=1,NPX
            POLY(IPTS,35)=POLY(IPTS,10)**2
 431       CONTINUE
          END IF
         END IF
C
         IF((MODDEN.EQ.1).AND.(IPOT.EQ.1))THEN
          DO ISPN=1,NSPN
           NGGA=1
           IF(GGA)NGGA=10
           DO II=1,NGGA
            DO IPTS=1,NPX
             ADDR(IPTS,II)=0.0D0
            END DO
           END DO 
           DO MPM=1,35
            IF (ABS(RHOC(MPM,IP,ISPN)).GT.1.0D-14) THEN
             DO IPTS=1,NPX
              ADDR(IPTS,1)=ADDR(IPTS,1)+POLY(IPTS,MPM)*RHOC(MPM,IP,ISPN)
             END DO
            END IF
           END DO
           IF(GGA)THEN
            CALL GTTIME(TMGG1)
            DO MPM=1,20
             DO IX=2,4
              IF(ABS(RHOD(MPM,IX,IP,ISPN)).GT.1.0D-14)THEN
               DO IPTS=1,NPX
                ADDR(IPTS,IX)=ADDR(IPTS,IX)
     &                       +POLY(IPTS,MPM)*RHOD(MPM,IX,IP,ISPN)
               END DO
              END IF
             END DO
            END DO
            DO MPM=1,10
             DO IX=5,10
              IF(ABS(RHOD(MPM,IX,IP,ISPN)).GT.1.0D-14)THEN
               DO IPTS=1,NPX
                ADDR(IPTS,IX)=ADDR(IPTS,IX)
     &                       +POLY(IPTS,MPM)*RHOD(MPM,IX,IP,ISPN)
               END DO
              END IF
             END DO
            END DO
            D1X4=DLT(IP)*4.0D0
            D2X4=DLT(IP)*DLT(IP)*4.0D0 
            D1X2=DLT(IP)*2.0D0
            DO IPTS=1,NPX
             ADDR(IPTS,5)=ADDR(IPTS,5)-D1X4*POLY(IPTS,4)*ADDR(IPTS,2)
     &                   +(D2X4*POLY(IPTS,10)-D1X2)*ADDR(IPTS,1)
            END DO
            DO IPTS=1,NPX
             ADDR(IPTS,6)=ADDR(IPTS,6)-D1X4*POLY(IPTS,3)*ADDR(IPTS,3)
     &                   +(D2X4*POLY(IPTS, 7)-D1X2)*ADDR(IPTS,1)
            END DO
            DO IPTS=1,NPX
             ADDR(IPTS,7)=ADDR(IPTS,7)-D1X4*POLY(IPTS,2)*ADDR(IPTS,4)
     &                   +(D2X4*POLY(IPTS, 5)-D1X2)*ADDR(IPTS,1)
            END DO
            DO IPTS=1,NPX
            ADDR(IPTS,8)=ADDR(IPTS,8)-D1X2*POLY(IPTS,3)*ADDR(IPTS,2) !YX
     &                               -D1X2*POLY(IPTS,4)*ADDR(IPTS,3) !XY
     &                  + D2X4*POLY(IPTS,9)*ADDR(IPTS,1)             !XY
            END DO
            DO IPTS=1,NPX
            ADDR(IPTS,9)=ADDR(IPTS,9)-D1X2*POLY(IPTS,2)*ADDR(IPTS,2) !ZX
     &                               -D1X2*POLY(IPTS,4)*ADDR(IPTS,4) !XZ
     &                  + D2X4*POLY(IPTS,8)*ADDR(IPTS,1)             !XZ
            END DO
            DO IPTS=1,NPX
           ADDR(IPTS,10)=ADDR(IPTS,10)-D1X2*POLY(IPTS,3)*ADDR(IPTS,4)!YZ
     &                                -D1X2*POLY(IPTS,2)*ADDR(IPTS,3)!ZY
     &                  + D2X4*POLY(IPTS,6)*ADDR(IPTS,1)             !YZ
            END DO
            DO IPTS=1,NPX
             ADDR(IPTS,2)=ADDR(IPTS,2)-D1X2*POLY(IPTS,4)*ADDR(IPTS,1)
            END DO
            DO IPTS=1,NPX
             ADDR(IPTS,3)=ADDR(IPTS,3)-D1X2*POLY(IPTS,3)*ADDR(IPTS,1)
            END DO
            DO IPTS=1,NPX
             ADDR(IPTS,4)=ADDR(IPTS,4)-D1X2*POLY(IPTS,2)*ADDR(IPTS,1)
            END DO
            CALL GTTIME(TMGG2)
            TIMGGA=TIMGGA+(TMGG2-TMGG1)
            IF((TIMGGA.GT.TIMGGN).AND.DEBUG)THEN
             RATGGA=TIMGGA/TMGG2
             TIMGGN=3.0D0*TIMGGA
             write(6,*)'TIMGGA, RATGGA:',TIMGGA,RATGGA
            END IF
           END IF
C
           CALL GTTIME(TMGG1)
           NDMGGA=1
           IF (GGA) NDMGGA=10
           DO II=1,NDMGGA
            DO IPTS=1,NPX
             ADDR(IPTS,II)=CHGD(IPTS)*ADDR(IPTS,II)
            END DO
           END DO
           DO II=1,NDMGGA
            DO IPTS=1,NPX
             RHOG(JPTR(LPTS_BEG+IPTS),II,ISPN)=
     &        RHOG(JPTR(LPTS_BEG+IPTS),II,ISPN)+ADDR(IPTS,II)
            END DO
           END DO
           CALL GTTIME(TMGG2)
           TIMGGA=TIMGGA+(TMGG2-TMGG1)
          END DO
         END IF
C
C POLY(1 ) 1.0D0
C POLY(2 ) Z
C POLY(3 ) Y
C POLY(4 ) X
C POLY(5 ) ZZ
C POLY(6 ) YZ
C POLY(7 ) YY
C POLY(8 ) XZ
C POLY(9 ) XY
C POLY(10) XX
C POLY(11) ZZZ
C POLY(12) YZZ
C POLY(13) YYZ
C POLY(14) YYY
C POLY(15) XZZ
C POLY(16) XYZ
C POLY(17) XYY  
C POLY(18) XXZ 
C POLY(19) XXY 
C POLY(20) XXX 
C POLY(21) ZZZZ
C POLY(22) YZZZ
C POLY(23) YYZZ
C POLY(24) YYYZ
C POLY(25) YYYY
C POLY(26) XZZZ
C POLY(27) XYZZ
C POLY(28) XYYZ
C POLY(29) XYYY
C POLY(30) XXZZ
C POLY(31) XXYZ
C POLY(32) XXYY
C POLY(33) XXXZ
C POLY(34) XXXY
C POLY(35) XXXX
C
         CALL GTTIME(TIME2)
         IF(IPOT.EQ.1)TPOLYS=TPOLYS+TIME2-TIME1
         IF(MAXFMT.EQ.5)THEN
          CALL FFFMTC(MPX,MPX,XX,CHGD,SSS)
         ELSE IF(MAXFMT.EQ.4)THEN
          CALL FFFMT4(MPX,MPX,XX,CHGD,SSS)
         ELSE IF(MAXFMT.EQ.3)THEN
          CALL FFFMT3(MPX,MPX,XX,CHGD,SSS)
         END IF
         CALL GTTIME(TIME3)
         IF(IPOT.EQ.1)TFMTTM=TFMTTM+TIME3-TIME2
         DO 433 I=1,3
         DO 433 IPTS=1,NPX
          ERP(IPTS,I)=COEFICIENT(I,IP)
  433    CONTINUE
         DO 434 I=4,5
         DO 434 IPTS=1,NPX
          ERP(IPTS,I)=0.0D0
  434    CONTINUE
         DO 45 I=2,4
         COEFTMP=COEFICIENT(I+2,IP)
         DO 45 IPTS=1,NPX
          ERP(IPTS,2)=ERP(IPTS,2)+POLY(IPTS,I)*COEFTMP
   45    CONTINUE
C
C POLYNOMIALS OF DEGREE:  1
C
         DO 94 I=2,4 
          IF(CALC(I).AND.SSCL(3))THEN
           COEFTMP=COEFICIENT(I+5,IP)
           DO 93 IPTS=1,NPX
            ERP(IPTS,3)=ERP(IPTS,3)+POLY(IPTS,I)*COEFTMP
   93      CONTINUE
          END IF
   94    CONTINUE
C
C POLYNOMIALS OF DEGREE:  2
C
         DO 96 I=5,10
          IF(CALC(I).AND.SSCL(3))THEN
           COEFTMP=COEFICIENT(I+5,IP)
           DO 95 IPTS=1,NPX
            ERP(IPTS,3)=ERP(IPTS,3)+POLY(IPTS,I)*COEFTMP
   95      CONTINUE
          END IF
   96    CONTINUE
C
         DO 101 I=5,10
          IF(CALC(I).AND.SSCL(4))THEN
           COEFTMP=COEFICIENT(I+11,IP)
           DO 100 IPTS=1,NPX
            ERP(IPTS,4)=ERP(IPTS,4)+POLY(IPTS,I)*COEFTMP
  100      CONTINUE
          END IF
  101    CONTINUE
C
C POLYNOMIALS OF DEGREE:  3
C
         IF(NP_TOT.GT.2)THEN
          DO 106 I=11,20
           IF(CALC(I).AND.SSCL(4))THEN
            COEFTMP=COEFICIENT(I+11,IP)
            DO 105 IPTS=1,NPX
             ERP(IPTS,4)=ERP(IPTS,4)+POLY(IPTS,I)*COEFTMP
  105       CONTINUE
           END IF
  106     CONTINUE
C
C POLYNOMIALS OF DEGREE:  4
C
          IF(NP_TOT.GT.3)THEN
           DO 111 I=21,35
            IF(CALC(I).AND.SSCL(5))THEN
             COEFTMP=COEFICIENT(I+11,IP)
             DO 110 IPTS=1,NPX
              ERP(IPTS,5)=ERP(IPTS,5)+POLY(IPTS,I)*COEFTMP
  110        CONTINUE
            END IF
  111      CONTINUE
          END IF
         END IF
C
         DO 497 IPTS=1,NPX
          ADDP(IPTS)=0.0D0
 497     CONTINUE
         DO 498 IFMT=1,5
          IF(SSCL(IFMT))THEN
           DO IPTS=1,NPX
            ADDP(IPTS)=ADDP(IPTS)+ERP(IPTS,IFMT)*SSS(IPTS,IFMT)
           END DO
          END IF
 498     CONTINUE
         CALL GTTIME(WASTE1)
C
         IF(IPOT.EQ.1)THEN
          DO 499 IPTS=1,NPX
           POTNL(LPTS_BEG+IPTS)=POTNL(LPTS_BEG+IPTS)+ADDP(IPTS)
C          CHGDN(LPTS_BEG+IPTS)=CHGDN(LPTS_BEG+IPTS)+ADDR(IPTS,1)
 499      CONTINUE
         ELSE
          DO 500 IPTS=1,NPX
           ASYMPTOT(LPTS_BEG+IPTS)=ASYMPTOT(LPTS_BEG+IPTS)+ADDP(IPTS)
 500      CONTINUE
         END IF
         CALL GTTIME(TIME4)
         IF(IPOT.EQ.1)TRUNIT=TRUNIT+TIME4-TIME3
 550    CONTINUE
        LPTS_BEG=LPTS_BEG+NPX
        IF(LPTS_BEG.LT.NPTS) GOTO 250
        CALL GTTIME(TPOT2)
        IF(IPOT.EQ.2)ATIME=ATIME+TPOT2-TPOT
 600   CONTINUE
C
C MOVE TO COULOMB:
C
       IF(DEBUG)THEN
        DO LPTS=1,NMSH
         COULOMB(LPTS)=COULOMB(LPTS)+POTNL(LPTS)
        END DO
       ELSE
C
C THIS LOOP CAN BE VECTORIZED!
C
CDIR$ IVDEP
C
        DO LPTS=1,MPTS_CLS
         COULOMB(JPTR(LPTS))=COULOMB(JPTR(LPTS))+POTNL(LPTS)
        END DO
       END IF
c
c Is this used for anything anymore ?
c
       VLONG=0.0D0
       DO IANG=1,NANG
         VLONG=VLONG+DOMEGA(IANG)*ASYMPTOT(IANG)*RCUT/PI4
       END DO
       ACHRG=ACHRG+VLONG
c
C
       DO IANG=1,NANG
        Q(1,IANG)=ANG(1,IANG)*RCUT
        Q(2,IANG)=ANG(2,IANG)*RCUT
        Q(3,IANG)=ANG(3,IANG)*RCUT 
       END DO
c
c Q coordinates must be normalized for harmonics routine
c
      HNORM=SQRT(Q(1,1)**2+Q(2,1)**2+Q(3,1)**2)
      RECIPR(0)=HNORM
      HNORM=1.0D0/HNORM
      DO I=1,NANG
        Q(1,I)=HNORM*Q(1,I)
        Q(2,I)=HNORM*Q(2,I)
        Q(3,I)=HNORM*Q(3,I)
      END DO
      CALL HARMONICS(MPX,NANG,LM,Q,POLY,NPOLY)
C
C INTERPOLATE ONTO CLUSTER MESH:
C
       DO MOM=1,(LM+1)**2
        RMOMENT(MOM)=0.0D0
        DO IANG=1,NANG
         RMOMENT(MOM)=RMOMENT(MOM)
     &               +DOMEGA(IANG)*POLY(IANG,MOM)*ASYMPTOT(IANG)
        END DO
       END DO
C
      DO L=1,LM
        RECIPR(L)=RECIPR(0)*RECIPR(L-1)
      END DO
      MOM=0
      DO L=0,LM
       DO M=-L,L
        MOM=MOM+1
        RMOMENT(MOM)=RMOMENT(MOM)*RECIPR(L)
       END DO
      END DO
c
       NPTS=MPTS_FAR
       IF(NPTS.EQ.0)RETURN
       CALL GTTIME(TIME1)
       NTIMES=NPTS/MPX
       IF((NTIMES*MPX).LT.NPTS)NTIMES=NTIMES+1
       LPTS_BEG=0
       DO 650 II=1,NTIMES
       NPX=MIN(MPX,NPTS-LPTS_BEG)
       KPTS_MIN=IPTR(LPTS_BEG+1)  
       KPTS_MAX=IPTR(LPTS_BEG+NPX)
       DO 700 IPTS=1,NPX
        Q(1,IPTS)=RMSH(1,IPTR(LPTS_BEG+IPTS))+CGRAVITY(1)
        Q(2,IPTS)=RMSH(2,IPTR(LPTS_BEG+IPTS))+CGRAVITY(2)
        Q(3,IPTS)=RMSH(3,IPTR(LPTS_BEG+IPTS))+CGRAVITY(3)
 700   CONTINUE
c
c Q coordinates must be normalized for harmonics routine
c locate smallest r in the process
c
      RI1=0.0D0
      DO I=1,NPX
        HNORM=Q(1,I)**2+Q(2,I)**2+Q(3,I)**2
        HNORM=1.0D0/SQRT(HNORM)
        RECIP(I,0)=HNORM
        IF(RECIP(I,0).GT.RI1)RI1=RECIP(I,0)
        Q(1,I)=HNORM*Q(1,I)
        Q(2,I)=HNORM*Q(2,I)
        Q(3,I)=HNORM*Q(3,I)
      END DO
c
c test to see if we need to do any calculations
c
      DELTAV=1.0D-12
      DO L=0,LM
        RTEST(L)=RI1**(L+1)
      END DO
      LM_MAX=0
      MOM=0
      DO L=0,LM
        DO M=-L,L
          MOM=MOM+1
          SKIPTEST=(2*L+1)*RMOMENT(MOM)*RTEST(L)
          SKIPTEST=ABS(SKIPTEST)
          IF(SKIPTEST.LT.DELTAV)THEN
            SKIP(MOM)=.TRUE.
          ELSE
            SKIP(MOM)=.FALSE.
            LM_MAX=L
          END IF
        END DO
      END DO
c
c only get harmonics up to lm_max because the higher ones won't
c be needed
c
      CALL HARMONICS(MPX,NPX,LM_MAX,Q,POLY,NPOLY)
c 
c Generalized for arbitrary L 
c
c
      DO L=1,LM_MAX
        DO IPTS=1,NPX
          RECIP(IPTS,L)=RECIP(IPTS,0)*RECIP(IPTS,L-1)
        END DO
      END DO
      MOM=0
c
c generalized for arbitrary l

C$ OMP PARALLEL
C$ OMP DO
      DO IPTS=1,NPX
        POTLR(IPTS)=0.0D0
      END DO
C$ OMP ENDDO
C$ OMP END PARALLEL
       
      MOM=0

c DCP I am unrolling some of this for Mark (although it seems
c to have no effect on performance even on the IBMs)
c
c     do l=0,lm_max
c       do m=-l,l
c         mom=mom+1
c         if(.not.skip(mom))then
c           do ipts=1,npx
cc             potlr(ipts)=potlr(ipts)+rmoment(mom)*poly(ipts,mom)
c             potlr(ipts)=potlr(ipts)+rmoment(mom)*poly(ipts,mom)
c    &                    *recip(ipts,l)
c           end do
c         end if
c       end do
c     end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l=0
c
      IF(LM_MAX.GE.0)THEN
        L=0
        MOM=MOM+1
        IF(.NOT.SKIP(MOM))THEN
          DO IPTS=1,NPX
            POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
          END DO
        END IF
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l=1
c
      IF(LM_MAX.GE.1)THEN 
        L=1
c m=-1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=0
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l=2
c
      IF(LM_MAX.GE.2)THEN
        L=2
c m=-2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=0
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
      END IF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l=3
c
      IF(LM_MAX.GE.3)THEN
        L=3
c m=-3
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=0
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=3
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
      END IF
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l=4
c
      IF(LM_MAX.GE.4)THEN
        L=4
c m=-4
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-3
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=0
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=3
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=4
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
      END IF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l=5
c
      IF(LM_MAX.GE.5)THEN
        L=5
c m=-5
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-4
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-3
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=-1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=0
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=1
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=2
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=3
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=4
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
c m=5
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
      END IF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c l>5
c
      DO L=6,LM_MAX
        DO M=-L,L
          MOM=MOM+1
          IF(.NOT.SKIP(MOM))THEN
            DO IPTS=1,NPX
              POTLR(IPTS)=POTLR(IPTS)+RMOMENT(MOM)*POLY(IPTS,MOM)
     &                    *RECIP(IPTS,L)
            END DO
          END IF
        END DO
      END DO
c
       CALL GTTIME(WASTE1)
       IF(DEBUG)THEN
        DO IPTS=1,NPX
         DIFF=ABS(POTLR(IPTS)-POTNL(IPTR(LPTS_BEG+IPTS)))
         TDIFF=TDIFF+DIFF
         IF(DIFF.GT.DIFF_MAX)DIFF_MAX=DIFF
        END DO
       ELSE
        IF(KPTS_MAX-KPTS_MIN+1.EQ.NPX)THEN
         DO 720 IPTS=1,NPX
          COULOMB(KPTS_MIN+IPTS-1)=COULOMB(KPTS_MIN+IPTS-1)
     &                            +POTLR(IPTS)
  720    CONTINUE
        ELSE
C
CDIR$ IVDEP
C
         DO 725 IPTS=1,NPX
          COULOMB(IPTR(LPTS_BEG+IPTS))=COULOMB(IPTR(LPTS_BEG+IPTS))
     &                                +POTLR(IPTS)
  725    CONTINUE
        END IF
       END IF
       LPTS_BEG=LPTS_BEG+NPX
       CALL GTTIME(WASTE2)
  650  CONTINUE
       CALL GTTIME(TIME2)
       ATIME=ATIME+TIME2-TIME1
       IF(DEBUG)THEN
        DAVG=TDIFF/MPTS_FAR
        FLAG='    '
        IF(DIFF_MAX .GT. 0.0001)FLAG='<==='
        WRITE(73,900)MPTS_FAR,MPTS_CLS,DIFF_MAX,DAVG,FLAG
       END IF
  900  FORMAT(' ',2I7,2G15.6,A6) 
c
c I don't think we need this continue 950 do we ?
c
  950  CONTINUE
C
C DEALLLOCATE LOCAL ARRAYS
C
       DEALLOCATE(POTNL,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'POISSON1:ERROR DEALLOCATING POTNL'
       DEALLOCATE(IPTR,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'POISSON1:ERROR DEALLOCATING IPTR'
       DEALLOCATE(JPTR,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'POISSON1:ERROR DEALLOCATING JPTR'
       RETURN
       END
