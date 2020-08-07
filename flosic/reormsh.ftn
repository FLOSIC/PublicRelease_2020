C UTEP Electronic Structure Lab (2020)
C
       SUBROUTINE REORMSH
C
C WRITTEN BY MARK R PEDERSON 3-FEB 1998
C
       use mesh1,only : wmsh,rmsh,nmsh
       use debug1
       use common3,only : RMAT, NGRP
       use common6,only : TSPH, LIMSPH, NSPHERES
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:57 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MXBOX, I, IBEG, IBOX, IERR, IGP, INBOX, IPTS, IX, IY,
     & IZ, J, K, KAVG, KMAX, KMIN, KPTS, L1, L2, LBOX, LPTS, MAXBOX,
     & MBOX, MM, N7, NBOX, NX
       REAL*8 :: SYMBOL , BX, DX, EPS, ERRX, FACT, PAVG, PMAX, PMIN,
     & RAVG, RMAX, RMIN, RR, S, SAVG, SIDE, SIDMIN, SMAX, SMIN, SN,
     & SX1, SX2, SY1, SY2, SZ1, SZ2, VAVG, VMAX, VMIN, VOL, WASP
       SAVE
       PARAMETER (MXBOX=MAX_PTS/20)
       DIMENSION SIDE(2,3),S(3),SN(3,3),NX(3)
       REAL*8 :: BOX(2,3,MXBOX)
       INTEGER :: MPTS(MXBOX),IXNEW(MXBOX)
       INTEGER,ALLOCATABLE :: KBOX(:)
       DATA EPS/1.0D-5/
       DATA MAXBOX/400/ 
       DATA SIDMIN/1.0D0/
C
C FIND RANGE OF POINTS TO DETERMINE ORIGINAL BOX
C
       ALLOCATE(KBOX(MAX_PTS),STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'REORMSH:ERROR ALLOCATING KBOX'

       DO I=1,3
        S(I)=0.0D0
       END DO
       DO IPTS=1,NMSH
        DO I=1,3
         IF (ABS(RMSH(I,IPTS)).GT.S(I)) S(I)=ABS(RMSH(I,IPTS))
        END DO
       END DO
       DO I=1,3
        SIDE(1,I)= -(S(I)+0.1D0)
        SIDE(2,I)=   S(I)+0.1D0
       END DO
C
C TRY TO REDUCE BOX SIZE BY CHECKING FOR X/Y/Z REFLECTIONS 
C
       DO IGP=2,NGRP
        ERRX=ABS(RMAT(1,2,IGP))+ABS(RMAT(1,3,IGP))
     &      +ABS(RMAT(2,1,IGP))+ABS(RMAT(2,3,IGP))
     &      +ABS(RMAT(3,1,IGP))+ABS(RMAT(3,2,IGP))
        IF (ERRX.LT.EPS) THEN
         ERRX=ABS(RMAT(1,1,IGP)+1.0D0)+
     &        ABS(RMAT(2,2,IGP)-1.0D0)+
     &        ABS(RMAT(3,3,IGP)-1.0D0)
         IF (ERRX.LT.EPS) SIDE(1,1)=0.0D0
         ERRX=ABS(RMAT(1,1,IGP)-1.0D0)+
     &        ABS(RMAT(2,2,IGP)+1.0D0)+
     &        ABS(RMAT(3,3,IGP)-1.0D0)
         IF (ERRX.LT.EPS) SIDE(1,2)=0.0D0
         ERRX=ABS(RMAT(1,1,IGP)-1.0D0)+
     &        ABS(RMAT(2,2,IGP)-1.0D0)+
     &        ABS(RMAT(3,3,IGP)+1.0D0)
         IF (ERRX.LT.EPS) SIDE(1,3)=0.0D0
        END IF
       END DO
C
C MOVE ALL POINTS INTO SMALLER BOX
C
       DO IPTS=1,NMSH
        DO IGP=1,NGRP
         DO J=1,3
          S(J)=0.0D0
          DO K=1,3
           S(J)=S(J)+RMAT(J,K,IGP)*RMSH(K,IPTS)
          END DO
         END DO
         INBOX=0
         DO I=1,3
          IF (S(I).LT.SIDE(1,I)) INBOX=1
          IF (S(I).GT.SIDE(2,I)) INBOX=1
         END DO
         IF (INBOX.EQ.0) GOTO 10
        END DO
        write(6,*)'REORMSH: PROBLEM 1'
        CALL STOPIT
 10     CONTINUE
        DO K=1,3
         RMSH(K,IPTS)=S(K)
        END DO
       END DO
C
C SETUP FOR FIRST BOX
C
       DO I=1,3
        DO J=1,2
         BOX(J,I,1)=SIDE(J,I)
        END DO
       END DO
       DO I=1,NMSH
        KBOX(I)=1
       END DO
       MPTS(1)=NMSH
C
C BEGIN SUBDIVISON
C                
       NBOX=1
   50  CONTINUE
!YY SIC insert begins
C
C GET RID OF BOXES DEVOID OF POINTS
C
        LPTS=0
        LBOX=0
        DO IBOX=1,NBOX
         LPTS=LPTS+MPTS(IBOX)
         IXNEW(IBOX)=0
         IF (MPTS(IBOX).NE.0) THEN
          LBOX=LBOX+1
          IXNEW(IBOX)=LBOX
          MPTS(LBOX)=MPTS(IBOX)
          DO I=1,3
           DO J=1,2
            BOX(J,I,LBOX)=BOX(J,I,IBOX)
           END DO
          END DO
         END IF
        END DO
        IF (LBOX .NE. NBOX) THEN
         DO I=1,NMSH
          KBOX(I)= IXNEW(KBOX(I))
         END DO
        END IF
        NBOX=LBOX
C
C FIND LARGEST BOX THAT SHOULD BE SLICED
C
        MBOX=NBOX+1
        DO IBOX=NBOX,1,-1
         SMAX=0.0D0
         DO I=1,3
          SMAX=MAX(SMAX,BOX(2,I,IBOX)-BOX(1,I,IBOX))
         END DO
         IF (SMAX.GT.SIDMIN) THEN
          IF (MPTS(IBOX).GT.MAXBOX) MBOX=IBOX
         END IF
        END DO
C
C DETERMINE HOW MANY SLICES
C
        IF (MBOX.LE.NBOX) THEN
         SAVG=0.0D0
         DO I=1,3
          S(I)=BOX(2,I,MBOX)-BOX(1,I,MBOX)
          SAVG=SAVG+S(I)
         END DO
         SAVG=SAVG/3.0D0
         DO I=1,3
          S(I)=S(I)/SAVG
         END DO
         DO I=1,3
          NX(I)=1
          IF (S(I).GT.0.8) NX(I)=2
         END DO
         N7=NX(1)*NX(2)*NX(3)-1
         IF (NBOX+N7 .GT. MXBOX) THEN
          PRINT *,'REORMSH: MXBOX MUST BE AT LEAST: ',NBOX+N7
          CALL STOPIT
         END IF
C
C MAKE ROOM FOR NEW BOXES
C
         DO IBOX=NBOX,MBOX+1,-1
          DO I=1,3
           DO J=1,2
            BOX(J,I,IBOX+N7)=BOX(J,I,IBOX)
           END DO
          END DO
          MPTS(IBOX+N7)=MPTS(IBOX)
         END DO
         DO I=1,NMSH
          IF (KBOX(I).GT.MBOX) KBOX(I)=KBOX(I)+N7
         END DO
         NBOX=NBOX+N7
C
C REPARTITION BOX
C
         DO I=1,3
          BX=BOX(1,I,MBOX)
          DX=(BOX(2,I,MBOX)-BOX(1,I,MBOX))/NX(I)
          DO J=1,NX(I)+1
           SN(J,I)=BX
           BX=BX+DX
          END DO
         END DO
C
C CREATE NEW BOXES
C
         MBOX=MBOX-1
         DO IX=1,NX(1)
          SX1=SN(IX  ,1)
          SX2=SN(IX+1,1)
          DO IY=1,NX(2)
           SY1=SN(IY  ,2)
           SY2=SN(IY+1,2)
           DO IZ=1,NX(3)
            SZ1=SN(IZ  ,3)
            SZ2=SN(IZ+1,3)
            MBOX=MBOX+1
            BOX(1,1,MBOX)=SX1
            BOX(2,1,MBOX)=SX2
            BOX(1,2,MBOX)=SY1
            BOX(2,2,MBOX)=SY2
            BOX(1,3,MBOX)=SZ1
            BOX(2,3,MBOX)=SZ2
           END DO
          END DO
         END DO
C
C NOW RELOCATE THE MESH POINTS
C
         DO IBOX=MBOX,MBOX-N7,-1
          MPTS(IBOX)=0
         END DO
         DO IPTS=1,NMSH
          IF (KBOX(IPTS).EQ.MBOX-N7) THEN
           DO IBOX=MBOX,MBOX-N7,-1
C
C TRY TO MOVE THIS POINT TO IBOX
C
            DO IGP=1,NGRP
             DO J=1,3
              S(J)=0.0D0
              DO K=1,3
               S(J)=S(J)+RMAT(J,K,IGP)*RMSH(K,IPTS)
              END DO
             END DO
             INBOX=0
             DO I=1,3
              IF (S(I).LT.BOX(1,I,IBOX)) INBOX=1
              IF (S(I).GT.BOX(2,I,IBOX)) INBOX=1
             END DO
             IF (INBOX.EQ.0) THEN
              DO J=1,3
               RMSH(J,IPTS)=S(J)
              END DO
              MPTS(IBOX)=MPTS(IBOX)+1
              KBOX(IPTS)=IBOX
              GOTO 20
             END IF
            END DO
           END DO
          END IF
 20       CONTINUE
         END DO
         GOTO 50
        END IF
       CONTINUE
!YY SIC insert end
C
C MOVE POINTS TOGETHER
C
       LPTS=0     
       DO IBOX=1,NBOX
        IBEG=LPTS+1
        KPTS=0
        DO IPTS=IBEG,NMSH
         IF (KBOX(IPTS).EQ.IBOX) THEN
          LPTS=LPTS+1
          KPTS=KPTS+1
          DO J=1,3
           RR=RMSH(J,IPTS)
           RMSH(J,IPTS)=RMSH(J,LPTS)
           RMSH(J,LPTS)=RR
          END DO
          RR=WMSH(IPTS)
          WMSH(IPTS)=WMSH(LPTS)
          WMSH(LPTS)=RR
          MM=KBOX(IPTS)
          KBOX(IPTS)=KBOX(LPTS)
          KBOX(LPTS)=MM
         END IF
        END DO
       END DO
C
C CHECK TO SEE IF POINTS ARE ORDERED
C
C      LPTS=0
C      IERR=0
C      DO IBOX=1,NBOX
C       DO IPTS=1,MPTS(IBOX)
C        LPTS=LPTS+1
C        IF (KBOX(LPTS).NE.IBOX) IERR=IERR+1
C       END DO
C      END DO
C      IF (IERR.NE.0) THEN
C       write(6,*)'REORMSH: ARRANGEMENT ERROR:',IERR
C       CALL STOPIT
C      END IF
C
C PRINT OUT INFORMATION ABOUT SPHERES
C
       IF (NBOX.GT.MX_SPH) THEN
        write(6,*)'REORMSH: MX_SPH MUST BE ATLEAST:',NBOX
        CALL STOPIT
       END IF
       RMIN=1.0D30
       RMAX=0.0D0
       RAVG=0.0D0
       VMIN=1.0D30
       VMAX=0.0D0
       VAVG=0.0D0
       PMIN=1.0D30
       PMAX=0.0D0
       PAVG=0.0D0
       KMIN=NMSH
       KMAX=0
       WASP=1.0D0
       FACT=4.0D0*ATAN(1.0D0)*4.0D0/3.0D0
       IF (DEBUG) write(6,*)NBOX,' TOTAL NUMBER OF SPHERES'
C
       LPTS=1
       DO IBOX=1,NBOX
        RR=0.0D0
        DO I=1,3
         S(I)=(BOX(1,I,IBOX)+BOX(2,I,IBOX))/2.0D0
         RR=RR+(S(I)-BOX(1,I,IBOX))**2
        END DO
        RR=SQRT(RR)
C
        VOL=FACT*RR**3 
        IF (RR.GT.RMAX) RMAX=RR
        IF (RR.LT.RMIN) RMIN=RR
        RAVG=RAVG+RR 
        IF (VOL.GT.VMAX) VMAX=VOL
        IF (VOL.LT.VMIN) VMIN=VOL
        VAVG=VAVG+VOL 
        IF (MPTS(IBOX).GT.KMAX) KMAX=MPTS(IBOX)
        IF (MPTS(IBOX).LT.KMIN) KMIN=MPTS(IBOX)
        VOL=MPTS(IBOX)/VOL
        IF (VOL.GT.PMAX) PMAX=VOL
        IF (VOL.LT.PMIN) PMIN=VOL
        PAVG=PAVG+VOL
C
        TSPH(1,IBOX)=S(1)
        TSPH(2,IBOX)=S(2)
        TSPH(3,IBOX)=S(3)
        TSPH(4,IBOX)=RR
C
C CONSTRUCT ASPECT RATIOS
C
        SMAX=0.0D0
        SMIN=1.0D30
        DO I=1,3
         S(I)=(BOX(2,I,IBOX)-BOX(1,I,IBOX))
         IF (S(I).GT.SMAX) SMAX=S(I)
         IF (S(I).LT.SMIN) SMIN=S(I)
        END DO
        WASP=MAX(WASP,SMAX/SMIN)
        DO I=1,3
         S(I)=S(I)/SMAX
        END DO
        L1=LPTS
        L2=LPTS+MPTS(IBOX)
        LIMSPH(1,IBOX)=L1
        LIMSPH(2,IBOX)=L2-1
        LPTS=LPTS+MPTS(IBOX)
        IF (DEBUG) PRINT 1010,IBOX,L1,L2,
     &                   (TSPH(J,IBOX),J=1,4),(S(I),I=1,3)
 1010   FORMAT(' ',3I7,4F14.9,3F10.2)
       END DO
C
       RAVG=RAVG/NBOX
       VAVG=VAVG/NBOX
       PAVG=PAVG/NBOX
       KAVG=NMSH/NBOX
       PRINT '(A,I5)','TOTAL NUMBER OF SPHERES: ',NBOX
       PRINT '(2A)','SPHERE INFORMATION: RADIUS, VOLUME, ',
     &              'NUMBER OF POINTS, POINT DENSITY'
       PRINT 1020,'MINIMUM:',RMIN,VMIN,KMIN,PMIN
       PRINT 1020,'MAXIMUM:',RMAX,VMAX,KMAX,PMAX
       PRINT 1020,'AVERAGE:',RAVG,VAVG,KAVG,PAVG
 1020  FORMAT(A,1X,F10.3,1X,F10.3,1X,I8,1X,F10.3)
       PRINT '(A,F10.3,A)','WORST ASPECT RATIO: ',WASP,' : 1'
C
C CHECK TO SEE IF MESH POINTS ARE IN ADVERTISED PLACES I
C
       IERR=0
       DO IPTS=1,NMSH
        LBOX=KBOX(IPTS)
        INBOX=0
        DO I=1,3
         IF (RMSH(I,IPTS).LT.BOX(1,I,LBOX)) INBOX=1
         IF (RMSH(I,IPTS).GT.BOX(2,I,LBOX)) INBOX=1
        END DO
        IF (INBOX.EQ.1) IERR=IERR+1
       END DO
       IF (IERR.NE.0) THEN
        write(6,*)'REORMSH: BOX ASSIGNMENT ERROR:',IERR  
        CALL STOPIT
       END IF
C
C CHECK TO SEE IF MESH POINTS ARE IN ADVERTISED PLACES II
C
       LPTS=0
       IERR=0
       NSPHERES=NBOX
       LPTS=0
       DO IBOX=1,NBOX
        DO IPTS=LIMSPH(1,IBOX),LIMSPH(2,IBOX)
         LPTS=LPTS+1
         S(1)=RMSH(1,IPTS)-TSPH(1,IBOX)
         S(2)=RMSH(2,IPTS)-TSPH(2,IBOX)
         S(3)=RMSH(3,IPTS)-TSPH(3,IBOX)
         RR=S(1)*S(1)+S(2)*S(2)+S(3)*S(3)
         RR=SQRT(RR)
         IF (RR.GT.TSPH(4,IBOX)) IERR=IERR+1
        END DO
       END DO
       IF (IERR.NE.0) THEN
        write(6,*)'REORMSH: SPHERE ASSIGNMENT ERROR:',IERR  
        CALL STOPIT
       END IF
       IF (LPTS.NE.NMSH) THEN
        write(6,*)'REORMSH: LPTS AND NMSH DIFFER'   
        CALL STOPIT
       END IF

       DEALLOCATE(KBOX,STAT=IERR)
       IF(IERR.NE.0) WRITE(6,*)'REORMSH:ERROR DEALLOCATING KBOX'

       RETURN
      END
