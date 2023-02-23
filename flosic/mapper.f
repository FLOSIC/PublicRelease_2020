C
C ***************************************************************
c   KAJ 9-2-2022
c   MAPPER DETERMINES HOW TO OBTAIN THE VALUE OF A SYMMETRY
C   EQUIVALENT FLO FROM THE VALUE OF THE CORRESPONDING 
C   SYMMETRY INEQUIVALENT FLO ON UNSYMMETRIZED MESH
C
c    VARIABLES PASSED TO MAPPER:
c      MEQV -- number of ineq FLOs for the current FLO
c      IFND -- array of FLO indices
c    VARIABLES PASSED BACK FROM MAPPER:
c      MGRP -- number of wedges needed to get info
c      MAP  -- ARRAY WITH WEDGE INDICES
       SUBROUTINE MAPPER(MGRP,idim,MEQV,IFND,MAP)
       INCLUDE 'PARAMS'
       INCLUDE 'commons.inc'
       DIMENSION IFND(MX_GRP)
c
c  KAJ  FLO info in common blocks
c
       dimension ibeg(3), iend(3)
       common/tmp4/pot(max_pts)
       COMMON/FORRx/NFLO,KSPX,TMAT(NDH,NDH,2)
c      DIMENSION OG(40*MX_CNT,MAX_OCC),MAP(MX_GRP,1000)    !doesn't work 
       DIMENSION OG(MAX_PTS,MAX_OCC)
       DIMENSION MAP(MX_GRP,1)    
c      DIMENSION OG(MAX_PTS,MAX_OCC),MAP(MX_GRP,1000)    ! works
c      DIMENSION OG(40*MX_CNT,mx_grp),MAP(MX_GRP,MAX_OCC) ! should work

c  MAX_OCC could be smaller -- MX_GRP, but need to track FLO indices
c
c  KAJ  dimension arrays needed for map grid
c
c  
       dimension RAD(12)
       dimension RANG(3,10)
       DATA RAD/0.05d0,0.1d0,0.2d0,0.5d0,1.0d0,2.0d0,3.0d0,
     &          4.0d0,5.0d0,6.0d0,7.0d0,8.0d0/
C  randomize angular points
       DATA RANG/0.57735d0, 0.57735d0, 0.57735d0,
     &          -0.57735d0, 0.57735d0, 0.57735d0,
     &          -0.67735d0, 0.27735d0,-0.37735d0,
     &           0.77735d0,-0.57735d0,-0.17735d0,
     &           1.00d0   , 0.00d0   , 0.00d0   ,
     &          -1.00d0   , 0.00d0   , 0.00d0   ,
     &           0.00d0   , 1.00d0   , 0.00d0   ,
     &           0.00d0   ,-1.00d0   , 0.00d0   ,
     &           0.00d0   , 0.00d0   , 1.00d0   ,
     &           0.00d0   , 0.00d0   ,-1.00d0   /
c
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA ZED/1.0D-30/
       GAUSS_CUT=1.0D30
       IF(40*MX_CNT.gt.MAX_PTS) then
          print *, 'MAPPER 40*MX_CNT bigger than MAX_PTS'
          call stopit
       END IF
c
c  LMSH = # of unsymmetrized of points for MAPPER
c 
       NF1=1          
       NF2=MEQV        
       MGRP=0
C 
C  KAJ  Create mapping grid 
c  
c  Loop over symmetry independent atoms
c
         NMSH = 0
       ISHELLA = 0
       DO IFNCT = 1,NFNCT
       DO IPOS = 1,N_POS(IFNCT) 
         ISHELLA = ISHELLA + 1
                  NRAD = 10
                  NANG = 4  !10
                   do irad = 1,nrad
                    radx = RAD(irad)
                    do iang = 1,nang
                      nmsh = nmsh + 1
                      do j = 1,3
                        rmsh(j,nmsh) = ridt(j,ishella) + 
     &                                     rang(j,iang)*rad(irad) 
                      end do
                      print *, 'nmsh, rmsh', nmsh, (rmsh(j,nmsh),j=1,3)
                    end do
                   end do
       END DO
       END DO
       LMSH=NMSH ! NUMBER OF POINTS IN THE FUNDAMENTAL WEDGE

c
C  Unroll new grid using group operations
c
                  NMSH_SAV=NMSH
                  MMSH=NMSH
                  DO IMSH=1,NMSH
                  WMSH(IMSH)=WMSH(IMSH)/NGRP
                  END DO
                  DO IGP=2,NGRP 
                  DO IMSH=1,NMSH
                      MMSH=MMSH+1
                      WMSH(MMSH)=WMSH(IMSH)
                    DO J=1,3
                    RMSH(J,MMSH)=0.0D0
                    DO L=1,3
                    RMSH(J,MMSH)=RMSH(J,MMSH)+RMAT(J,L,IGP)*RMSH(L,IMSH)
                    END DO
                    END DO
                  END DO
                  END DO
                  NMSH=MMSH !TOTAL NUMBER OF POINTS
                  PRINT*,'POINTS BEFORE FLONASE:',NMSH
c
c  OBTAIN FLOS ON THE MAPPER GRID
C
       DO JFLO=NF1,NF2 
          IFLO=IFND(JFLO)
          NFLO=-IFLO
       POT=0.0D0
       PRINT*,'NFLO:',NFLO,KSPX,IFLO,NF1,NF2
       CALL FLONASE(TIME)
          TDD=0.0D0
          DO IMSH=1,MMSH
          TDD=TDD+ABS(POT(IMSH))
          OG(IMSH,IFLO)=POT(IMSH)
          END DO
          IF(ABS(TDD).LT.1.0D-10)THEN
              print *, 'TDD zero in MAPPER'
              STOP
          END IF
       END DO
c
c  Perform mapping function -  FLO indices stored in IFND
c
          MAP=NGRP+1
          DO JFLO=NF1,NF2
             IFLO=IFND(JFLO)
                  DO IGP=1,NGRP
                  MAP(IGP,IFLO)=NGRP+1
                  END DO
           DO JGP=1,NGRP              
           DO IGP=1,NGRP
             ERROR=0.0D0
             DO IMSH=1,LMSH
             ERROR=ERROR+ABS(OG(IMSH+(IGP-1)*LMSH,IFND(NF1))  
     &                      -OG(IMSH+(JGP-1)*LMSH,IFLO))
             END DO 
             IF(ERROR.LT.1.0D-5)THEN
             MAP(JGP,IFLO)=MIN(IGP,MAP(JGP,IFLO))
             MGRP=MAX(MAP(JGP,IFLO),MGRP)
             PRINT 610,LMSH,IFLO,IGP,JGP,ERROR
             END IF
           END DO  !IGP
           END DO  !JGP
c         PRINT 615,IFLO,(MAP(JGP,IFLO),JGP=1,NGRP),IFND(NF1)
 610         FORMAT(4I5,G15.6)
 615         FORMAT('MAP IFLO:',I4,48I3)
c
c  KAJ 6/30/22  CHECK OF MAP -- are all inequivalent FLO's mapped?
c
             DO IGP = 1,NGRP
              IF(MAP(IGP,IFLO).eq.NGRP+1) then
                print *, 'MAPPER FAILED FOR IFLO,IGP',IFLO,IGP
                CALL STOPIT
              END IF
             END DO  !IGP
          END DO     !FLOs
c
c  KAJ 7/4/22  RESTORE VARIATIONAL MESH 
c
           PRINT*,'NMSH:',NMSH
           OPEN(99,FILE='VMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
           REWIND(99)
           READ(99) NMSH,JCALC
           READ(99)((RMSH(J,I),J=1,3),I=1,NMSH)
           READ(99)(WMSH(I),I=1,NMSH)
           CLOSE(99)
         print *, 'in mapper nmsh', nmsh

       RETURN
       END
