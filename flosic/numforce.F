C UTEP Electronic Structure Lab (2020)
C
C ******************************************************************
C
      SUBROUTINE NUMFORCE
      use global_inputs,only : LIBXC1
      use mesh1,only : nmsh
      use common2,only : N_POS, NFNCT, FTOT, FRC1, FRC2, NSPN,
     &                   over1, ek, dftV, allH !<<<<
      use common3,only : NGRP
      use common5,only : NWF, NWFS
      use mpidat1,only : NPROC, NCALLED
      use SICMAT,only  : SIC, DERSIC !<<<<
      use FRM,only     : LFRM !<<<<
      use SICFLAG,only : LSICF
      use xtmp2a,only  : MIXIN,MIXINS,ISMGGA
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!     INCLUDE  'PARAMAS'  
      INCLUDE  'PARAMA2'  
      INTEGER :: I, I_POS, IERR, IFNCT, IID, IOFS, IOFS1, IOFS2, IPTS, 
     &           IX, LOFS, LPTS, LPTS_BEG, MPTS,
     &           NMAST, NMAX, NTID, NWFTOT,
     &           ISPN, IWF, JWF !<<<<
      REAL*8 :: SYMBOL 
      SAVE
      PARAMETER (NMAX=MPBLOCK)
C
      IF(.NOT. LSICF) THEN
       NWFTOT=0
       DO I=1,MXSPN
        NWFTOT=NWFTOT+NWFS(I)
       END DO
       IF (NWF.NE.NWFTOT) THEN
        write(6,*)'NUMFORCE: NWF AND NWFS ARE INCONSISTENT'
        CALL STOPIT
       END IF
      ELSE
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
c Added from Force0.0
c  NWFTOT = total of spin up and down occupied orbitals
c
       NWFTOT = LFRM(1) + LFRM(2)
c
c  print Lagrange multiplier matrix as test
c  Here LMM = <phi_i|VSIC_i|phi_j>
c
c      print *, 'numforce nwftot', nwftot
c      print *, 'SIC -- LAG MULT Matrix: up'
c      do iwf = 1, lfrm(1)
c       print 83, (sic(jwf,iwf,1),jwf=1,lfrm(1))
c      end do 
c      print *, 'SIC -- LAG MULT Matrix: down'
c      do iwf = 1, lfrm(2)
c       print 83, (sic(jwf,iwf,2),jwf=1,lfrm(2))
c      end do 
c
c
c zero test arrays
c
       do iwf = 1, nwftot
        do jwf = iwf, nwftot
         over1(iwf,jwf) = 0.0d0
         over1(jwf,iwf) = 0.0d0
         ek(iwf,jwf) = 0.0d0
         ek(jwf,iwf) = 0.0d0
         dftV(iwf,jwf)=0.0d0
         dftV(jwf,iwf)=0.0d0
         allH(iwf,jwf)=0.0d0
         allH(jwf,iwf)=0.0d0
        end do
       end do

c
c  LAGMULT_SLV adds the DFT part to the lagrange multiplier matrix.
c  After LAGMULT_SLV:  SIC(i,j,ispn) = <phi_i|HDFT + VSIC_i|phi_j>
c
c
c
c      KAI THIS SHOULD BE PARALLELIZED
c
       LPTS=0
   25  CONTINUE
        LPTS_BEG=LPTS
        MPTS=MIN(NMAX,NMSH-LPTS)
        CALL LAGMULT_SLV(MPTS,LPTS)

        LPTS = LPTS_BEG + MPTS
       if(LPTS.LT.NMSH) GOTO 25
C
c  done with LAGMULT
c
c
c  print test integrals
c
c      print *, 'in NUMFORCE:  test integrals'
c      print *, 'over 1'
c      do iwf = 1, nwftot
c       print 83, (over1(jwf,iwf),jwf=1,nwftot)
c      end do 
c      print  *, 'ek'
c      do iwf = 1, nwftot
c       print 83, (ek(jwf,iwf),jwf=1,nwftot)
c      end do 
c      print *, 'dftV'
c      do iwf = 1, nwftot 
c       print 83, (dftV(jwf,iwf),jwf=1,nwftot)
c      end do 
c      print *, 'allH'
c      do iwf = 1, nwftot
c       print 83, (allH(jwf,iwf),jwf=1,nwftot)
c      end do 
c      print *, 'SIC -- LAG MULT Matrix: up'
c      do iwf = 1, nwftot
c       print 83, (sic(jwf,iwf,1),jwf=1,nwftot)
c      end do 
c      if(LFRM(2) > 0) then
c       print *, 'SIC -- LAG MULT Matrix: down'
c       do iwf = 1, nwftot
c        print 83, (sic(jwf,iwf,2),jwf=1,nwftot)
c       end do
c      end if
83     format(1x,10(1x,f10.4))
c
c  SYMMETRIZE LM matrix 
c
       do ispn = 1,nspn
        do iwf = 1,lfrm(ispn)
         do jwf = iwf,lfrm(ispn)
          sic(jwf,iwf,ispn) = (sic(jwf,iwf,ispn)+
     &                   sic(iwf,jwf,ispn))/2.0d0
          sic(iwf,jwf,ispn) = sic(jwf,iwf,ispn)
         end do
        end do
       end do
c
c zero test arrays again -- to test frcslv
c
       do iwf = 1, nwftot
        do jwf = iwf, nwftot
         over1(iwf,jwf) = 0.0d0
         over1(jwf,iwf) = 0.0d0
         ek(iwf,jwf) = 0.0d0
         ek(jwf,iwf) = 0.0d0
         dftV(iwf,jwf)=0.0d0
         dftV(jwf,iwf)=0.0d0
         allH(iwf,jwf)=0.0d0
         allH(jwf,iwf)=0.0d0
        end do
       end do
      END IF
!>>>>>>>>>>>>>>>>>>>
c
#ifdef MPI
c
c MPI: send data to workers
c This will zero out FRC1 and FRC2 on workers
c
      IF (NPROC .GT. 0) THEN
       CALL SENDDATA(105)
       IF(LSICF) CALL SENDDATA(205) !Broadcast SIC matrix !<<<<
                 !You also want to broadcast LSICF
       NMAST= NMAX/(4*NPROC)
       IF(LIBXC1.OR.ISMGGA) NMAST = NMAX
        !This is needed to keep the correct libxc offset
      END IF
c
#endif
C
C LB OPEN AND READ LIBXCFRC FILE
C
      IF(LIBXC1.OR.ISMGGA)THEN
        ALLOCATE(MIXIN(4,NMSH*MXSPN),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'NUMFORCE:ERROR ALLOCATING MIXIN'
        ALLOCATE(MIXINS(4,MPBLOCK*MXSPN),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'NUMFORCE:ERROR ALLOCATING MIXINS'

        MIXIN(:,:) = 0.0d0

        open(12001,FILE='LIBXCFRC',FORM='UNFORMATTED',STATUS='OLD')
        rewind(12001)
        do IPTS=1,NMSH*NSPN
          read(12001) MIXIN(1,IPTS),MIXIN(2,IPTS),
     &                MIXIN(3,IPTS),MIXIN(4,IPTS)
        end do
        close(12001)
      ENDIF

C
C LOOP OVER ALL POINTS
C
      IOFS1=(NSPN-1)*NMSH   ! Spin down offset of the entire mesh
!      IOFS2=(NSPN-1)*MPTS   ! Spin down offset of chopped mesh
      LPTS=0
   35 CONTINUE
       LPTS_BEG=LPTS
       MPTS=MIN(NMAX,NMSH-LPTS)
       IOFS2=(NSPN-1)*MPTS   ! Spin down offset of chopped mesh
C
C LB TRANSFER FROM MIXIN ARRAY TO MIXINS
C
       IF(LIBXC1.OR.ISMGGA)THEN
         MIXINS(:,:)=0.0d0

         DO IPTS=1,MPTS
! GGA
!Spin up mixin elems
           MIXINS(1,IPTS)=MIXIN(1,LPTS+IPTS)
           MIXINS(2,IPTS)=MIXIN(2,LPTS+IPTS)
           MIXINS(3,IPTS)=MIXIN(3,LPTS+IPTS)
!Spin dn mixinelems
           MIXINS(1,IPTS+IOFS2)=MIXIN(1,LPTS+IPTS+IOFS1)
           MIXINS(2,IPTS+IOFS2)=MIXIN(2,LPTS+IPTS+IOFS1)
           MIXINS(3,IPTS+IOFS2)=MIXIN(3,LPTS+IPTS+IOFS1)

! meta-GGA
!Spin up mixin elems
           MIXINS(4,IPTS)=MIXIN(4,LPTS+IPTS)
!Spin dn mixinelems
           MIXINS(4,IPTS+IOFS2)=MIXIN(4,LPTS+IPTS+IOFS1)
         ENDDO
       ENDIF
c
#ifndef MPI
c
       IF(LSICF) THEN
#ifdef ATOMFORCE
        CALL FRCSLV0(MPTS,LPTS)
#else
        PRINT *,"SKIPPING ATOMIC FORCE"
        PRINT *,"RECOMPILE WTH ATOMFORCE TO USE IT"
#endif
       ELSE
        CALL FRCSLV(MPTS,LPTS)
       ENDIF
c
#else
c
c MPI: if all workers are busy, let manager do some work
c      but do not allow packets smaller than 10 points
c
       IF (NPROC .EQ. 0) THEN
        IF(LSICF) THEN
#ifdef ATOMFORCE
         CALL FRCSLV0(MPTS,LPTS)
#else
c        Moved the warning message to apotnl_dft
c        PRINT *,"SKIPPING ATOMIC FORCE"
c        PRINT *,"RECOMPILE WITH ATOMFORCE TO USE IT"
#endif
        ELSE
         CALL FRCSLV(MPTS,LPTS)
        END IF
       ELSE
        CALL CKWORKER(1,NTID)
        IF (NCALLED .NE. NPROC) THEN
         CALL PAMFPUL(1,MPTS,LPTS)
        ELSE
         MPTS= MIN(NMAST,NMSH-LPTS) !This statement will break libxc/mgga
                                    !offset unless NMAST == NMAX
         IF(LIBXC1.OR.ISMGGA) THEN
         !For meta-GGA, it is faster when the manager proc is free all time.
          MPTS= 0
          CALL CKWORKER(2,NTID)
         ELSE IF (MPTS .LE. 10) THEN
          MPTS= 0
          CALL CKWORKER(2,NTID)
         ELSE
          IF(LSICF) THEN
           CALL FRCSLV0(MPTS,LPTS)
          ELSE
           CALL FRCSLV(MPTS,LPTS)
          END IF
         ENDIF
        END IF
       END IF
c
#endif
c
       LPTS=LPTS_BEG+MPTS
       IF(LPTS.LT.NMSH) GOTO 35
      CONTINUE
c
c MPI: collect data
c
#ifdef MPI
c
      CALL PAMFPUL(2,0,0)
c
#endif
C
C DEALLOCATE LOCAL ARRAYS
C
      IF(LIBXC1.OR.ISMGGA)THEN
        DEALLOCATE(MIXIN,STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'NUMFORCE:ERROR DEALLOCATING MIXIN'
        DEALLOCATE(MIXINS,STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'NUMFORCE:ERROR DEALLOCATING MIXINS'
      ENDIF
c
C FORCES ARE TOO BIG BY A FACTOR OF N(G)
C
c Debug print for Force 0.0
c     print *, 'numforce nfnct npos', nfnct, 
c    &     (n_pos(ifnct),ifnct=1,nfnct)
      IID=0
      DO IFNCT=1,NFNCT
       DO I_POS=1,N_POS(IFNCT)
        IID=IID+1
        DO IX=1,3
         FRC1(IX,IID)=FRC1(IX,IID)/NGRP
         FRC2(IX,IID)=FRC2(IX,IID)/NGRP
        END DO
c       print *, 'frc1', (frc1(ix,iid),ix=1,3)
c       print *, 'frc2', (frc2(ix,iid),ix=1,3)
       END DO
      END DO
      RETURN
      END
