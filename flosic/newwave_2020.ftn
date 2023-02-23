       SUBROUTINE WF_CONTROL(LUN,MODE)
       IF(MODE.EQ.0)THEN
       WRITE(LUN,*)'CONTROLS:'
       WRITE(LUN,*)'NO-FLIP -- Do Nothing' 
       WRITE(LUN,*)'ZERO    -- Set Eigenvector to ZERO'
       WRITE(LUN,*)'FLIP    -- Move Eigenvector to opposite spin'
       WRITE(LUN,*)'POLARIZE-- Force spin-polarized calculation'
       WRITE(LUN,*)'PROJECT -- move state up by 1 au           '
       END IF
       END
C      SUBROUTINE NEWWAVE(NITER,TRACE)
C      [ ] NON-EXECUTABLE LINE
C      CALL NEWAVE_2020(NITER,TRACE,ISUCCEED)
C      IF(ISUCCEEED.EQ.1)RETURN
C
C
       SUBROUTINE NEWWAVE_2020(NITER,TRACE,ISUCCEED)
C
C WRITTEN BY MARK R PEDERSON (1986-1989)
c
       use global_inputs,only : mixing1
       use for_diag1
       use hstor1, only : hstor
       use mixpot1,only : POTIN,POT=>POTOUT
       use debug1
       use common2,only : E_UP, E_DN, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : ISTSCF, PSI_COEF, OCCUPANCY,
     &   N_OCC, NWF, NWFS, EFERMI, EVLOCC 
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT, LDMREP

!      INCLUDE 'PARAMS'
       INCLUDE 'PARAMA2'
!      INCLUDE 'commons.inc'
       INTEGER :: I,I_REP,i1,IB,IBAS,IERR,i2,ii,ind,IPTS,IR,IREC,IREP,
     &            IS,ISAV,ISUCCEED,ITER,itimes,itry,itt,
     &            jb,JBAS,jnd,JSPN,jwf,jwi,
     &            k,KBAS,KNT,LBAS,lst,MAX_TOT,MDRN,MREC,MSPN,MXX,
     &            MXXX,N_VIRT,NBAS,NHTOT,NITER,NREC,NSAV,NSPTS,
     &            ntimes

       REAL*8 :: AVG,CUTOCC,DSINGV,ef,electrons,emax,emin,etarget,
     &           EVALSAV,EVLO,GAP,OCCV,
     &           spfact,TEMP,TIME1,TIMEA,TIMEB,TRACE

       SAVE
       PARAMETER (MAX_TOT=NDH*MAX_REP)
       LOGICAL EXIST,FERMISTAT
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM,sichammix
       CHARACTER*80 FLINE
       CHARACTER*12 EVALSTR
       CHARACTER*7 NAMES
       character*1 ud(2)
       character*3 ch
       character*10 control_wave
       character*25 comment
       DIMENSION NAMES(3)
C       DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION EVALSAV(MAX_TOT*MXSPN)
C       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
C     &  NSAV(MAX_REP,MXSPN)
       DIMENSION NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
C       DIMENSION DIAG(NDH,MAX_REP)
C       DIMENSION NTEMP(MAX_TOT*MXSPN)
       INTEGER,ALLOCATABLE :: NDEG(:),INDREP(:),NTEMP(:)
       REAL*8, ALLOCATABLE :: OCCTMP(:),DIAG(:,:)
       DIMENSION EVLO(NDH,MAX_REP,MXSPN),OCCV(NDH,MAX_REP,MXSPN)
C       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN) 
       DIMENSION KNT(3,MAX_TOT*MXSPN)          
C
!      COMMON/SICSHIFT/SICSCI(NDH_TOT,2)
! additional
       INTEGER :: J,KSPN
       REAL*8 :: SICSCI(NDH_TOT,2)
C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
       DATA NAMES/'BROYDEH','KBROY3','KBROY4'/
       POTIN =0.0D0
       POT   =0.0D0

!      hammixing option H_DFA only or H_DFA+SIC
       sichammix=.true. 
!reading SICSCI  !if sichammix=.true., add it in the sichammixdrv
       CALL CHECK_INPUTS
       if(sichammix .and. (MIXING1.EQ.2)) then
        continue !do nothing
       else
        SICSCI=0.0D0
        INQUIRE(FILE='SICSCI',EXIST=EXIST)
        IF(EXIST)THEN
         OPEN(50,FILE='SICSCI',FORM='UNFORMATTED')
         DO KSPN=1,NSPN
          READ(50)(SICSCI(J,KSPN),J=1,NDH_TOT)
         END DO
         CLOSE(50)
        END IF
       control_wave='NO-FLIP'
        ISUCCEED=1
        IF(MAX_VIRT_PER_SYM.NE.NDH)THEN
        PRINT*,'IF USING NEWWAVE, MAX_VIRT_PER_SYM MUST BE:',NDH
        ISUCCEED=0
        RETURN
        END IF 
       end if
       OPEN(97,FILE='TMPTRE')
       READ(97,*,END=10)TEMP
 10    CONTINUE
       REWIND(97)
       WRITE(97,*)TEMP,' TEMPERATURE'
       CLOSE(97)
       IF (N_REP.GT.MAX_REP) THEN
        PRINT *,'NEWWAVE: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
       END IF
       IF (NSPN.GT.MXSPN) THEN
        PRINT *,'NEWWAVE: MXSPN MUST BE AT LEAST: ',NSPN
        CALL STOPIT
       END IF
       TRACE=0.0D0
C
C REMOVE OLDER EVALXXX FILES IF NITER=1
C
       IF (NITER.EQ.1) THEN
        CALL SYSTEM('rm GOTOLD')
        ITER=0
   15    ITER=ITER+1
         WRITE(EVALSTR,'(A,I3.3)')'EVAL',ITER 
         INQUIRE(FILE=EVALSTR,EXIST=EXIST)
         IF (EXIST) THEN
          OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='OLD')
          CLOSE(98,STATUS='DELETE')
          GOTO 15
         END IF  
        CONTINUE
       END IF 
       FLINE='NORM ... CHOICES ARE NORM, OCCU OR FIXM'  
       OPEN(47,FILE='CONTROL_WFOUT',FORM='FORMATTED')
       OPEN(97,FILE='EVALUES')
       READ(97,'(A)',END=60,ERR=60)FLINE
 60    CONTINUE
                                                  MDRN=1
       if(FLINE(1:4).EQ.'FIXM'.OR.FLINE.EQ.'fixm')MDRN=2
       if(FLINE(1:4).EQ.'OCCU'.OR.FLINE.EQ.'occu')MDRN=3
       if(niter.eq.1.and.ISTSCF.eq.3)then
        if(fline(1:4).ne.'OCCU')MDRN=2
       end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if(mdrn.eq.3)THEN
        occv=0.0d0
        do ispn=1,nspn
         do irep=1,n_rep
          read(97,*)n_occ(irep,ispn)
          read(97,*)(occv(j,irep,ispn),j=1,n_occ(irep,ispn))
         end do
        end do
       end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C open FLOSIC scissor shift file:         
C      CALL SCISSORSHIFT
C      INQUIRE(FILE='SCISSORP',EXIST=EXIST)
C      IF(EXIST)THEN
C       OPEN(57,FILE='SCISSORP',FORM='UNFORMATTED')
C       REWIND(57)
C      END IF 
       ISPN=1
C      IF(ISPN.EQ.1)CALL OVERLAP(2)
       CALL GTTIME(TIMEA)
C
C Dynamic allocation
       NHTOT=0
       DO ISPN=1,NSPN
          DO I_REP=1,N_REP
           NHTOT=NHTOT+NS_TOT(I_REP)
          END DO
       END DO
C       write(6,*)'newwave:NHTOT=',NHTOT,NSPN,N_REP,NS_TOT(1) 

       NHTOT=NHTOT*MXSPN
       ALLOCATE(NDEG(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating NDEG'
       ENDIF
       ALLOCATE(INDREP(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating INDREP'
       ENDIF
       ALLOCATE(NTEMP(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating NTEMP'
       ENDIF
       ALLOCATE(OCCTMP(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating OCCTMP'
       ENDIF

C       write(6,*)'newwave:NHTOT=',NHTOT
C       write(6,*)'newwave:MAX_TOT*MXSPN=',MAX_TOT*MXSPN

C       DO IR=1,MAX_TOT*MXSPN
       DO IR=1,NHTOT
           OCCTMP(IR)=0.0D0
       END DO
C
CCCCCCCCCCCC NEW PART CCCCCCCCCCCCCCCCC
       CALL OVERLAP(1)
       
       !Below is done inside of Hammixdriver
       !IREC=0
       !ALLOCATE(DIAG(NDH,MAX_REP))
       !DO IREP=1,N_REP 
       ! DO I=1,NS_TOT(IREP)
       !  DO J=I,NS_TOT(IREP)
       !   IREC=IREC+1
       !   IF(J.EQ.I)THEN
       !    DIAG(I,IREP)=1.0D0/SQRT(HSTOR(IREC,1))
       !    IF(DIAG(I,IREP).NE.DIAG(I,IREP))THEN
       !     PRINT*,'DIAG IS NOT A NUMBER'
       !     CALL STOPIT
       !    END IF
       !   END IF
       !  END DO
       ! END DO
       !END DO

! CHECK TO SEE IF WE ARE DOING HAMILTONIAN MIXING
       CALL CHECK_INPUTS
       IF(MIXING1.EQ.2)THEN
         DO ISPN=1,NSPN
           !YY. Calculate hamiltonian. This will create HAMOLD. 
           CALL OVERLAP(2)
         END DO
         OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED')
         READ(99)NREC,MXX
         CLOSE(99) 
         IF(sichammix) THEN
          CALL SICHAMMIXDRV(NREC)
         ELSE
          CALL HAMMIXDRV(NREC)
         END IF
       ENDIF
!       DO ISPN=1,NSPN
!         CALL OVERLAP(2)
!       END DO
!       IF(NITER.NE.1)THEN  
!        OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED')
!        REWIND(99)
!        READ(99)NREC,MXX
!        DO JSPN=1,NSPN            
!         READ(99)(HSTOR(IREC,JSPN),IREC=1,NREC)
!C        DO IREC=1,NREC
!          IREC=0
!          DO IREP=1,N_REP
!           DO I=1,NS_TOT(IREP)
!            DO J=I,NS_TOT(IREP)
!             IREC=IREC+1
!             HSTOR(IREC,JSPN)=(HSTOR(IREC,JSPN)-SICSCI(IREC,JSPN))*    
!     &        DIAG(I,IREP)*DIAG(J,IREP)
!            END DO
!           END DO
!          END DO
!C        END DO
!        END DO
!C CALL MIXING..
!C00000000000000000000000000000000000000000000000000000000000000000000000
!        NSPTS=0
!        DO JSPN=1,NSPN
!         IREC=0
!         DO IREP=1,N_REP
!          DO I=1,NS_TOT(IREP)
!C          TIMES=1.0D0
!           DO J=I,NS_TOT(IREP)
!            NSPTS=NSPTS+1
!            IREC=IREC+1
!            POT(NSPTS)=HSTOR(IREC,JSPN)!*TIMES
!C           TIMES=1.0D0 
!           END DO
!          END DO
!         END DO
!        END DO
!        PRINT*,'IREC,NSPTS:',IREC,NSPTS
!        AVERAGE=.TRUE.
!        AVG= 0.15
!        IF(AVERAGE)THEN
!         call system('echo -1 T > AVRGDAT')
!         OPEN(12,FILE='GOTOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
!         REWIND(12)
!         ITER=0
!         READ(12,END=500)MREC,MSPN,ITER
!         IF(ITER.EQ.0)GOTO 500
!         IF((MREC.NE.NREC).OR.(MSPN.NE.NSPN))THEN
!          PRINT *,'NEWWAVE2020: FILE GOTOLD IS UNUSABLE'
!          CALL STOPIT
!         END IF
!         READ(12)(POTIN(IPTS), IPTS=1,NREC*NSPN)
!         IF(AVG.GT.0.0D0)THEN
!          PRINT '(A)','BROYDEN MIXING OF POTENTIAL',NSPTS
!          NSPTS=NREC*NSPN
!          CALL MIXING(ITER,AVG,NSPTS,NAMES)
!         ELSE
!          PRINT '(A)','SIMPLE LINEAR MIXING OF POTENTIAL',NSPTS
!          DO IPTS=1,NSPTS
!           POTIN(IPTS)=(1.0D0-AVG)*POTIN(IPTS)+AVG*POT(IPTS)
!          END DO  
!         END IF
!         DO 495 IPTS=1,NSPTS
!          POT(IPTS)=POTIN(IPTS)
!  495    CONTINUE
!  500    ITER=ITER+1
!         REWIND(12)  
!         WRITE(12)NREC,NSPN,ITER
!         WRITE(12)(POT(IPTS), IPTS=1,NREC*NSPN)
!         CLOSE(12)
!        ELSE
!         ITER=0
!         OPEN(12,FILE='GOTOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
!         REWIND(12)
!         WRITE(12)NREC,NSPN,ITER
!         WRITE(12)(POT(IPTS), IPTS=1,NREC*NSPN)
!         CLOSE(12)
!        END IF
!        NSPTS=0
!        DO JSPN=1,NSPN
!         IREC=0
!C        DO IREC=1,NREC
!          DO IREP=1,N_REP
!           DO I=1,NS_TOT(IREP)
!c           TIMES=1.0D0
!            DO J=I,NS_TOT(IREP)
!             NSPTS=NSPTS+1
!             IREC=IREC+1
!             HSTOR(IREC,JSPN)=POT(NSPTS)/DIAG(I,IREP)/DIAG(J,IREP)
!c    & /TIMES
!c            TIMES=1.0D0
!            END DO
!           END DO
!          END DO
!C        END DO
!        END DO
!        PRINT*,'A,IREC,NSPTS:',IREC,NSPTS,NREC
!C00000000000000000000000000000000000000000000000000000000000000000000000
!C CCCCCCCCCCCC
!        REWIND(99)
!        WRITE(99)NREC,MXX
!        DO JSPN=1,NSPN            
!         WRITE(99)(HSTOR(IREC,JSPN),IREC=1,NREC)
!        END DO
!        CLOSE(99)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!       END IF
!       DEALLOCATE(DIAG)   

       DO 240 ISPN=1,NSPN
        IF(MIXING1.NE.2) THEN
          CALL OVERLAP(2)
        ELSE
          PRINT*,'ISPN:',ISPN
          OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED')
          REWIND(99)
          READ(99)NREC,MXXX
          DO JSPN=1,ISPN            
            READ(99)(HSTOR(IREC,2),IREC=1,NREC)
          END DO
          CLOSE(99)
        ENDIF
        PRINT*,'MADE IT HERE'
        CALL OVERLAP(1)
C HSTOR(IREC,1)=OVERLAP
C HSTOR(IREC,2)=HAMILTONIAN FOR SPIN 1
        KBAS=0
        DO 130 IREP=1,N_REP
         N_VIRT(IREP,ISPN)=0
         NBAS=NS_TOT(IREP)
         IF (NBAS.GT.NDH) THEN
          PRINT *,'NEWWAVE: NDH MUST BE AT LEAST: ',NBAS
          CALL STOPIT
         END IF
         ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error allocating Ham'
         ENDIF
         ALLOCATE(AOVER(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error allocating Overlap'
         ENDIF
         ALLOCATE(AEVAL(NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error allocating Eval'
         ENDIF
C
         AOVER=0.0D0 
         LBAS=NBAS
c        IF(EXIST)THEN
c        READ(57)LBAS
c        READ(57,END=56)((AOVER(J,I),J=1,LBAS),I=1,LBAS)
c        PRINT*,'NEWWAVE SCISP'
c        DO I=1,LBAS
c         PRINT '(10F10.4)',(AOVER(I,J),J=1,LBAS)
c        END DO
c        END IF
         IF(NBAS.NE.LBAS)THEN
          PRINT*, 'LBAS.NE.NBAS',NBAS,LBAS,IREP,N_REP
          CALL STOPIT
         END IF
 56      CONTINUE
         DO 80 IBAS=1,NBAS
          DO 70 JBAS=IBAS,NBAS
           KBAS=KBAS+1
           if(sichammix .and. (MIXING1.EQ.2)) then
            AHAM (JBAS,IBAS)=HSTOR(KBAS,2)
            AHAM (IBAS,JBAS)=HSTOR(KBAS,2)
           else
            AHAM (JBAS,IBAS)=HSTOR(KBAS,2)-SICSCI(KBAS,ISPN)
            AHAM (IBAS,JBAS)=HSTOR(KBAS,2)-SICSCI(KBAS,ISPN)
           end if
           AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
           AOVER(IBAS,JBAS)=HSTOR(KBAS,1)
   70     CONTINUE
   80    CONTINUE
C ASSUMING GETSIC HAS CONSTRUCTED SCISSOR:
C
C
         OPEN(55,FILE='4TEST',FORM='UNFORMATTED')
         WRITE(55)NBAS
         WRITE(55)((AOVER(J,I),J=1,NBAS),I=1,NBAS)
         !CALL DIAGGE(NDH,NBAS,AHAM,AOVER,AEVAL,SC1,1) 
         CALL DIAGGE(NBAS,NBAS,AHAM,AOVER,AEVAL,SC1,1)
         WRITE(55)((AHAM(J,I),J=1,NBAS),I=1,43)
         CLOSE(55)
         IF(.false.)THEN
          PRINT*,'NEWWAVE EVALS:'
          DO I=1,NBAS
           PRINT*,AEVAL(I)
          END DO
         ENDIF
C WRITE HAM TO TEMP_PSI_COEF
         DO 120 ISAV=1,NBAS                   
          IF(.false.) PRINT*,ISPN,ISAV,AEVAL(ISAV)
          EVLO(ISAV,IREP,ISPN)=AEVAL(ISAV)
          DO 110 IB=1,NBAS
           PSI_COEF(IB,ISAV,IREP,ISPN)=AHAM(IB,ISAV)
C          PSI_COEF(IB,ISAV,IREP)=AHAM(IB,ISAV)
  110     CONTINUE
  120    CONTINUE
         IF(allocated(AHAM))THEN
          DEALLOCATE(AHAM,STAT=IERR)
          IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error deallocating Ham'
          ENDIF
         END IF
         DEALLOCATE(AOVER,STAT=IERR)
         IF(IERR.NE.0)THEN
          WRITE(6,*)'newwave_serial:Error deallocating Overlap'
         ENDIF
         DEALLOCATE(AEVAL,STAT=IERR)
         IF(IERR.NE.0)THEN
          WRITE(6,*)'newwave_serial:Error deallocating Eval'
         ENDIF
  130   CONTINUE
  240  CONTINUE    
       IF(EXIST)CLOSE(57)

       efermi(1)=-1.0d30
       efermi(2)=-1.0d30
       if(mdrn.eq.3)then
        do is=1,nspn
         do ir=1,n_rep
          do ib=1,n_occ(ir,is)
           if(occv(ib,ir,is).gt.0.00001)then
            efermi(is)=max(efermi(is),evlo(ib,ir,is))
           end if
          end do 
         end do 
        end do 
       else
        if(mdrn.eq.1)ntimes=1
        if(mdrn.eq.2)ntimes=2
        if(nspn.eq.1)ntimes=1
        spfact=dfloat(2/nspn)
        do 250 itimes=1,ntimes
         i1=itimes
         i2=itimes
         if(ntimes.eq.1.and.nspn.eq.2)i2=nspn
         etarget=e_up+e_dn
         if(ntimes.eq.2.and.itimes.eq.1)etarget=e_up        
         if(ntimes.eq.2.and.itimes.eq.2)etarget=e_dn       
         if(nspn.eq.1)then
          etarget=e_up+e_dn   
         end if
         emin= 1.0d30
         emax=-1.0d30
         do is=i1,i2
          do ir=1,n_rep
           do ib=1,ns_tot(ir)
            emin=min(evlo(ib,ir,is),emin)
            emax=max(evlo(ib,ir,is),emax)
           end do
          end do
         end do
         emin=emin-0.5
         emax=emax+0.5
         ef=(emin+emax)/2.0d0
         do itry=1,50 
          electrons=0.0d0
          do is=i1,i2
           do ir=1,n_rep
            do ib=1,ns_tot(ir)
             occv(ib,ir,is)=1.0d0/(1.0d0+exp((evlo(ib,ir,is)-ef)/temp))
             electrons=electrons
     &                +occv(ib,ir,is)*(ndmrep(ir)/ldmrep(ir))*spfact
             efermi(is)=ef
            end do
           end do
          end do
          if(electrons.le.etarget)then
           emin=ef
           ef=(emax+ef)/2.0d0
          else
           emax=ef
           ef=(emin+ef)/2.0d0
          end if
         end do
  250   continue
        if(efermi(nspn).lt.-1.0d15)efermi(nspn)=efermi(1)
 260    format(A)
 270    format(4I5,35x,' # of occupied states for spin,rep,dim')
 290    format(4I5,45x,' (Spin, Rep, Nbas, Deg)')
 280    format(8g12.4)
        n_occ=0
        do is=1,nspn
         do ir=1,n_rep
          n_occ(ir,is)=min(ns_tot(ir),2)
          do ib=1,ns_tot(ir)
           if(occv(ib,ir,is).gt.0.000001d0)then
            n_occ(ir,is)=min(ib+2,ns_tot(ir))
           end if 
          end do
         end do 
        end do
       end if
       rewind(97)
       write(97,'(A)')FLINE
       do is=1,nspn
        do ir=1,n_rep
         write(97,270)n_occ(ir,is),is,ir,ndmrep(ir)
         write(97,280)(occv(ib,ir,is),ib=1,n_occ(ir,is))
        end do 
       end do
       do is=1,nspn
        do ir=1,n_rep
         write(97,"(43A1,' NEW REP ',43A1)")('*',i=1,86)
         write(97,290)is,ir,ns_tot(ir),ndmrep(ir) 
         write(97,*)' '
         write(97,280)(evlo(ib,ir,is),ib=1,ns_tot(ir))
        end do
       end do
       write(97,'(95A1)')('*',i=1,95)
       write(97,"(30x,F12.6,' Fermi Levels',F12.6)")
     &                   (efermi(JSPN),JSPN=1,NSPN)
       write(97,'(95A1)')('*',i=1,95)
       call wf_control(97,0)
       ud(1)='u'
       ud(2)='d'
       ch='<=='
       jwf=0
       do is=1,nspn
        do ir=1,n_rep
         do ib=1,ns_tot(ir)
          jwf=jwf+1
          knt(3,jwf)=ib
          knt(2,jwf)=ir
          knt(1,jwf)=is
         end do
        end do
       end do
       do ind=1     ,jwf 
        do jnd=ind+1,jwf  
         if(evlo(knt(3,ind),knt(2,ind),knt(1,ind)).gt.
     &      evlo(knt(3,jnd),knt(2,jnd),knt(1,jnd)))then
          do k=1,3
           jwi=knt(k,ind)
           knt(k,ind)=knt(k,jnd)
           knt(k,jnd)=jwi
          end do
         end if
        end do
        if(occv(knt(3,ind),knt(2,ind),knt(1,ind)).gt.0.000001)lst=ind
       end do
       do ind=1,lst+10*NSPN
        ib=knt(3,ind)
        ir=knt(2,ind)
        is=knt(1,ind)
        if(ind.eq.lst.or.(ind.lt.lst.and.occv(ib,ir,is).le.0.999)
     & .or.ind.eq.lst+1)then
         if(ind.eq.lst)then
          GAP=evlo(ib,ir,is)
          comment='HOMO'
         else if(ind.eq.lst+1)then
          GAP=GAP-evlo(ib,ir,is)
          write(comment(1:23),'(A11,F12.6)')'LUMO, GAP =',GAP
C               comment='LUMO, GAP = ',evlo(ib,ir,is)-evlo(ib,ir,is-1)
         else
          comment='HOLE'
         end if
         write(97,92)ind,ir,ndmrep(ir),ud(is),evlo(ib,ir,is)
     & ,occv(ib,ir,is),ch(1:3),comment,control_wave,ib,ir,is
         write(47,47)ib,ir,is,ns_tot(ir),evlo(ib,ir,is),occv(ib,ir,is)
         write(47,47)(psi_coef(jb,ib,ir,is),jb=1,ns_tot(ir))
         write(47,*)' '
 47      format(12g15.6)
         if(ind.eq.lst)write(97,*)' '
        else
         comment=' '
         write(97,92)ind,ir,ndmrep(ir),ud(is),evlo(ib,ir,is)
     &,occv(ib,ir,is),ch(1:3),comment,control_wave,ib,ir,is
         write(47,47)ib,ir,is,ns_tot(ir),evlo(ib,ir,is),occv(ib,ir,is)
         write(47,47)(psi_coef(jb,ib,ir,is),jb=1,ns_tot(ir))
         write(47,*)' '
        end if
       end do
 92    format(i5,' Rep:',I3,' Deg:',I3,' Spin: ',A1,
     &' Energy',F15.6,' Occ:',F15.6,' ',A3,A25,A10,I6,I3,I3)
  300  FORMAT(15G15.3)
       close(97)
       close(47)
       write(fline(1:18),'(A15,I3.3)')'cp EVALUES EVAL',NITER
       call system(fline(1:18))
       nwfs(1)=0
       nwfs(2)=0
       nwf    =0
       itt    =0
       e_up=0.0d0
       e_dn=0.0d0
       do is=1,nspn
        do ir=1,n_rep
         do ib=1,n_occ(ir,is)
          itt=itt+1
          indrep(itt)=ir  
            ndeg(itt)=ndmrep(ir)
           ntemp(itt)=is
         evalsav(itt)=evlo(ib,ir,is)
          occtmp(itt)=occv(ib,ir,is)
          if(is.eq.1)then
           e_up=e_up+occv(ib,ir,is)*ndmrep(ir)/ldmrep(ir)
          else
           e_dn=e_dn+occv(ib,ir,is)*ndmrep(ir)/ldmrep(ir)
          end if
          occupancy(itt)=occtmp(itt)/ldmrep(ir)
          if(itt.gt.max_virt_per_sym*max_rep*mxspn)then
           print*,'error, itt, newwave_2020'
           call stopit
          end if
          trace=trace+(2/nspn)*evlo(ib,ir,is)*occv(ib,ir,is)
     &         *ndmrep(ir)/ldmrep(ir)
          do ii=1,ndmrep(ir)
           nwf=nwf+1
           if(nwf.gt.max_occ)then
            print*,'max_occ must be >',nwf
            call stopit
           end if
           nwfs(is)=nwfs(is)+1
           evlocc(nwf)=evlo(ib,ir,is)
          end do
         end do
        end do
       end do
       print*,'Number of Spin Up Electrons:',e_up,' ELECTRONS'
       print*,'Number of Spin Dn Electrons:',e_dn,' ELECTRONS'
       if(nspn.eq.1)e_dn=e_up
       print*,nwf,nwfs,lst+10*nspn
       if(nwf.gt.max_occ)then
        print*,'max_occ must be >',nwf
        call stopit
       end if
       CALL GTTIME(TIMEB)
       PRINT*,'TIME IN NEW NEWWAVE:',TIMEB-TIMEA

C deallocate local arrays that were created for diagonalization

       DEALLOCATE(NDEG,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating NDEG'
       ENDIF
       DEALLOCATE(INDREP,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating INDREP'
       ENDIF
       DEALLOCATE(NTEMP,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating NTEMP'
       ENDIF
       DEALLOCATE(OCCTMP,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating OCCTMP'
       ENDIF
       RETURN
       END
C      FUNCTION PSI_COEF(ARG1,ARG2,ARG2,ARG4)

