! UTEP Electronic Structure Lab (2020)
!> YY. Made this to get kinetic density
!
!> This code is based on densold subroutine:
!> DENSOLD BASED ON OLD VERSION OF APOTNL BY M. PEDERSON AND D. POREZAG
!> ATTENTION: FIRST TWO ARRAYS OF COMMON BLOCK TMP1 MUST BE IDENTICAL IN 
!> DENSOLD AND APOTNL SINCE THEY ARE USED TO PASS DENSITY AND COULOMB POT

! I also need to get LPTS to get a certan tau(MPTS) value
subroutine gettau_par(LPTS,tau,LDFTF,MXXD,KXXS,KXXO)
!use xtmp1,only : rhog
!use xtmp1a,only : phig
use debug1
use mesh1,only : rmsh,nmsh
use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT, ISPN, NSPN
use common5,only : PSI, NWF, NWFS
!USE XTMP2A,ONLY : TAU=>TAUTOT
!include 'PARAMS'
include 'PARAMA2'
integer :: ICON, IERR, IFNCT, IGR, ILOC, IPTS, ISHDUM, ISHELLA, ISIZE, &
     &     IWF, I_POS, JBEG, JLOC, JPTS, JWF, J_POS, KPTS, LI, LMAX1,  &
     &     LPV, L_NUC, MPTS, MU, M_NUC, NDERV, NGRAD, NMAX, NPV
real*8 ::  FACTOR, TIMEGORB, TTIME1, TTIME2
save
parameter (NMAX=MPBLOCKXC)


! RHOG(IPTS,1, 1)= rho_up   
! RHOG(IPTS,2, 1)= d rho_up/dx
! RHOG(IPTS,3, 1)= d rho_up/dy
! RHOG(IPTS,4, 1)= d rho_up/dz
! RHOG(IPTS,1, 2)= rho_dn   
! RHOG(IPTS,2, 2)= d rho_dn/dx
! RHOG(IPTS,3, 2)= d rho_dn/dy
! RHOG(IPTS,4, 2)= d rho_dn/dz

real*8, allocatable :: PSIG(:,:,:),PTS(:,:),GRAD(:,:,:,:,:),RVECA(:,:)
logical, allocatable :: ICOUNT(:,:)

logical LGGA,IUPDAT
dimension ISIZE(3)

integer,intent(in)  :: LPTS
real*8, intent(out) :: tau(MPBLOCKXC,NSPN)   !kinetic density
!real*8 :: RHOG(NMSH,4,NSPN)   !rhog(:,1,:) is for debugging. I should get the same density
logical, intent(in) :: LDFTF  !T:DFT, F:SIC
integer, intent(in) :: MXXD,KXXS,KXXO
data ISIZE/1,3,6/

!TIMEGORB=0.0D0
allocate(PSIG(NMAX,10,MAX_OCC),STAT=IERR);   if(IERR/=0)write(6,*)'DENSOLD:ERROR ALLOCATING PSIG'
allocate(PTS(NSPEED,3),STAT=IERR);           if(IERR/=0)write(6,*)'DENSOLD:ERROR ALLOCATING NSPEED'
allocate(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR);   if(IERR/=0)write(6,*)'DENSOLD:ERROR ALLOCATING GRAD'
allocate(RVECA(3,MX_GRP),STAT=IERR);         if(IERR/=0)write(6,*)'DENSOLD:ERROR ALLOCATING RVECA'
allocate(ICOUNT(MAX_CON,3),STAT=IERR);       if(IERR/=0)write(6,*)'DENSOLD:ERROR ALLOCATING ICOUNT'


!This is for meta-GGA. Hence no need to check IGGA
LGGA= .TRUE.
NGRAD=4 !To get kinetic density, I only need IGR=2-4

! LOOP OVER ALL POINTS
TTIME1=0.0D0
TTIME2=0.0D0
!LPTS=0
!10 continue
if(LPTS+NMAX.LT.NMSH)THEN
 MPTS=NMAX
else
 MPTS=NMSH-LPTS
end if

! INITIALIZE PSIG AND RHOB

PSIG(1:MPTS,1:NGRAD,1:NWF)=0.0d0
tau(1:MPTS,1:NSPN)=0.0d0
ISHELLA=0

! FOR ALL CENTER TYPES

do IFNCT=1,NFNCT
 LMAX1=LSYMMAX(IFNCT)+1

! FOR ALL POSITIONS OF THIS CENTER

 do I_POS=1,N_POS(IFNCT)
  ISHELLA=ISHELLA+1

! GET SYMMETRY INFO

  call OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
  if(NWF.GT.MAX_OCC)THEN
   write(6,*)'GETTAU: MAX_OCC MUST BE AT LEAST:',NWF
   call STOPIT
  end if

! FOR ALL EQUIVALENT POSITIONS OF THIS ATOM

  do J_POS=1,M_NUC

! UNSYMMETRIZE 

!YY sic needs unravel2, dft needs unravel here. 3/19/18
   IF(LDFTF)THEN
    call UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),RVECA,L_NUC,1)
   ELSE
! These three input parameters (MXXD, KXXS, and KXXO) are same as coupot_sic.
!   MXXD=0 - 0: SIC -1:DFT
!   KXXS LSPN, spin index
!   KXXO IORB (FOD index). Shifted by +NWFX(1) for spin down. 
    CALL UNRAVEL2(MXXD,KXXS,KXXO,IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),RVECA,L_NUC,1)
   !CALL UNRAVEL2(MXXD,KXXS,KXXO,IFNCT,ISHELLA,K_SITEI,RIDT(1,ISHELLA),RVECI,LST,1)
   END IF
   if(L_NUC.NE.M_NUC)THEN
    write(6,*)'GETTAU: PROBLEM IN UNRAVEL'
    call STOPIT
   end if

! FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK

   KPTS=0
   do JPTS=1,MPTS,NSPEED
    NPV=MIN(NSPEED,MPTS-JPTS+1)
    do LPV=1,NPV
     KPTS=KPTS+1
     PTS(LPV,1)=RMSH(1,LPTS+KPTS)-RVECA(1,J_POS)
     PTS(LPV,2)=RMSH(2,LPTS+KPTS)-RVECA(2,J_POS)
     PTS(LPV,3)=RMSH(3,LPTS+KPTS)-RVECA(3,J_POS)
    end do

! GET ORBITS AND DERIVATIVES

    NDERV=0
    if (LGGA) NDERV=2
!   call GTTIME(TIME3)
    call GORBDRV(NDERV,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
!   call GTTIME(TIME4)
!   TIMEGORB=TIMEGORB+TIME4-TIME3

! UPDATING ARRAY PSIG

    if (IUPDAT) THEN
     IPTS=JPTS-1
     ILOC=0
     do LI=1,LMAX1
      do MU=1,ISIZE(LI)
       do ICON=1,N_CON(LI,IFNCT)
        ILOC=ILOC+1
        if (ICOUNT(ICON,LI)) THEN    ! <=== This will ignore gauss == 0
         do IWF=1,NWF
          FACTOR=PSI(ILOC,IWF,1)
          if(abs(FACTOR) .GT. 1.0d-10) then
           do IGR=1,NGRAD
            do LPV=1,NPV
             PSIG(IPTS+LPV,IGR,IWF)=PSIG(IPTS+LPV,IGR,IWF) +FACTOR*GRAD(LPV,IGR,MU,ICON,LI)
             ! GRAD is stored as LI and MU combinations
             ! LI=1 MU=1, LI=2 MU=1-3, LI=3 MU=1-6
            end do
           end do
          end if  
         end do  
        end if
       end do  
      end do  
     end do
    end if
   end do
  end do
 end do
end do


! UPDATING RHOG, START WITH DENSITY 

!do ISPN=1,NSPN
! JBEG= (ISPN-1)*NWFS(1) 
! do JWF=1,NWFS(ISPN)
!  JLOC=JWF+JBEG
!  do IPTS=1,MPTS
!   RHOG(LPTS+IPTS,1,ISPN)=RHOG(LPTS+IPTS,1,ISPN) +PSIG(IPTS,1,JLOC)**2
!  end do
! end do
!end do

! UPDATE DERIVATIVES IF GGA CALCULATION
         
if (LGGA) THEN
 do 96 ISPN=1,NSPN
  JBEG= (ISPN-1)*NWFS(1)
  do 94 JWF=1,NWFS(ISPN)
   JLOC=JWF+JBEG

! GRADIENT 

!  do IGR=2,4
!   do IPTS=1,MPTS
!    RHOG(LPTS+IPTS,IGR,ISPN)=RHOG(LPTS+IPTS,IGR,ISPN)+2.0d0*PSIG(IPTS,1,JLOC)*PSIG(IPTS,IGR,JLOC)
!    !PHIG(LPTS+IPTS,ISPN)    =PHIG(LPTS+IPTS,ISPN)    +PSIG(IPTS,IGR,JLOC)*PSIG(IPTS,IGR,JLOC)
!   end do
!  end do

! Kinetic Density
   do IPTS=1,MPTS
    tau(IPTS,ISPN)=tau(IPTS,ISPN)+0.5d0*( PSIG(IPTS,2,JLOC)**2 + PSIG(IPTS,3,JLOC)**2 + PSIG(IPTS,4,JLOC)**2 )
   end do

  94 continue
 96 continue
end if


!LPTS=LPTS+MPTS
!if (LPTS .LT. NMSH) goto 10
!continue

deallocate(PSIG,STAT=IERR);   if(IERR/=0)write(6,*)'DENSOLD:ERROR DEALLOCATING PSIG'
deallocate(PTS,STAT=IERR);    if(IERR/=0)write(6,*)'DENSOLD:ERROR DEALLOCATING NSPEED'
deallocate(GRAD,STAT=IERR);   if(IERR/=0)write(6,*)'DENSOLD:ERROR DEALLOCATING GRAD'
deallocate(RVECA,STAT=IERR);  if(IERR/=0)write(6,*)'DENSOLD:ERROR DEALLOCATING RVECA'
deallocate(ICOUNT,STAT=IERR); if(IERR/=0)write(6,*)'DENSOLD:ERROR DEALLOCATING ICOUNT'

!CALL TIMOUT('GETTAU TOTAL OBINFO:                ',TTIME1)
!CALL TIMOUT('GETTAU TOTAL UNRAVEL:               ',TTIME2)
!CALL TIMOUT('GETTAU GORBDRV:                     ',TIMEGORB)
RETURN
end subroutine
