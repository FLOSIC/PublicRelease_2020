! UTEP Electronic Structure Lab (2020)
module SICFLAG

 LOGICAL :: LSICF,MESH_SAV

end module SICFLAG


module LOCORB
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

 real*8  :: TMAT(MAX_OCC,MAX_OCC,MXSPN)
! real*8, allocatable :: TMAT(:,:,:)
 integer :: MORB(2)
 real*8  :: ZSIC
 integer :: IRBSIC
end module LOCORB


module MOCORB
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

 real*8  :: SLAT(MAX_OCC,MAX_OCC,MXSPN)
! real*8, allocatable :: SLAT(:,:,:)
 integer :: NFRM(2)
 real*8  :: ZTZL
 integer :: JJJJJJ
end module MOCORB


module SICMAT
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

!  real*8 :: SIC(MAX_OCC,MAX_OCC,MXSPN)
 real*8 :: DERSIC(3,MAX_OCC,MX_CNT)
! real*8, allocatable :: SIC(:,:,:)
#ifdef MPI_3
 real*8,pointer :: SIC_COL(:)
#else
 real*8, allocatable :: SIC_COL(:)
#endif
 real*8, allocatable :: SIC(:,:,:)
 REAL*8 :: FMAT(MAX_OCC,MAX_OCC,4,2)
 real*8, allocatable :: ZPOT(:,:,:)
 real*8, allocatable :: ZMGGA(:,:,:) !4*NMSH*ORBITALS
 real*8, allocatable :: ZMGGAS(:,:,:)

end module SICMAT


module FRM
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

 real*8  :: BFRM(3,MAX_OCC,MXSPN)
 real*8  :: RESULTS(13,MAX_OCC,MXSPN)
 integer :: LFRM(MXSPN)
 real*8  :: DEBDAX(3,MAX_OCC,MXSPN)
end module FRM


module HMATSIC
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

 real*8 :: OVTM(MAX_OCC,MAX_OCC,2)
 real*8 :: HMTM(MAX_OCC,MAX_OCC,2)
end module HMATSIC


module FOCENT
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

 real*8 :: CFRM(3,MAX_OCC,MXSPN)
end module FOCENT


module DIRECTIONS
 integer :: NSPN_SKIP
end module DIRECTIONS


module ORBENG
! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  

 real*8 :: EVALOCC(MAX_OCC)
end module ORBENG

module DIAGV1

integer            :: NORB
real*8,allocatable :: PHIRES(:,:),PHIRES_TMP(:,:)

end module DIAGV1

module scaledpzsic
!Flag to turn on scaled SIC methods (LSIC,OSIC,GSIC,sdSIC). 
!Also used to revert to PZSIC.
logical ::  scaledsic = .FALSE. !T: scaled sic, F: pzsic

!LSIC of J. Chem. Phys. 151, 214108 (2019)
logical ::  LSICON   = .FALSE. 

!OrbSIC of J. Chem. Phys. 124, 094108 (2006) adapted to FLOSIC as J. Chem. Phys. 152, 174112 (2020)
logical ::  ORBSICON = .FALSE. 

!GSIC of J. P. Perdew
logical ::  GSICON   = .FALSE. 

!Slater averaging of SIC potential
logical ::  AVGSICON = .FALSE. 

!sdSIC of J. Chem. Phys. 152, 214109 (2020)
logical ::  SDSICON  = .FALSE. 

!Perturbative or quasi-SCF (scaled potential) used with the scaled SIC methods.
!quasi-SCF is supported for LSIC, OSIC, and sdSIC. GSIC is perturbative only.
logical ::  oneshot  = .FALSE. !T: oneshot, F: quasi-SCF 

!Arrays used for scaled SIC methods.
real*8, allocatable :: AVGSIC(:) !For average sic potential
REAL*8,ALLOCATABLE :: SICEXCS(:), SICEXC(:) !SIC_XC energy density
REAL*8 :: SDSICSAV(50,2)
end module scaledpzsic
