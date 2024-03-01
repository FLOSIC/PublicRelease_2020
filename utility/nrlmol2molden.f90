! nrlmol2molden.f90  (version 0.1)
!
! A utility code for NRLMOL/FLOSIC codes
! Convert WFOUT/XMOL.DAT/ISYMGEN/EVALUES to a molden file
!
! compilation: gfortran nrlmol2molden.f90 -o nrlmol2molden.bin
! usage      : copy the input files in a directory and run the binary
! inputs     : WFOUT, XMOL.DAT, ISYMGEN, EVALUES
! output     : converted.molden  
!
! Memo: Tested for non symmetry case.
!
! Yoh Y. 2024-2-29

program main

implicit none

! psuedo PARAM
integer,parameter :: MAX_VIRT_PER_SYM = 2100
integer,parameter :: NDH              = 1000
integer,parameter :: max_bare         = 25
integer,parameter :: max_con          = 15
integer,parameter :: extrabas_on      = 0 ! 1: include extra basis 0: not include.

logical :: ldo

integer :: ISPN, NSPN, IWF, NWF, NWFS(1:2), NBASF(2,3), I_REP, N_REP
integer :: I, IL, ICOUNT, I_OCC, ITOT, J, JAT, JATOMS, JJ, JSET
integer :: iatoms , noatoms, zatom, L, NALP, NCOUNT, NFIND, NSETS, NSUM

real*8  :: pos(3), ALPSET(max_bare), CONSET(max_bare,max_con,3)
real*8  :: FJL, charge_e, charge_n, EFERMI, tempe

integer,allocatable :: N_OCC(:,:), NS_TOT(:)
real*8, allocatable :: EVLOCC(:), OCCUPANCY(:,:), PSI_COEF(:,:,:,:)

character*80 :: strbuffer
character*8  :: strread(3)
character*2  :: lett
character*2  :: z2lett !function
character*1  :: NSPD(3) 
DATA NSPD/'s','p','d'/


! ********************************************************************

! Read WFOUT
print *,'Processing WFOUT....'
OPEN(98,FILE='WFOUT',FORM='UNFORMATTED',STATUS='UNKNOWN')
REWIND(98)
READ(98) NSPN
READ(98) NWF,(NWFS(ISPN),ISPN=1,NSPN)
allocate(EVLOCC(NWF))
READ(98) (EVLOCC(IWF), IWF=1,NWF)
READ(98) N_REP
allocate(N_OCC(N_REP,2), NS_TOT(N_REP))
allocate(OCCUPANCY(MAX_VIRT_PER_SYM*N_REP,2))
allocate(PSI_COEF(NDH,MAX_VIRT_PER_SYM,N_REP,2))

DO ISPN=1,NSPN
  ITOT=0
  DO I_REP=1,N_REP
    READ(98) N_OCC(I_REP,ISPN),NS_TOT(I_REP)
    READ(98) (OCCUPANCY(I_OCC+ITOT,ISPN), I_OCC=1,N_OCC(I_REP,ISPN))
    ITOT=ITOT+N_OCC(I_REP,ISPN)
    DO IWF=1,N_OCC(I_REP,ISPN)
      READ(98) (PSI_COEF(I,IWF,I_REP,ISPN), I=1,NS_TOT(I_REP))
    END DO
  END DO
END DO

CLOSE(98)
! End reading WFOUT


!Read FERMI level
open(19,file='EVALUES')
rewind(19)
ldo=.true.
do while(ldo)
  read(19,'(A)') strbuffer
  !print *, strbuffer(1:20)
  if(strbuffer(2:6) == 'FERMI') then
    read(strbuffer,*) strread(1), strread(2), EFERMI, strread(3), tempe
    print *,'EF found in EVALUES', EFERMI
    ldo=.false.
  end if
end do
close(19)


! Outpuf file
open(98,FILE='converted.molden',FORM='FORMATTED',STATUS='UNKNOWN')
rewind(98)


! Nuclear coordinate block

write(98,*) '[Molden Format]'
write(98,*) '[Atoms] Angs'

print *, 'Processing XMOL.DAT....'
open(18,file='XMOL.DAT')
rewind(18)
read(18,*) noatoms
read(18,*) strbuffer
do iatoms=1, noatoms
  read(18,*) zatom, pos(1:3)
  write(98,'(2X, A5, 2I7, 3F17.6 )') z2lett(zatom), iatoms, zatom, pos(1:3)
end do
close(18)


! Basis set block

WRITE(98,*)'[6D]'
WRITE(98,*)'[10F]'
WRITE(98,*)'[GTO]'

!Read ISYMGEN
print *, 'Processing ISYMGEN....'
OPEN(20,file='ISYMGEN')
REWIND(20)
ICOUNT = 0
NCOUNT = 0
READ(20,*)NSETS
DO JSET = 1, NSETS
  READ(20,*)charge_e,charge_n
  READ(20,*)
  READ(20,*)NFIND  !no. of atoms of type xx
  DO JATOMS = 1,NFIND
    READ(20,*)     !loop over ALL-OXY00? etc.
  END DO
  READ(20,*)
  READ(20,*)NALP
  READ(20,*)(NBASF(1,L),L=1,3)
  READ(20,*)(NBASF(2,L),L=1,3)
  READ(20,*)(ALPSET(J),J=1,NALP)
  READ(20,*)
  DO L = 1, 3
    DO IL = 1, NBASF(1,L)+NBASF(2,L)
      READ(20,*)(CONSET(J,IL,L),J=1,NALP)
      READ(20,*)
    END DO
  END DO

  DO JAT = 1, NFIND !NUPP
    ICOUNT = ICOUNT + 1
    WRITE(98,*)NCOUNT+JAT,'0'
    DO L = 1, 3
      DO IL = 1, NBASF(1,L)+NBASF(2,L)*extrabas_on
        NSUM = 0
        DO JJ = 1, NALP
          IF(CONSET(JJ,IL,L).NE.0) THEN
            NSUM = NSUM + 1
          END IF
        END DO
        WRITE(98,'(1X,A1,3X,I2,2X,A4)')NSPD(L),NSUM,'1.00'

        DO J = 1, NALP
          IF(CONSET(J,IL,L).NE.0) THEN
            !--- NORMALIZE BASIS SETS ---
            FJL = 1.0/(3.14**(0.75))
            FJL = FJL*ALPSET(J)**(0.75+0.5*(L-1))
            FJL = FJL*(SQRT(2.0)**(L-1))
            WRITE (98,'(2D18.10)')ALPSET(J),CONSET(J,IL,L)/FJL
          END IF
        END DO
      END DO
    END DO
    WRITE(98,*)' '
  END DO
  NCOUNT = ICOUNT
END DO

READ(20,*)
READ(20,*)
CLOSE(20)


! Molecular orbital block

WRITE(98,*) '[MO]'

print *,'Writing MOs....'
DO ISPN=1,NSPN
  ITOT=0
  DO I_REP=1,N_REP
    !WRITE(98,*) N_OCC(I_REP,ISPN),NS_TOT(I_REP)
    !WRITE(98,*) (OCCUPANCY(I_OCC+ITOT,ISPN), I_OCC=1,N_OCC(I_REP,ISPN))
    !ITOT=ITOT+N_OCC(I_REP,ISPN)
    DO IWF=1,N_OCC(I_REP,ISPN)
      !Symmetry label
      write(strread(1),'(I5)') IWF
      write(strread(2),'(I5)') I_REP
      write(strread(3),'(I5)') ISPN
      strbuffer =  'Sym= '//trim(adjustl(strread(1)))//' - '//trim(adjustl(strread(2)))//' - '//trim(adjustl(strread(3)))
      !WRITE(98,*) 'Sym= ', IWF,'-', I_REP, '-', ISPN 
      WRITE(98,'(A)') strbuffer
      WRITE(98,'(A4,2F20.8)') 'Ene=', EVLOCC(IWF), EFERMI !EFERMI(ISPN)
      IF(ISPN==1) THEN
        WRITE(98,'(A11)') 'Spin= Alpha'
      ELSE
        WRITE(98,'(A10)') 'Spin= Beta'
      END IF
      WRITE(98,'(A6, F20.8)') 'Occup=', OCCUPANCY(IWF+ITOT,ISPN)
      DO I=1,NS_TOT(I_REP)
       WRITE(98,'(I10, F20.8)') I, PSI_COEF(I,IWF,I_REP,ISPN)
      END DO
    END DO

    ITOT=ITOT+N_OCC(I_REP,ISPN)
  END DO
END DO

close(98)

deallocate(PSI_COEF)
deallocate(OCCUPANCY)
deallocate(N_OCC, NS_TOT)
deallocate(EVLOCC)

print *,'Completed!'

end program main
!***********************************************************


! A function to convert atomic number Z to the atom letter
function z2lett(input)

integer     :: input
character*2 :: z2lett
character*2 :: tableletter(0:112)
DATA tableletter  /"X", "H", "He", "Li", "Be", "B", "C",         &
     "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",&
     "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",   &
     "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", &
     "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",  &
     "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",  &
     "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", &
     "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",  &
     "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", &
     "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",  &
     "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", &
     "Hs", "Mt", "Ds", "Rg", "Cn" /


z2lett = tableletter(input)


end function z2lett
