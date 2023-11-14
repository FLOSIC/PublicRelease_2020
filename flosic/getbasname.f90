! UTEP Electronic Structure Lab (2020)
SUBROUTINE GETBASNAME(BASIS,FILENAME)
  use global_inputs, only : basis_filename
  IMPLICIT NONE
  INTEGER :: BASIS
  CHARACTER*40 :: FILENAME
  SELECT CASE(BASIS)
  CASE(2)
    FILENAME='/3-21G.basis'
  CASE(3)
    FILENAME='/3-21Gs.basis'
  CASE(4)
    FILENAME='/3-21GSP.basis'
  CASE(5)
    FILENAME='/3-21Gs_Polarization.basis'
  CASE(6)
    FILENAME='/3-21ppG.basis'
  CASE(7)
    FILENAME='/3-21ppGs.basis'
  CASE(8)
    FILENAME='/6-311G2df2pd.basis'
  CASE(9)
    FILENAME='/6-311G.basis'
! Removed
!  CASE(10)
  CASE(11)
    FILENAME='/6-311Gss.basis'
  CASE(12)
    FILENAME='/6-311Gss_Polarization.basis'
  CASE(13)
    FILENAME='/6-311pGs.basis'
  CASE(14)
    FILENAME='/6-311ppG2d2p.basis'
  CASE(15)
    FILENAME='/6-311ppG3df3pd.basis'
  CASE(16)
    FILENAME='/6-311ppGss.basis'
  CASE(17)
    FILENAME='/6-31G_3df3pd.basis'
  CASE(18)
    FILENAME='/6-31G-Blaudeau.basis'
  CASE(19)
    FILENAME='/6-31Gs.basis'
  CASE(20)
    FILENAME='/6-31Gs-Blaudeau.basis'
  CASE(21)
    FILENAME='/6-31Gs_Polarization.basis'
  CASE(22)
    FILENAME='/6-31Gss.basis'
  CASE(23)
    FILENAME='/6-31Gss_Polarization.basis'
  CASE(24)
    FILENAME='/6-31pGs.basis'
  CASE(25)
    FILENAME='/6-31ppG.basis'
  CASE(26)
    FILENAME='/6-31ppGs.basis'
  CASE(27)
    FILENAME='/6-31ppGss.basis'
  CASE(28)
    FILENAME='/6-32G.basis'
  CASE(29)
    FILENAME='/Ahlrichs_Coulomb_Fitting.basis'
  CASE(30)
    FILENAME='/Ahlrichs_Polarization.basis'
  CASE(31)
    FILENAME='/Ahlrichs_pVDZ.basis'
  CASE(32)
    FILENAME='/Ahlrichs_TZV.basis'
  CASE(33)
    FILENAME='/Ahlrichs_VDZ.basis'
  CASE(34)
    FILENAME='/Ahlrichs_VTZ.basis'
  CASE(35)
    FILENAME='/DZVP2.basis'
  CASE(36)
    FILENAME='/DZVP.basis'
  CASE(37)
    FILENAME='/GAMESS_PVTZ.basis'
  CASE(38)
    FILENAME='/GAMESS_VTZ.basis'
  CASE(39)
    FILENAME='/Huzinaga_polarization.basis'
  CASE(40)
    FILENAME='/IGLO-II.basis'
  CASE(41)
    FILENAME='/IGLO-III.basis'
  CASE(42)
    FILENAME='/McLeanChandler_VTZ.basis'
  CASE(43)
    FILENAME='/MIDI_Huzinaga.basis'
  CASE(44)
    FILENAME='/MINI_Huzinaga.basis'
  CASE(45)
    FILENAME='/Partridge_Uncontracted_1.basis'
  CASE(46)
    FILENAME='/Partridge_Uncontracted_2.basis'
  CASE(47)
    FILENAME='/Partridge_Uncontracted_3.basis'
  CASE(48)
    FILENAME='/Partridge_Uncontracted_4.basis'  
  CASE(49)
    FILENAME='/Saddlej.basis'
  CASE(50)
    FILENAME='/STO-2G.basis'
  CASE(51)
    FILENAME='/STO-3G.basis'
  CASE(52)
    FILENAME='/STO-3Gs.basis'
  CASE(53)
    FILENAME='/STO-3Gs_Polarization.basis'
  CASE(54)
    FILENAME='/STO-6G.basis'
  CASE(55)
    FILENAME='/TZVP.basis'
  CASE(56)
    FILENAME='/stuttgart_rlc_ecp.basis'
  CASE(57)
    FILENAME='/lanl2dz_ecp.basis'
  CASE(58)
    FILENAME='/lanl2dzdp_ecp.basis'
  CASE(59)
    FILENAME='/UGBS.basis'
  CASE(100)
    WRITE(FILENAME,'(3A)') '/', trim(basis_filename), '.basis'
  END SELECT

END SUBROUTINE GETBASNAME
