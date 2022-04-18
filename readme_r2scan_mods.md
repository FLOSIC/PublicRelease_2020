# FLO-SIC with r<sup>2</sup>SCAN

## Maintainers: James W. Furness , Aaron D. Kaplan
## Validation: Chandra Shahi, Raj Sah, Pradeep Bhetwal

This repo is a fork of the FLO-SIC 2020 public release (https://github.com/FLOSIC/PublicRelease_2020)
containing files implementing the novel and efficient r<sup>2</sup>SCAN meta-GGA:
>J. W. Furness, A. D. Kaplan, J. Ning, J. P. Perdew, and J. Sun,
> "Accurate and Numerically Efficient r<sup>2</sup>SCAN Meta-Generalized Gradient Approximation", \
>J. Phys. Chem. Lett. **11**, 8208-8215 (2020). [DOI: 10.1021/acs.jpclett.0c02405] \
>https://doi.org/10.1021/acs.jpclett.0c02405.

r<sup>2</sup>SCAN subroutines written by James Furness (@JFurness1, james.w.furness.1@gmail.com). Licensed CC0 (without any access restrictions).
For further information, please see the metadata in *r2scan.f* or the public r<sup>2</sup>SCAN code repository \
https://gitlab.com/dhamil/r2scan-subroutines.

Two extant files have been modified:\
*setdftyp.ftn* \
*subvlxc_libxc.F91* \
The new keywords are *R2SCAN* for accessing r2 SCAN, and *SDBG*, which performs SCAN calculations using the subroutines in *r2scan.f*.
(Please note that these subroutines permit self consistent calculations with SCAN, rSCAN, r++SCAN, r<sup>2</sup>SCAN and r<sup>4</sup>SCAN via modification of integer flags therein.)

*SDBG* can be used to validate the implementation of the subroutines, as has been done in the **R2SCAN_tests** directory.
*scan_comps_AE6.xslx* compares the existing implementation of SCAN in FLO-SIC with *SDBG*, and shows agreement is essentially exact.

It appears that Mesh size 1 (intermediate size, intended for rSCAN) can be used with r<sup>2</sup>SCAN.
See the tests in *r2scan-test-N-Ne-apr12.xlsx* for the atoms N, Ne, and the AE6 set.

**R2SCAN_tests/NaCl_dissociation** contains a mesh test of the dissociation of NaCl.
This again confirms that Mesh 1 can be used for small atoms with r<sup>2</sup>SCAN.

A few auxiliary files are included in **R2SCAN_tests**: *test_atoms.py* and *cbs_extrap.py*.
These are PySCF scripts used to generate reference exchange-only and exchange-correlation total energies with PySCF.
These reference energies can then be compared against FLO-SIC, as is done in \
**R2SCAN_tests/atom_tests/atom_calcs_cart_sp_norm_w_cbs.xlsx**.
Agreement is generally to 1e-6 Hartree or better.
