The directory where this files reside must be assigned to a system variable called NRLMOL_BASIS so that the program has access to the .basis files

You may set it in your .bashrc file by adding this line:

export NRLMOL_BASIS=/whatever/basis/directory

In C shell, add this line to .cshrc (or .cshrc.ext)

setenv NRLMOL_BASIS /whatever/basis/directory

In both cases the /whatever/basis/directory is the path to where the .basis file are located

There is an extra directory here called 'original' which contains the original .dat files that were downloaded from the website https://bse.pnl.gov/bse/portal

Some later basis were downloaded from http://www.emsl.pnl.gov/forms/basisform-orig.html

In order to process them, the extnbas2.c file must be compiled with:

gcc extnbas2.c -o process_basis -lm

When run, that binary processes a single .dat file that is given as a command line argument, creating the .basis file from the .dat file.

example: ./process_basis some_basis.dat

will create a file calle some_basis.basis

In order to process all the files, run the python script process_bases.py which will run the process_basis for all .dat files in the directory where it resides.

The .basis files are written in ISYMGEN format wich is obtained from the alphas and coefficients from the original .dat file, the entry in the .basis file is calculated as entry=coef*alpha^(0.75+(l-1)/2) where:
l=1 for S functions
l=2 for P functions
l=3 for D functions

