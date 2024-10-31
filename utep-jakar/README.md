# A simple instruction on setting up and running the FLOSIC code on UTEP-Jakar cluster

**The files in the utep-jakar directory is a part of the University of Texas at El Paso (UTEP) special topics class course materials.**

## Setting up required libraries and compiling the source code

First, on UTEP Jakar you need to load the modules with a command as follows,

> module load gnu12 mpich openblas 

Optional: You can append the module load command in your ~/.bashrc file, such that these modules are loaded automatically each time you log in to Jakar.

Installing lapack to your home directory can be done with the included setup_jakar.sh script. Run it with the following command. You need to do this only once.

> bash setup_jakar.sh

Optional: You can set up gnuplot with gnuplot_install script. Simply run the as follows. You need to do this only once.

> bash gnuplot_install.sh  
 
For the FLOSIC (PR2020) code compilation, copy of makefile for Jakar is included (Makefile.jakar). Copy it to the PublicRelease_2020/flosic/ directory before the code compilation.

> cp Makefile.jakar ../flosic/makefile
> cd ../flosic/
> make

The make command should be done in the source code directory (i.e. PublicResease_2020/flosic/ directory).

Once the compilation finishes, you will find a binary file named *nrlmol_exe*. Copy this binary file to a directory where you will run your calculation.

## Parallel calculation using MPI

For a large molecule especially SIC calculations, you may want to run your calculation using a parallelized code. For this, simply set *MPI=Y* in your Makefile if not already done so. Then you can recompile the code with

> make clean -f *Makefile.jakar* && make -f *Makefile.jakar*
 

## Running a job on Jakar using the slurm job scheduler

To submit/queue a job to computing nodes, first, prepare a directory with input files, binary, and job script. An included job.sl is a minimum job script that uses 10 CPU cores on general partition. For using multicores, your code must be compiled with MPI. It can be customized for your calculation. Once you have your job script ready, run

> sbatch job.sl

You can check the status of your jobs with the squeue command.

> squeue -u *username*

In the status, you see R (running), CF (configuring), CG (completed), etc. 
If you need to cancel a job, do

> scancel *job-id-number*
