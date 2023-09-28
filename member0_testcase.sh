#!/bin/sh

### ********************************************************************
### ********************************************************************
### **                                                                **
### **  Example shell script to run CiTG-TU Delft/Deltares SWAN       **
### **  executable in parallel with OpenMP via SGE on linux-cluster.  **
### **  c 2016 CiTG-TU Delft/Deltares                                 **
### **  author: Menno Genseberger                                     **
### **                                                                **
### ********************************************************************
### ********************************************************************

### The "-cwd" requests execution of this script from the current
### working directory; without this, the job would be started from the
### user's home directory.
#$ -cwd

### The "-q" requests one of the possible queues of SGE.
### At Deltares, for quad-core H6 e3 nodes
### (bare metal, Intel quad-core Xeon E3-1276 v3 nodes, 3.6 GHz core)
### this is the "normal-e5" queue. To select this queue remove the first
### two # in the next line, the third # should remain.
#$ -q normal-e3-c7
### At Deltares, for quad-core H6 e5 nodes
### (virtual, Intel quad-core Xeon E5-2667 v3 nodes, 3.2 GHz core)
### this is the "normal-e5" queue. To select this queue remove the first
### two # in the next line, the third # should remain.
###$ -q normal-e5-c7
###$ -q test-c7

### 
### To use this script replace your_job_name by a clear name (for
### instance  my_job), this name will appear in the queue of the linux-
### cluster.
### The name of this SGE job is explicitly set to another name;
### otherwise the name of the SGE script itself would be used. The name
### of the job also determines how the jobs output files will be called. 
#$ -N swn_reference_test

### To run the executable swan_4091AB_8_del_l64_i11_omp.exe
### (4091AB_8 = 40.91AB_8 version number, del = Deltares version,
###  l64 = linux 64 bits platform, i11 = Intel 11 Fortran compiler, and
###  omp = OpenMP) references to the appropriate libraries for Intel 11
### compiler and precompiled NetCDF libraries have to be set.
### Note: these are specific setting for the Deltares H6 linux-cluster.

module load swan/41.31A.1_intel18.0.3_timers
swan_omp_exe=swan_omp.exe

export OMP_NUM_THREADS=4

### Some general information available via SGE.
echo ----------------------------------------------------------------------
echo Run of
echo $swan_omp_exe
echo with OpenMP on linux-cluster.
echo SGE_O_WORKDIR : $SGE_O_WORKDIR
echo HOSTNAME : $HOSTNAME
echo OMP_NUM_THREADS : $OMP_NUM_THREADS
echo ----------------------------------------------------------------------
cp member0_testcase.swn INPUT

$swan_omp_exe

cp PRINT member0_testcase.prt
rm -f PRINT
rm -f INPUT
rm -f norm_end