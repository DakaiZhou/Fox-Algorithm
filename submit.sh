#!/bin/bash
#
#SBATCH --nodes=7
#SBATCH --ntasks=49
#SBATCH --exclusive
#SBATCH --partition=NODE2008


MPIDIR=/cluster/mpi/openmpi/1.6.5-gcc4.8.2/
GCCDIR=/cluster/gcc/4.8.2/
export PATH=$PATH:${MPIDIR}/bin:${GCCDIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPIDIR}/lib:${GCCDIR}/lib:${GCCDIR}/lib64

cat A.txt B.txt > AB.txt
mpirun ./ex7 < AB.txt