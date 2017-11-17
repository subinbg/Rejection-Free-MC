#!/bin/sh

#PBS -N MonteCarlo
#PBS -l nodes=xb003

cd ${PBS_O_WORKDIR}
rm coord.dat
#mpirun -n 8 ./LJ
mpirun ./LJ
sleep 20
exit
