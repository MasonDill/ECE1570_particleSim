#!/bin/bash -l


#SBATCH --time=00:20:00					# adjust the amount of time as necessary
#SBATCH --job-name=serial
#SBATCH --account=ece1570_2022f
#SBATCH --output auto-omp.%j.stdout
#SBATCH --error auto-omp.%j.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cluster=smp
#SBATCH --partition=smp


module purge
module load gcc/5.4.0

rm openmp.txt
srun -n 1 -c 1 ./serial -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=1
srun -n 1 -c 1 ./qtree_omp -p 1 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=2
srun -n 1 -c 2 ./qtree_omp -p 2 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=4
srun -n 1 -c 4 ./qtree_omp -p 4 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=6
srun -n 1 -c 6 ./qtree_omp -p 6 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=12
srun -n 1 -c 12 ./qtree_omp -p 12 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=18
srun -n 1 -c 18 ./qtree_omp -p 18 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=24
srun -n 1 -c 24 ./qtree_omp -p 24 -n 500 -no -s openmp.txt
export OMP_NUM_THREADS=2
srun -n 1 -c 2 ./qtree_omp -p 2 -n 1000 -no -s openmp.txt
export OMP_NUM_THREADS=4
srun -n 1 -c 4 ./qtree_omp -p 4 -n 2000 -no -s openmp.txt
export OMP_NUM_THREADS=6
srun -n 1 -c 6 ./qtree_omp -p 6 -n 3000 -no -s openmp.txt
export OMP_NUM_THREADS=12
srun -n 1 -c 12 ./qtree_omp -p 12 -n 6000 -no -s openmp.txt
export OMP_NUM_THREADS=18
srun -n 1 -c 18 ./qtree_omp -p 18 -n 9000 -no -s openmp.txt
export OMP_NUM_THREADS=24
srun -n 1 -c 24 ./qtree_omp -p 24 -n 12000 -no -s openmp.txt
./autograder -v openmp -s openmp.txt
