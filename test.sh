for i in {1..12}
do
    let a="2*i"
    export OMP_NUM_THREADS=$a
    ./serial -n 10000 -s summary_omp.txt
done