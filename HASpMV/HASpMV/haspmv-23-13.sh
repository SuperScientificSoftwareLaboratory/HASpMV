export OMP_NUM_THREADS=32
export GOMP_CPU_AFFINITY="0-31"
input="../matrices_path/examples.csv"
{
	read
	while IFS=',' read -r num str		
	do
		./haspmv $str 850 | tee -a ../DATA/haspmv-23-13.csv
	done
} < "$input"
