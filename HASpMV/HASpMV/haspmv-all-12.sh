export OMP_NUM_THREADS=24
export GOMP_CPU_AFFINITY="0-23"
input="../matrices_path/MM.csv"
{
	read
	while IFS=',' read -r num str		
	do
		./haspmv $str 750 | tee -a ../DATA/haspmv-all-12.csv
	done
} < "$input"
