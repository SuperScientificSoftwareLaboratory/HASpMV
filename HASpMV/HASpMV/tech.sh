export OMP_NUM_THREADS=24
export GOMP_CPU_AFFINITY="0-23"
input="../matrices_path/examples.csv"
{
	read
	while IFS=',' read -r num str		
	do
		./naive $str | tee -a ../DATA/naive.csv
		./partition $str 750 | tee -a ../DATA/partition.csv
		./reorder $str 750 | tee -a ../DATA/reorder.csv
		./avx2 $str 750 | tee -a ../DATA/avx2.csv
		./loop $str 750 | tee -a ../DATA/loop.csv
	  i='exp $i + 1'
	done
} < "$input"



