input="../matrices_path/MM.csv"
{
	read
	i=1
	while IFS=',' read -r num str		
	do
	  ./mkl_spmv $str | tee -a ../DATA/mkl_all.csv
	  i='exp $i + 1'
	 done
} < "$input"
