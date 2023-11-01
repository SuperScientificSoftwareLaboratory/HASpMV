input="../matrices_path/examples.csv"
{
	read
	i=1
	while IFS=',' read -r num str		
	do
	  ./spmv $str | tee -a ../DATA/CSR5_23.csv
	  i='exp $i + 1'
	 done
} < "$input"
