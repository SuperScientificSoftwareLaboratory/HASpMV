input="../matrices_path/examples.csv"
{
	read
	i=1
	while IFS=',' read -r num str		
	do
	  ./_cpu_spmv_driver $str | tee -a ../DATA/Merge_23.csv
	  i='exp $i + 1'
	 done
} < "$input"
