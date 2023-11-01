cd HASpMV/
python3 genmake.py 12
make
sh haspmv-all-12.sh
cd ../CSR5_avx2
make
sh csr5-all.sh
cd ../merge-spmv
make cpu_spmv
sh Merge-all.sh
cd ../mkl
make
sh mkl-all.sh
# 下面要写画图步骤