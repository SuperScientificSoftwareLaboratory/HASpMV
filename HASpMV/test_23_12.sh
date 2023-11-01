cd HASpMV/
python3 genmake.py 12
make clean
make
# sh haspmv-23-12.sh
sh tech.sh
cd ../CSR5_avx2
make
sh csr5-23.sh
cd ../merge-spmv
make cpu_spmv
sh Merge-23.sh
cd ../mkl
make
sh mkl-23.sh
# 下面要写画图步骤