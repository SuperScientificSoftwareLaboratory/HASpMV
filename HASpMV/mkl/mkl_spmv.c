#include "common.h"
#include "biio.h"
#include <assert.h>
#include <math.h>
#include "mkl.h"
#include "mkl_spblas.h"
#include <omp.h>
#include <stdlib.h>

int cmp(const void* a,const void* b){
  double ret = *(double*)a - *(double*)b;
  if(ret > 0){
    return 1;
  }
  else if(ret < 0){
    return -1;
  }
  return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Run the code by './mkl_spmv matrix.mtx'.\n");
        return 0;
    }
    struct timeval t1, t2, t3, t4;
    // printf("--------------------------------!!!!!!!!------------------------------------\n");
    char  *filename;
    filename = argv[1];
    //printf("MAT: -------------- %s --------------\n", filename);

    // int nthreads=8;
    // omp_set_num_threads(nthreads);
    
    
    int m;
    int n;
    int nnz;
    int isSymmetric;
    MAT_PTR_TYPE *csrRowPtr;
    int *csrColIdx;
    MAT_VAL_TYPE *csrVal;
    read_Dmatrix_32(&m, &n, &nnz, &csrRowPtr, &csrColIdx, &csrVal,&isSymmetric, filename);
  //  printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnz);
    char *name = (char*)malloc(30*sizeof(char));
    for(int i = strlen(filename);i > 0;i--){
        if(filename[i] == '/'){
            // printf("%s\n",filename + i);
            strncpy(name,filename + i + 1,strlen(filename) - i - 5);
            name[strlen(filename) - i - 5] = '\0';
            // name = filename + i + 1;
            break;
        }
    }
    MAT_VAL_TYPE *x = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * n);
    for (int i = 0; i < n; i++)
    x[i] = rand() % 10 + 1;

    MAT_VAL_TYPE *y = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * m);
    memset(y, 0, sizeof(MAT_VAL_TYPE) * m);
    MAT_VAL_TYPE *y_ref = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * m);
    memset(y_ref, 0, sizeof(MAT_VAL_TYPE) * m);
    for (int i = 0; i < m; i++)
    {
    	y_ref[i] = 0;
        for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
        {
            y_ref[i] += csrVal[j] * x[csrColIdx[j]];
        }
    }
    
    struct matrix_descr descrA;
    sparse_matrix_t csrA;
    sparse_status_t status;
    int exit_status = 0;
    
    gettimeofday(&t1, NULL);
    for(int i = 0;i < 1;i++){
        status = mkl_sparse_d_create_csr ( &csrA,
                                        SPARSE_INDEX_BASE_ZERO,
                                        m,  // number of rows
                                        n,  // number of cols
                                        csrRowPtr,
                                        csrRowPtr+1,
                                        csrColIdx,
                                        csrVal );
        if (status != SPARSE_STATUS_SUCCESS) {
            printf(" Error in mkl_sparse_d_create_csr: %d \n", status);
            exit_status = 1;
            goto exit;
        }

        descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
        status = mkl_sparse_set_mv_hint(csrA, SPARSE_OPERATION_NON_TRANSPOSE, descrA, BENCH_REPEAT );
        if (status != SPARSE_STATUS_SUCCESS && status != SPARSE_STATUS_NOT_SUPPORTED) {
            printf(" Error in set hints for mv: mkl_sparse_set_mv_hint: %d \n", status);
        }

        status = mkl_sparse_optimize ( csrA );
        if (status != SPARSE_STATUS_SUCCESS) {
            printf(" Error in mkl_sparse_optimize: %d \n", status);
            exit_status = 1;
            goto exit;
        }
    }
    gettimeofday(&t2, NULL);
    double pretime_mkl  = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
    // pretime_mkl /= 40;

    //printf("------------------------------oneAPI-MKL--------------------------------\n");
    //kernel
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    double time_sum = 0;
    double time_min = 0x7ffffff;
    double time;
    for (int k = 0; k < REPEAT; k++){
        gettimeofday(&t1, NULL);
    	status = mkl_sparse_d_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1, csrA,descrA, x, 0, y );
        gettimeofday(&t2, NULL);
        time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        time_sum += time;
        time_min = time_min > time ? time : time_min;
    }
    time_sum /= REPEAT;
    double GFlops = (nnz * 2) / (time_min * 1000000);
    double GFlops_avg = (nnz * 2) / (time_sum * 1000000);
    //printf("mkl_max,%d,%.3f\n",nnz,GFlops);
    printf("mkl,%s,%d,%.3lf,%.3lf,%.3lf,%.3lf\n",name,nnz,GFlops,GFlops_avg,pretime_mkl,time_sum);
    // printf("%.3lf,%.3lf\n",pretime_mkl,time_sum);
   
    
    if (status != SPARSE_STATUS_SUCCESS) {
        printf(" Error in Task 1 mkl_sparse_d_mv: %d \n", status);
        exit_status = 1;
        goto exit;
    }

    // validate x
    double accuracy = 1e-4;
    double ref = 0.0;
    double res = 0.0;

    for (int i = 0; i < m; i++)
    {
        ref += abs(y_ref[i]);
        res += abs(y[i] - y_ref[i]);
        // printf("index = %d  %8.1lf %8.1lf\n", i, x[i], x_ref[i]);
    }

    res = ref == 0 ? res : res / ref;
    
    /*if (res < accuracy)
    {
        flag = 1;
        printf("MKL_y check passed! |vec-vecref|/|vecref| = %8.2e\n", res);
        
    }
    else
        printf("MKL_y check _NOT_ passed! |vec-vecref|/|vecref| = %8.2e\n", res);*/

    // FILE *fouttime = fopen("mkl_spmv_res.csv", "a");
    //fprintf(fouttime, "%s,%i,%i,%i,%f,%f,%f,%i\n", 
          //          filename, m, n, nnz, pretime_mkl, time_mkl, MKL_GFlops, flag);
    // fclose(fouttime);
    
    //printf("time = %.3lf\n", time_mkl);
    //printf("GFlops: %.3lf\n", MKL_GFlops);
    // FILE *fouttime = fopen("mklspmv_res.csv", "a");
    // fprintf(fouttime, "%s,%i,%i,%i,%f,%f\n", 
    //                   filename, m, n, nnz, time_mkl, MKL_GFlops);
    // fclose(fouttime);

    free(x);
    free(y);
    free(y_ref);

    free(csrRowPtr);
    free(csrColIdx);
    free(csrVal);

    exit:
    // Release matrix handle and deallocate matrix
    mkl_sparse_destroy ( csrA );

    return exit_status;
}
