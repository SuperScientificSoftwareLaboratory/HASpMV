#include "omp.h"
#include "common.h"
#include "biio.h"
int main(int argc, char **argv){
    if (argc < 1)
    {
        printf("Run the code by './bind_spmv matrix.cbd 90'.\n");
        return 0;
    }
	char  *filename;
    filename = argv[1];
    // int test = atoi(argv[3]);
    int m;
    int n;
    int nnz;
    int ls;
    MAT_PTR_TYPE *csrRowPtr;
    MAT_PTR_TYPE *csrColIdx;
    MAT_VAL_TYPE *csrVal;
    MAT_VAL_TYPE *x,*y_ref,*y;
    struct timeval t1, t2, t3, t4;
    read_Dmatrix_32(&m, &n, &nnz, &csrRowPtr, &csrColIdx, &csrVal, &ls, filename);
    // int *row_map = (int*)malloc(m * sizeof(m));//old -> new
    int *cost_sum = (int*)malloc((m + 1) * sizeof(int));
    memset(cost_sum,0,sizeof cost_sum);
    // printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnz);
    x = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * n);
    for(int i = 0; i < n; i++)
    x[i] = rand() % 10 + 1;
    y_ref = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * m);
    memset(y_ref, 0, sizeof(MAT_VAL_TYPE) * m);


    //reference
    double time_min,time = 0;
    for(int re = 0;re < REPEAT;re++){
        gettimeofday(&t3, NULL);
        #pragma omp parallel for
        for (int i = 0; i < m; i++)
        {
            y_ref[i] = 0;
            for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
            {
                y_ref[i] += csrVal[j] * x[csrColIdx[j]];
            }
        }
        gettimeofday(&t4, NULL);
        double time_tp = (t4.tv_sec - t3.tv_sec) * 1000.0 + (t4.tv_usec - t3.tv_usec) / 1000.0;
        if(re == 0) time_min = time_tp;
        else time_min = time_min < time_tp ? time_min : time_tp;
        time += time_tp;
    }
    
    time /= REPEAT;
    double GFlops_avg = (nnz * 2) / (time * 1000000);
    double GFlops = (nnz * 2) / (time_min * 1000000);


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
    
    printf("naive,%s,%.3lf,%.3lf\n",name,GFlops,GFlops_avg);
    free(name);


    free(y_ref);
    free(x);

    free(csrRowPtr);
    free(csrColIdx);
    free(csrVal);

    return 1;
}