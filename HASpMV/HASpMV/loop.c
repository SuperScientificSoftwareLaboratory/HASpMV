#include "omp.h"
#include "common.h"
#include "biio.h"
int binary_search_cost(int *cost_sum,int left,int right,int cost_bound){//find row
    int mid = left + right >> 1;
    while(left < right){
        if(cost_sum[mid] > cost_bound) right = mid;
        else left = mid + 1;
        mid = left + right >> 1;
    }
    return mid;
}
int search_nnz(MAT_PTR_TYPE *csrColIdx,int col_begin,int col_end,int bound){//find nnz
    int ben = -1;
    int cost = 0;
    for(int j = col_begin;j < col_end;j++){
        // cost += 1;
        if(cost >= bound){
            return j;
        }
        if(csrColIdx[j] / 8 > ben){
            cost ++;
            ben = csrColIdx[j] / 8;
            if(cost == bound){
                int p = j + 1;
                while(csrColIdx[p] / 8 == ben){
                    p++;
                }
                return p - 1;
            }
        }
    }
    return col_end;
}
double avx2_kernel_p(int begin,int end,MAT_PTR_TYPE *csrColIdx,MAT_VAL_TYPE *csrVal,MAT_VAL_TYPE *x){
    int length = end - begin;
    if(length < 4){
        double sum = 0;
        for(int j = begin;j < end;j++)
            sum += csrVal[j] * x[csrColIdx[j]];
        return sum;
    }
    int remainder = 0; 
    double res_y;
    __m256d res,val,vec;
    res = _mm256_setzero_pd();
    // printf("%d %d %d\n",remainder,begin,end - remainder);
    int j = begin;
    if(length < 125){
        remainder = length % 4;
        int loop = length / 4;
        for(int i = 0;i < loop;i++){
            val = _mm256_loadu_pd(&csrVal[j]);
            vec = _mm256_set_pd(x[csrColIdx[j + 3]],x[csrColIdx[j + 2]],x[csrColIdx[j + 1]],x[csrColIdx[j]]);
            res = _mm256_fmadd_pd(val,vec,res);
            j += 4;
        }
    }
    else{
        int loop = length / 8;
        remainder = length % 8;
        for(int i = 0;i < loop;i++){
            val = _mm256_loadu_pd(&csrVal[j]);
            vec = _mm256_set_pd(x[csrColIdx[j + 3]],x[csrColIdx[j + 2]],x[csrColIdx[j + 1]],x[csrColIdx[j]]);
            res = _mm256_fmadd_pd(val,vec,res);
            j += 4;
            val = _mm256_loadu_pd(&csrVal[j]);
            vec = _mm256_set_pd(x[csrColIdx[j + 3]],x[csrColIdx[j + 2]],x[csrColIdx[j + 1]],x[csrColIdx[j]]);
            res = _mm256_fmadd_pd(val,vec,res);
            j += 4;
        }
    }
    
    res = _mm256_hadd_pd(res,res);
    res_y = res[0] + res[2];
    for(int jj = end - remainder;jj < end;jj++){
        res_y += csrVal[jj] * x[csrColIdx[jj]];
    }
    return res_y;
}
double avx2_kernel_e(int begin,int end,MAT_PTR_TYPE *csrColIdx,MAT_VAL_TYPE *csrVal,MAT_VAL_TYPE *x){
    int length = end - begin;
    if(length < 4){
        double sum = 0;
        for(int j = begin;j < end;j++)
            sum += csrVal[j] * x[csrColIdx[j]];
        return sum;
    }
    int remainder = 0; 
    double res_y;
    __m256d res,val,vec;
    res = _mm256_setzero_pd();
    // printf("%d %d %d\n",remainder,begin,end - remainder);
    int j = begin;
    if(length < 175){
        remainder = length % 4;
        int loop = length / 4;
        for(int i = 0;i < loop;i++){
            val = _mm256_loadu_pd(&csrVal[j]);
            vec = _mm256_set_pd(x[csrColIdx[j + 3]],x[csrColIdx[j + 2]],x[csrColIdx[j + 1]],x[csrColIdx[j]]);
            res = _mm256_fmadd_pd(val,vec,res);
            j += 4;
        }
    }
    else{
        int loop = length / 8;
        remainder = length % 8;
        for(int i = 0;i < loop;i++){
            val = _mm256_loadu_pd(&csrVal[j]);
            vec = _mm256_set_pd(x[csrColIdx[j + 3]],x[csrColIdx[j + 2]],x[csrColIdx[j + 1]],x[csrColIdx[j]]);
            res = _mm256_fmadd_pd(val,vec,res);
            j += 4;
            val = _mm256_loadu_pd(&csrVal[j]);
            vec = _mm256_set_pd(x[csrColIdx[j + 3]],x[csrColIdx[j + 2]],x[csrColIdx[j + 1]],x[csrColIdx[j]]);
            res = _mm256_fmadd_pd(val,vec,res);
            j += 4;
        }
    }
    
    res = _mm256_hadd_pd(res,res);
    res_y = res[0] + res[2];
    for(int jj = end - remainder;jj < end;jj++){
        res_y += csrVal[jj] * x[csrColIdx[jj]];
    }
    return res_y;
}
int main(int argc, char **argv){
    if (argc < 2)
    {
        printf("Run the code by './bind_spmv matrix.cbd 90'.\n");
        return 0;
    }
	char  *filename;
    filename = argv[1];
    int part = atoi(argv[2]);
    // int test = atoi(argv[3]);
    int pl[CORE_P * 2 + 1],pr[CORE_P * 2 + 1];
    int nnz_l[CORE_E + CORE_P];
    int nnz_r[CORE_E + CORE_P];
    int m;
    int n;
    int nnz;
    int ls;
    MAT_PTR_TYPE *csrRowPtr;
    MAT_PTR_TYPE *csrColIdx;
    MAT_VAL_TYPE *csrVal;
    MAT_VAL_TYPE *x,*y_ref,*y;
    MAT_PTR_TYPE *old_csrRowPtr;
    MAT_PTR_TYPE *old_csrColIdx;
    MAT_VAL_TYPE *old_csrVal;
    struct timeval t1, t2, t3, t4;
    read_Dmatrix_32(&m, &n, &nnz, &old_csrRowPtr, &old_csrColIdx, &old_csrVal,&ls, filename);
    // printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnz);
    csrRowPtr = (MAT_PTR_TYPE*)malloc((m + 1) * sizeof(MAT_PTR_TYPE));
    csrColIdx = (MAT_PTR_TYPE*)malloc(nnz * sizeof(MAT_PTR_TYPE));
    csrVal = (MAT_VAL_TYPE*)malloc(nnz * sizeof(MAT_VAL_TYPE));
    // int *row_map = (int*)malloc(m * sizeof(m));//old -> new
    int *cost_sum = (int*)malloc((m + 1) * sizeof(int));
    memset(cost_sum,0,sizeof cost_sum);
    // printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnz);
    x = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * n);
    for(int i = 0; i < n; i++)
    x[i] = rand() % 10 + 1;
    y = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * m);
    memset(y, 0, sizeof(MAT_VAL_TYPE) * m);
    y_ref = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * m);
    memset(y_ref, 0, sizeof(MAT_VAL_TYPE) * m);
    // y_tp = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * m);
    // memset(y_tp, 0, sizeof(MAT_VAL_TYPE) * m);
    int* row_begin_nnz = (int*)malloc(m * sizeof(int));
    
    gettimeofday(&t1,NULL);
    //re
    int front_row = 0;
    int tail_row = m - 1;
    int front_nnz = 0;
    int tail_nnz = nnz;
    csrRowPtr[0] = 0;
    csrRowPtr[m] = nnz;
    for(int i = 0;i < m;i++){
        int nnz_num = old_csrRowPtr[i + 1] - old_csrRowPtr[i];
        if(nnz_num < 1000){
            row_begin_nnz[front_row] = old_csrRowPtr[i];
            csrRowPtr[front_row+1] = front_nnz + nnz_num;
            front_nnz += nnz_num;
            front_row++;
        }
        else{
            row_begin_nnz[tail_row] = old_csrRowPtr[i];
            tail_nnz -= nnz_num;
            csrRowPtr[tail_row] = tail_nnz;
            tail_row--;
        }
    }
    #pragma omp parallel for
    for(int i = 0;i < m;i++){
        int bas1 = csrRowPtr[i];
        int bas2 = row_begin_nnz[i];
        int len = csrRowPtr[i + 1] - csrRowPtr[i];
        for(int j = 0;j < len;j++){
            csrColIdx[bas1 + j] = old_csrColIdx[bas2 + j];
            csrVal[bas1 + j] = old_csrVal[bas2 + j];
        }
    }
    int cost_x = 0;
    int ben = -1;

    for(int i = 0;i < m;i++){
        cost_x = 0;
        ben = -1;
        for(int j = csrRowPtr[i];j < csrRowPtr[i + 1];j++){
            if(csrColIdx[j] / 8 > ben){
                cost_x++;
                ben = csrColIdx[j] / 8;
            }
        }
        cost_sum[i] = cost_x;
    }
    for(int i = 1;i < m;i++){
        cost_sum[i] += cost_sum[i - 1];
    }

    int COST = cost_sum[m - 1];
    // printf("line:%d\n",COST);
    int costp = COST * (part / 1000.0);
    int coste = COST - costp;
    int bound_mid = binary_search_cost(cost_sum,0,m - 1,costp);//row
    // printf("%d %d\n",costp,cost_sum[bound_mid]);
    // printf("%d %d %d\n",costp,bound_mid,cost_sum[bound_mid]);
    int gapp = ceil(costp / CORE_P);
    int gape = ceil(coste / CORE_E);
    // printf("%d %d\n",gapp,gape);
    int bound = 0;
    pl[0] = 0;
    int row = 0;
    int bound_list[CORE_P + CORE_E];
    for(int i = 0;i < CORE_E + CORE_P;i++){
        if(i < CORE_P) bound += gapp;
        else bound += gape;
        bound_list[i] = bound;
    }
    for(int i = 0;i < CORE_E + CORE_P;i++){
        bound = bound_list[i];
        if(i < CORE_P) row = binary_search_cost(cost_sum,0,bound_mid,bound);//find row
        else row = binary_search_cost(cost_sum,bound_mid,m,bound);//find row
        // printf("%d %d %d ",row,cost_sum[row],bound);
        pr[i] = pl[i + 1] = row;
        nnz_l[i] = nnz_r[i - 1];
        nnz_r[i] = search_nnz(csrColIdx,csrRowPtr[row],csrRowPtr[row+1],bound - cost_sum[row - 1]);
    }
    nnz_l[0] = 0;
    pr[CORE_P + CORE_E - 1] = m;
    nnz_r[CORE_P + CORE_E - 1] = nnz;
    gettimeofday(&t2,NULL);
    // for(int i = 0;i < CORE_E + CORE_P;i++){
    //     printf("%d %d %d %d\n",nnz_l[i],nnz_r[i],pl[i],pr[i]);
    // }

    MAT_VAL_TYPE *extra_y = (MAT_VAL_TYPE*)malloc(sizeof(MAT_VAL_TYPE) * (CORE_E + CORE_P));
    double time_sum = 0;
    double time_min = 0x7ffffff;
    double time;
    double pre_time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

    //reference
    for (int i = 0; i < m; i++)
    {
        y_ref[i] = 0;
        for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
        {
            y_ref[i] += csrVal[j] * x[csrColIdx[j]];
        }
    }

    //warm
    for(int pp = 0;pp < REPEAT;pp++){
        #pragma omp parallel for
        for(int id = 0;id < CORE_E + CORE_P;id++){
            if(id < CORE_P){
                if(pl[id] == pr[id]){
                    extra_y[id] = avx2_kernel_p(nnz_l[id],nnz_r[id],csrColIdx,csrVal,x);
                    continue;
                }
                y[pl[id]] = avx2_kernel_p(nnz_l[id],csrRowPtr[pl[id] + 1],csrColIdx,csrVal,x);
                for(int i = pl[id] + 1;i < pr[id];i++){
                    y[i] = avx2_kernel_p(csrRowPtr[i],csrRowPtr[i + 1],csrColIdx,csrVal,x);
                }
                extra_y[id] = avx2_kernel_p(csrRowPtr[pr[id]],nnz_r[id],csrColIdx,csrVal,x);
            }
            else{
                if(pl[id] == pr[id]){
                    extra_y[id] = avx2_kernel_e(nnz_l[id],nnz_r[id],csrColIdx,csrVal,x);
                    continue;
                }
                y[pl[id]] = avx2_kernel_e(nnz_l[id],csrRowPtr[pl[id] + 1],csrColIdx,csrVal,x);
                for(int i = pl[id] + 1;i < pr[id];i++){
                    y[i] = avx2_kernel_e(csrRowPtr[i],csrRowPtr[i + 1],csrColIdx,csrVal,x);
                }
                extra_y[id] = avx2_kernel_e(csrRowPtr[pr[id]],nnz_r[id],csrColIdx,csrVal,x);
            }
        } 
        #pragma omp barrier
    }
    for(int pp = 0;pp < REPEAT;pp++){
        gettimeofday(&t1,NULL);
        #pragma omp parallel for
        for(int id = 0;id < CORE_P + CORE_E;id++){
            if(id < CORE_P){
                if(pl[id] == pr[id]){
                    extra_y[id] = avx2_kernel_p(nnz_l[id],nnz_r[id],csrColIdx,csrVal,x);
                    continue;
                }
                y[pl[id]] = avx2_kernel_p(nnz_l[id],csrRowPtr[pl[id] + 1],csrColIdx,csrVal,x);
                for(int i = pl[id] + 1;i < pr[id];i++){
                    y[i] = avx2_kernel_p(csrRowPtr[i],csrRowPtr[i + 1],csrColIdx,csrVal,x);
                }
                extra_y[id] = avx2_kernel_p(csrRowPtr[pr[id]],nnz_r[id],csrColIdx,csrVal,x);
            }
            else{
                if(pl[id] == pr[id]){
                    extra_y[id] = avx2_kernel_e(nnz_l[id],nnz_r[id],csrColIdx,csrVal,x);
                    continue;
                }
                y[pl[id]] = avx2_kernel_e(nnz_l[id],csrRowPtr[pl[id] + 1],csrColIdx,csrVal,x);
                for(int i = pl[id] + 1;i < pr[id];i++){
                    y[i] = avx2_kernel_e(csrRowPtr[i],csrRowPtr[i + 1],csrColIdx,csrVal,x);
                }
                extra_y[id] = avx2_kernel_e(csrRowPtr[pr[id]],nnz_r[id],csrColIdx,csrVal,x);
            }
        }
        
        #pragma omp barrier
        gettimeofday(&t2, NULL);
        time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        time_sum += time;
        time_min = time_min > time ? time : time_min;
    }

    gettimeofday(&t3, NULL);
    // #pragma omp parallel for
    for(int id = 0;id < CORE_E + CORE_P;id++){
        y[pl[id]] += extra_y[id - 1];
    }
    gettimeofday(&t4, NULL);

    // printf("4642: %d %d\n",csrRowPtr[row_map[4642]],csrRowPtr[row_map[4642] + 1]);

    double time1 = (t4.tv_sec - t3.tv_sec) * 1000.0 + (t4.tv_usec - t3.tv_usec) / 1000.0;
    time_sum /= REPEAT;
    time_sum += time1;
    time_min += time1;
    double GFlops = (nnz * 2) / (time_min * 1000000);
    double GFlops_avg = (nnz * 2) / (time_sum * 1000000);

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
    
    printf("loop,%s,%.3lf,%.3lf\n",name,GFlops,GFlops_avg);
    // printf("%s,%d,%.3lf,%.3lf\n",name,part/10,GFlops,GFlops_avg);
    // printf("%s,%.3lf,%.3lf\n",name,pre_time,time_sum);
    // printf("by_avx2_zhan_max_%d_%d,%d,%.3f\n",part,100 - part,nnz,GFlops);
    // printf("by_avx2_zhan_avg_%d_%d,%d,%.3f\n",part,100 - part,nnz,GFlops_avg);
    free(name);

  // validate x
    {
    double accuracy = 0;
    double ref = 0.0;
    double res = 0.0;

    for (int i = 0; i < m; i++)
    {
        ref += abs(y_ref[i]);
        res += abs(y[i] - y_ref[i]);
        if(res != 0){
            printf("i %d y[i] %lf y_ref[i] %lf\n",i,y[i],y_ref[i]);
            break;
        }
        // printf("index = %d  %8.1lf %8.1lf\n", i, x[i], x_ref[i]);
    }

    int flag = 0;
    res = ref == 0 ? res : res / ref;
    
    if (res == accuracy)
    {
        flag = 1;
        // printf("y check passed! |vec-vecref|/|vecref| = %8.2e\n", res);
    }
    else
        printf("y check _NOT_ passed! |vec-vecref|/|vecref| = %8.2e\n", res);

    free(y);
    free(y_ref);
    free(x);
    free(extra_y);
    free(cost_sum);
    free(row_begin_nnz);

    free(csrRowPtr);
    free(csrColIdx);
    free(csrVal);
    free(old_csrRowPtr);
    free(old_csrColIdx);
    free(old_csrVal);

    return 1;
    }
}