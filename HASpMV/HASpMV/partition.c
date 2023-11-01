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
    read_Dmatrix_32(&m, &n, &nnz, &csrRowPtr, &csrColIdx, &csrVal, &ls, filename);
    // printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnz);
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

     for(int pp = 0;pp < REPEAT;pp++){
        #pragma omp parallel for
        for(int id = 0;id < CORE_E + CORE_P;id++){
            // forceAffinity(id);
                extra_y[id] = 0;
                if(pl[id] == pr[id]){
                    for(int j = nnz_l[id];j < nnz_r[id];j++){
                        extra_y[id] += csrVal[j] * x[csrColIdx[j]];
                    }
                    continue;
                }
                y[pl[id]] = 0;
                for (int j = nnz_l[id]; j < csrRowPtr[pl[id] + 1]; j++)
                {
                    y[pl[id]] += csrVal[j] * x[csrColIdx[j]];
                }
                for (int i = pl[id] + 1; i < pr[id]; i++){
                    y[i] = 0;
                    for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
                    {
                        y[i] += csrVal[j] * x[csrColIdx[j]];
                    }
                }
                for (int j = csrRowPtr[pr[id]]; j < nnz_r[id]; j++)
                {
                    extra_y[id] += csrVal[j] * x[csrColIdx[j]];
                }
        }
        #pragma omp barrier
    }
    for(int pp = 0;pp < REPEAT;pp++){
        gettimeofday(&t1,NULL);
        #pragma omp parallel for
        for(int id = 0;id < CORE_E + CORE_P;id++){
            // forceAffinity(id);
                extra_y[id] = 0;
                if(pl[id] == pr[id]){
                    for(int j = nnz_l[id];j < nnz_r[id];j++){
                        extra_y[id] += csrVal[j] * x[csrColIdx[j]];
                    }
                    continue;
                }
                y[pl[id]] = 0;
                for (int j = nnz_l[id]; j < csrRowPtr[pl[id] + 1]; j++)
                {
                    y[pl[id]] += csrVal[j] * x[csrColIdx[j]];
                }
                for (int i = pl[id] + 1; i < pr[id]; i++){
                    y[i] = 0;
                    for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
                    {
                        y[i] += csrVal[j] * x[csrColIdx[j]];
                    }
                }
                for (int j = csrRowPtr[pr[id]]; j < nnz_r[id]; j++)
                {
                    extra_y[id] += csrVal[j] * x[csrColIdx[j]];
                }
        }
        #pragma omp barrier
        gettimeofday(&t2, NULL);
        time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        // printf("%lf\n",time);
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
    
    printf("partition,%s,%.3lf,%.3lf\n",name,GFlops,GFlops_avg);
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