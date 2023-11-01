#include <iostream>
#include <cmath>
#include <string>

#include "anonymouslib_avx2.h"

#include "mmio.h"

#include "biio.h"

using namespace std;

#ifndef VALUE_TYPE
#define VALUE_TYPE double
#endif

#ifndef NUM_RUN
#define NUM_RUN 1000
#endif


int call_anonymouslib(int m, int n, int nnzA,
                  int *csrRowPtrA, int *csrColIdxA, VALUE_TYPE *csrValA,
                  VALUE_TYPE *x, VALUE_TYPE *y, VALUE_TYPE alpha, char * name)
{
    int err = 0;

    memset(y, 0, sizeof(VALUE_TYPE) * m);

    double gb = getB<int, VALUE_TYPE>(m, nnzA);
    double gflop = getFLOP<int>(nnzA);

    anonymouslibHandle<int, unsigned int, VALUE_TYPE> A(m, n);
    err = A.inputCSR(nnzA, csrRowPtrA, csrColIdxA, csrValA);
    //cout << "inputCSR err = " << err << endl;

    err = A.setX(x); // you only need to do it once!
    //cout << "setX err = " << err << endl;

    VALUE_TYPE *y_bench = (VALUE_TYPE *)malloc(m * sizeof(VALUE_TYPE));

    int sigma = ANONYMOUSLIB_CSR5_SIGMA; //nnzA/(8*ANONYMOUSLIB_CSR5_OMEGA);
    A.setSigma(sigma);

    for (int i = 0; i < 5; i++)
    {
        err = A.asCSR5();
        err = A.asCSR();
    }

    anonymouslib_timer asCSR5_timer;
    asCSR5_timer.start();

    err = A.asCSR5();
    double pre = asCSR5_timer.stop();
    // cout<<pre<<",";
//    cout << "CSR->CSR5 time = " << asCSR5_timer.stop() << " ms." << endl;
    //cout << "asCSR5 err = " << err << endl;

    // check correctness by running 1 time
    err = A.spmv(alpha, y);
    //cout << "spmv err = " << err << endl;

    // warm up by running 50 times
    if (NUM_RUN)
    {
        double time_sum = 0;
        double time_min = 0x7ffffff;
        double time;
        for (int i = 0; i < 50; i++) {
            memset(y,0,sizeof(double )*m);
            err = A.spmv(alpha, y);
        }

        anonymouslib_timer CSR5Spmv_timer;

        for (int i = 0; i < NUM_RUN; i++) {
            CSR5Spmv_timer.start();
            err = A.spmv(alpha, y_bench);
            time = CSR5Spmv_timer.stop();
            time_sum += time;
            if(i == 0) time_min = time;
            else time_min = time_min > time ? time : time_min;
        }
        time_sum /= NUM_RUN;
        double GFlops = (nnzA * 2) / (time_min * 1000000);
        double GFlops_avg = (nnzA * 2) / (time_sum * 1000000);
        //intf("CSR5_max,%d,%.3f\n",nnzA,GFlops);
        printf("CSR5,%s,%d,%.3f,%.3lf,%.3f,%.3lf\n",name,nnzA,GFlops,GFlops_avg,pre,time_sum);

    }

    free(y_bench);

    A.destroy();

    return err;
}

int main(int argc, char ** argv)
{
    // report precision of floating-point
    //cout << "------------------------------------------------------" << endl;
    char  *precision;
    if (sizeof(VALUE_TYPE) == 4)
    {
        precision = "32-bit Single Precision";
    }
    else if (sizeof(VALUE_TYPE) == 8)
    {
        precision = "64-bit Double Precision";
    }
    else
    {
        cout << "Wrong precision. Program exit!" << endl;
        return 0;
    }

    //cout << "PRECISION = " << precision << endl;
    //cout << "------------------------------------------------------" << endl;

    int m, n, nnzA;
    int sl;
    int *csrRowPtrA;
    int *csrColIdxA;
    VALUE_TYPE *csrValA;

    //ex: ./spmv webbase-1M.mtx
    int argi = 1;

    char  *filename;
    if(argc > argi)
    {
        filename = argv[argi];
        argi++;
    }
    //cout << "--------------" << filename << "--------------" << endl;
	read_Dmatrix_32(&m, &n, &nnzA, &csrRowPtrA, &csrColIdxA, &csrValA, &sl, filename);
	
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
    filename = name;
	/*a lot of codes is deleted*/
    
    srand(time(NULL));

    // set csrValA to 1, easy for checking floating-point results
    for (int i = 0; i < nnzA; i++)
    {
        csrValA[i] = rand() % 10;
    }

    //cout << " ( " << m << ", " << n << " ) nnz = " << nnzA << endl;

    VALUE_TYPE *x = (VALUE_TYPE *)_mm_malloc(n * sizeof(VALUE_TYPE), ANONYMOUSLIB_X86_CACHELINE);
    for (int i = 0; i < n; i++)
        x[i] = rand() % 10;

    VALUE_TYPE *y = (VALUE_TYPE *)_mm_malloc(m * sizeof(VALUE_TYPE), ANONYMOUSLIB_X86_CACHELINE);
    VALUE_TYPE *y_ref = (VALUE_TYPE *)_mm_malloc(m * sizeof(VALUE_TYPE), ANONYMOUSLIB_X86_CACHELINE);

    double gb = getB<int, VALUE_TYPE>(m, nnzA);
    double gflop = getFLOP<int>(nnzA);

    VALUE_TYPE alpha = 1.0;

    // compute reference results on a cpu core
    anonymouslib_timer ref_timer;
    ref_timer.start();

    int ref_iter = 1;
    for (int iter = 0; iter < ref_iter; iter++)
    {
        for (int i = 0; i < m; i++)
        {
            VALUE_TYPE sum = 0;
            for (int j = csrRowPtrA[i]; j < csrRowPtrA[i+1]; j++)
                sum += x[csrColIdxA[j]] * csrValA[j] * alpha;
            y_ref[i] = sum;
        }
    }

    double ref_time = ref_timer.stop() / (double)ref_iter;
    double GFlops = (2 * nnzA) / (ref_time * 1000000);

    // launch compute
    call_anonymouslib(m, n, nnzA, csrRowPtrA, csrColIdxA, csrValA, x, y, alpha, name);

    // compare reference and anonymouslib results
    int error_count = 0;
    for (int i = 0; i < m; i++)
        if (abs(y_ref[i] - y[i]) > 0.01 * abs(y_ref[i]))
        {
            error_count++;

        }


    _mm_free(csrRowPtrA);
    _mm_free(csrColIdxA);
    _mm_free(csrValA);
    _mm_free(x);
    _mm_free(y);
    _mm_free(y_ref);

    return 0;
}

