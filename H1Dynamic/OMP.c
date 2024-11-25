#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <omp.h>

#define BLOCK_SIZE_CACHE 64 //trova quello giusto per il cluster - (cahce cluster/float)^0.5 dato che i blocchi son quadrati : (32K/4)^0.5 = 90

float time_diff(struct timeval *start, struct timeval *end) {
    return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}

void check(int n, float (*mat)[n], float(*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) if (mat[r][c] != tam[c][r]) {
            printf("\nError in matrix Trasposition\n");
            abort();
        }

    printf("\nMatrix Traspositions success\n");

}

void init_mat(int n, float(* mat)[n] ) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat[i][j] = rand() / (float) RAND_MAX;
            //mat[i][j] = 42;
        }
    }
}

int checkSymOMP(int n, float (*mat)[n]) { //casino avendo molti thread a causa di atomic
    int ret_val = 1;
#pragma omp parallel proc_bind(close)
    {
        int local_ret_val = 1;

#pragma omp for collapse(2) schedule(guided) nowait
        for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
            for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {
                for (int k = i; k < i + BLOCK_SIZE_CACHE && k < n; ++k) {
#pragma omp prefetch
                    for (int l = j; l < j + BLOCK_SIZE_CACHE && l < n; ++l) {
                        if (mat[k][l] != mat[l][k]) {
                            local_ret_val = 0;
                        }
                    }
                }
            }
        }

        if (!local_ret_val) {
#pragma omp atomic write
            ret_val = 0;
        }
    }

    return ret_val;

//    int ret_val = 1;
//#pragma omp parallel for collapse(2) schedule(static)
//    for (int r = 0; r < n; ++r) {
//      for (int c = 0; c < n; ++c) {
//        if(mat[r][c] != mat[c][r]){
//#pragma omp atomic write
//          ret_val = 0;
//        }
//      }
//    }
//    return ret_val;

//    int flag = 1;
//    int temp_flag;
//#pragma omp parallel for private(temp_flag)
//    for (int r = 0; r < n; ++r) {
//        {
//#pragma omp flush(flag)
//#pragma omp atomic read
//            temp_flag = flag;
//            if(temp_flag) {
//#pragma omp simd
//                for (int c = 0; c < n; ++c) {
//                    if (mat[r][c] != mat[c][r]) {
//#pragma omp flush
//#pragma omp atomic write
//                        flag = 0;
//#pragma omp flush(flag)
//                    }
//                }
//            }
//        }
//    }
//    return flag;
}

void matTransposeOMP(int n, float (*mat)[n], float (*tam)[n]) {

#pragma omp parallel proc_bind(close)
    {
#pragma omp for collapse(2) schedule(dynamic) nowait
        for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
            for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {

                int k_end = (i + BLOCK_SIZE_CACHE < n) ? i + BLOCK_SIZE_CACHE : n;
                int l_end = (j + BLOCK_SIZE_CACHE < n) ? j + BLOCK_SIZE_CACHE : n;

#pragma omp simd aligned(mat,tam:16)
                for (int k = i; k < k_end; ++k) {
#pragma omp prefetch
                    for (int l = j; l < l_end; ++l) {
                        tam[l][k] = mat[k][l];
                    }
                }
            }
        }
    }

//#pragma omp parallel for collapse(2) schedule(dynamic)
//#pragma omp parallel for simd collapse(2) schedule(static) aligned(mat, tam : 16)
//    for (int r = 0; r < n; ++r) {
//        for (int c = 0; c < n; ++c) {
//            tam[c][r] = mat[r][c];
//        }
//    }

}



int main(int argc, char *argv[]) {
      if (argc != 3) {
        printf("Uso: %s <matrix dim> <#threads>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int threads = atoi(argv[2]);

    srand(time(NULL));

    FILE *file = fopen("output.csv", "a");
    if (file == NULL) {
        perror("Errore nell'apertura del file");
        return 1;
    }

    //fprintf(file, "Threads,dim,checkSym,matTranspose,OMP\n");



    // matrix init
    float(*mat)[n];
    mat = (float(*)[n])malloc(sizeof(*mat) * n);

    float(*tam)[n];
    tam = (float(*)[n])malloc(sizeof(*tam) * n);

    if (!mat || !tam) {
        printf( "failed to allocate mat and/or tam\n");
        return EXIT_FAILURE;
    }

    init_mat(n, mat); //init it

    struct timeval start, end;
    float checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    double start_time, end_time, start_time2, end_time2;
    double checkSymOmpTime, matTransposeOmpTime;
    //float checkSymOmpTime, matTransposeOmpTime;


    start_time = omp_get_wtime();
    //gettimeofday(&start, NULL);
    int rr = checkSymOMP(n, mat);
    end_time = omp_get_wtime();
    //gettimeofday(&end, NULL);
    checkSymOmpTime = end_time - start_time;
    //checkSymOmpTime = time_diff(&start, &end);

    if(rr) {
        printf( "\n is sym\n");
    }else {
        printf("\n isn't sym\n");
    }

    start_time2 = omp_get_wtime();
    //gettimeofday(&start, NULL);
    matTransposeOMP(n, mat, tam);
    end_time2 = omp_get_wtime();
    //gettimeofday(&end, NULL);
    matTransposeOmpTime = end_time2 - start_time2;
    //matTransposeOmpTime = time_diff(&start, &end);


    check(n, mat, tam); //check the correctness of the transposition

    printf("\nN: %d\n", n);
    printf( "checkSymOMP [s]: %f\n", checkSymOmpTime);
    printf( "matTranposeOmp [s]: %f\n", matTransposeOmpTime);

    fprintf(file, "%d,%d,%f,%f,OMP\n",threads, n, checkSymOmpTime, matTransposeOmpTime);

    free(tam);
    free(mat);


    return EXIT_SUCCESS;


}
