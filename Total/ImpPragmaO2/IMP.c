#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <omp.h>





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

        }
    }
}



#pragma isolated_call(matTransposeImp)
int checkSymImp(int n, float(*mat)[n]) {  
    int ret_val = 1;
    for (int i = 0; i < n; ++i) {
#pragma unroll
        for (int j = 0; j < n; ++j) { 
#pragma execution_frequency(very_high)
            if (mat[i][j] != mat[j][i]) {

                ret_val = 0; // Non è simmetrica
            }
        }
    }
    return ret_val;
}




#pragma isolated_call(matTransposeImp)
void matTransposeImp(int n, float (*mat)[n], float (*tam)[n]) {  
    for (int r = 0; r < n; ++r){
#pragma unroll
        for (int c = 0; c < n; ++c) {
#pragma execution_frequency(very_high)
          tam[c][r] = mat[r][c];
         }
    }
}


int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Uso: %s <matrix dim>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);


    srand(time(NULL));

    FILE *file = fopen("outputPO2.csv", "a");
    if (file == NULL) {
        perror("Errore nell'apertura del file");
        return 1;
    }




    float(*mat)[n];
    mat = (float(*)[n])malloc(sizeof(*mat) * n);

    float(*tam)[n];
    tam = (float(*)[n])malloc(sizeof(*tam) * n);

    if (!mat || !tam) {
        printf( "failed to allocate mat and/or tam\n");
        return EXIT_FAILURE;
    }



    init_mat(n, mat);     // matrix init




    struct timeval start, end;
    float checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    double start_time, end_time, start_time2, end_time2;
    double checkSymOmpTime, matTransposeOmpTime;

    gettimeofday(&start, NULL);
#pragma inline
    int r = checkSymImp(n, mat);
    gettimeofday(&end, NULL);
    checkSymTimeImp = time_diff(&start, &end);

    if(r) {
        printf( "\n is sym\n");
    }else {
        printf("\n isn't sym\n");
    }

    gettimeofday(&start, NULL);
#pragma inline
    matTransposeImp(n, mat, tam);
    gettimeofday(&end, NULL);
    matTransposeTimeImp = time_diff(&start, &end);




    check(n, mat, tam); //check the correctness of the transposition

    printf("\nN: %d\n", n);
    printf( "checkSymImp [s]: %f\n", checkSymTimeImp);
    printf( "matTranposeImp [s]: %f\n", matTransposeTimeImp);

    fprintf(file, "1,%d,%f,%f,IMP\n", n, checkSymTimeImp, matTransposeTimeImp);

    free(tam);
    free(mat);



    return EXIT_SUCCESS;


}

