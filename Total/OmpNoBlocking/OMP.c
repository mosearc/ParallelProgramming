#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <omp.h>





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

int checkSymOMP(int n, float (*mat)[n]) { 

    int ret_val = 1;
#pragma omp parallel for collapse(2) schedule(static)
    for (int r = 0; r < n; ++r) {
      for (int c = 0; c < n; ++c) {
        if(mat[r][c] != mat[c][r]){
#pragma omp atomic write
          ret_val = 0;
        }
      }
    }
    return ret_val;

}

void matTransposeOMP(int n, float (*mat)[n], float (*tam)[n]) {

#pragma omp parallel for simd collapse(2) schedule(static) aligned(mat, tam : 64)
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            tam[c][r] = mat[r][c];
        }
    }

}



int main(int argc, char *argv[]) {
      if (argc != 3) {
        printf("Uso: %s <matrix dim> <#threads>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int threads = atoi(argv[2]);

    srand(time(NULL));

    FILE *file = fopen("outputNB.csv", "a");
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

    init_mat(n, mat); // matrix init

    struct timeval start, end;
    float checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    double start_time, end_time, start_time2, end_time2;
    double checkSymOmpTime, matTransposeOmpTime;



    start_time = omp_get_wtime();

    int rr = checkSymOMP(n, mat);
    end_time = omp_get_wtime();

    checkSymOmpTime = end_time - start_time;


    if(rr) {
        printf( "\n is sym\n");
    }else {
        printf("\n isn't sym\n");
    }

    start_time2 = omp_get_wtime();

    matTransposeOMP(n, mat, tam);
    end_time2 = omp_get_wtime();

    matTransposeOmpTime = end_time2 - start_time2;



    check(n, mat, tam); //check the correctness of the transposition

    printf("\nN: %d\n", n);
    printf( "checkSymOMP [s]: %f\n", checkSymOmpTime);
    printf( "matTranposeOmp [s]: %f\n", matTransposeOmpTime);

    fprintf(file, "%d,%d,%f,%f,OMP\n",threads, n, checkSymOmpTime, matTransposeOmpTime);

    free(tam);
    free(mat);


    return EXIT_SUCCESS;


}
