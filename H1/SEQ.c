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

#pragma option_override(checkSym, "opt(level, 0)")
int checkSym(int n, float(*mat)[n]) {
    int ret_val = 1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) { // Controlla solo sopra la diagonale
            if (mat[i][j] != mat[j][i]) {
                ret_val = 0;; // Non è simmetrica
            }
        }
    }
    return ret_val;
}

#pragma option_override(matTranspose, "opt(level, 0)")
void matTranspose(int n, float (*mat)[n], float (*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
      printf("Uso: %s <matrix dim>\n", argv[0]);
      return 1;
    }
    int n = atoi(argv[1]);

    srand(time(NULL));

    FILE *file = fopen("output.csv", "a");
    if (file == NULL) {
        perror("Errore nell'apertura del file");
        return 1;
    }

    //fprintf(file, "Threads,dim,checkSym,matTranspose,SEQ\n");

    // n init

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


    // double checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    struct timeval start, end;
    float checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    double start_time, end_time, start_time2, end_time2;
    double checkSymOmpTime, matTransposeOmpTime;


    // check simmetry
    gettimeofday(&start, NULL);
    int res = checkSym(n, mat);
    gettimeofday(&end, NULL);
    checkSymTime = time_diff(&start, &end);

    if(res) {
        printf( "\n is symmetric\n");
    }else {
        printf("\n isn't symmetirc\n");
    }

    // compute transpose
    gettimeofday(&start, NULL);
    matTranspose(n, mat, tam);
    gettimeofday(&end, NULL);
    matTransposeTime = time_diff(&start, &end);

    check(n, mat, tam);

    printf("\nN: %d\n", n);
    printf( "checkSym    [s]: %f\n", checkSymTime);
    printf( "matTranpose [s]: %f\n", matTransposeTime);

    fprintf(file, "1,%d,%f,%f,SEQ\n", n, checkSymTime, matTransposeTime);

    free(tam);
    free(mat);

    return EXIT_SUCCESS;

}