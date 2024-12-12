#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <omp.h>
#include <sys/param.h>





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

#pragma option_override(checkSym, "opt(level, 0)")
int checkSym(int n, float(*mat)[n]) {
    int ret_val = 1;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (mat[i][j] != mat[j][i]) {
                ret_val = 0;;
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

void matTransposeBlock(int N, float (*mat)[N], float (*tam)[N], int block_size) {
    // Make sure block_size divides N evenly
    // If not, we'll need to handle edge cases
    if (N % block_size != 0) {
        printf("Warning: Matrix size %d is not divisible by block size %d. Performance may be suboptimal.\n",
               N, block_size);
    }

    // Process the matrix in blocks
    for (int i = 0; i < N; i += block_size) {
        for (int j = 0; j < N; j += block_size) {
            // Process each block
            for (int bi = i; bi < MIN(i + block_size, N); bi++) {
                for (int bj = j; bj < MIN(j + block_size, N); bj++) {
                    tam[bj][bi] = mat[bi][bj];
                }
            }
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

    FILE *file = fopen("outputBvNB.csv", "a");
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

    init_mat(n, mat); //matrix init



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
    //matTranspose(n, mat, tam);
    matTransposeBlock(n, mat, tam, 2);
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
