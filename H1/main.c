#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>

#define OUTPUT_FILE "output.txt"
#define BLOCK_SIZE_CACHE 64 //trova quello giusto per il cluster
#define BLOCK_SIZE_SMID 4

// void print_mat(FILE* f, int n, double (*mat)[n]) {
//     for (int r = 0; r < n; ++r) {
//         fprintf(f, "[ ");
//         for (int c = 0; c < n; ++c) fprintf(f, "%3f ", mat[r][c]);
//         fprintf(f, "]\n");
//     }
//     fprintf(f, "\n");
// }

float time_diff(struct timeval *start, struct timeval *end) {
    return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}

void check(int n, float (*mat)[n], float(*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) if (mat[r][c] != tam[r][c]) {
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

int checkSym(int n, float(*mat)[n]) {
    for (int r = 1; r < n; ++r)
        for (int c = r; c < n; ++c)
            if (mat[r][c] != mat[c][r]) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int checkSymImp(int n, float(*mat)[n]) {  // see cache-oblivous alg or tiling
    for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
        for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {
            // Trasponi il blocco (i, j)
            for (int k = i; k < i + BLOCK_SIZE_CACHE && k < n; ++k) {
                for (int l = j; l < j + BLOCK_SIZE_CACHE && l < n; ++l) {

                    // // Prefetcha una riga successiva di mat, ad esempio 4 iterazioni in avanti
                    // if (l + 4 < n) {
                    //     __builtin_prefetch(&mat[k][l + 4], 0, 1);  // Prefetch per lettura
                    // }
                    // // Prefetcha una posizione di tam per la scrittura
                    // if (k + 4 < n) {
                    //     __builtin_prefetch(&tam[l][k + 4], 1, 1);  // Prefetch per scrittura
                    // }

                    //La distanza tra il dato prefetchato e l'elemento corrente (l + 4 e k + 4 in questo esempio) può variare.
                    //Distanze troppo grandi rischiano di sostituire in cache i dati ancora utili,
                    //mentre distanze troppo corte potrebbero non prefetchare con sufficiente anticipo
                    //in matrici piccole potrebbe addirittua peggiorare

                    if (mat[k][l] != mat[l][k]) return EXIT_FAILURE;
                }
            }
        }
    }
    return EXIT_SUCCESS;
}

int checkSymOMP(int n, float (*mat)[n]) {
#pragma omp parallel for collapse(2)
    for (int r = 1; r < n; ++r)
        for (int c = r; c < n; ++c)
            if (mat[r][c] != mat[c][r]) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}


void matTranspose(int n, float (*mat)[n], float (*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
}

void matTransposeImp(int n, float (*mat)[n], float (*tam)[n]) {  // see cache-oblivous alg or tiling
    for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
        for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {
            // Trasponi il blocco (i, j)
            for (int k = i; k < i + BLOCK_SIZE_CACHE && k < n; ++k) {
                for (int l = j; l < j + BLOCK_SIZE_CACHE && l < n; ++l) {

                    // // Prefetcha una riga successiva di mat, ad esempio 4 iterazioni in avanti
                    // if (l + 4 < n) {
                    //     __builtin_prefetch(&mat[k][l + 4], 0, 1);  // Prefetch per lettura
                    // }
                    // // Prefetcha una posizione di tam per la scrittura
                    // if (k + 4 < n) {
                    //     __builtin_prefetch(&tam[l][k + 4], 1, 1);  // Prefetch per scrittura
                    // }

                    //La distanza tra il dato prefetchato e l'elemento corrente (l + 4 e k + 4 in questo esempio) può variare.
                    //Distanze troppo grandi rischiano di sostituire in cache i dati ancora utili,
                    //mentre distanze troppo corte potrebbero non prefetchare con sufficiente anticipo
                    //in matrici piccole potrebbe addirittua peggiorare

                    tam[l][k] = mat[k][l];
                }
            }
        }
    }

    /*

//gcc -O3 -mavx -o transpose_avx matrix_transpose.c OCCHIO CHE È TUTTO CON I DOUBLE
    for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
        for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {
            // Trasponi un blocco 4x4 usando AVX
            for (int k = i; k < i + BLOCK_SIZE_CACHE&& k < n; k += 4) {
                for (int l = j; l < j + BLOCK_SIZE_CACHE && l < n; l += 4) {


                    //prefetch della prossima riga del blocco
                    __builtin_prefetch(&mat[k + 4][l], 0, 1);
                    __builtin_prefetch(&mat[k][l + 4], 0, 1);
                    __builtin_prefetch(&tam[l][k], 1, 1);



                    // Carica blocchi da 4 double alla volta
                    __m256d row0 = _mm256_loadu_pd(&mat[k][l]);
                    __m256d row1 = _mm256_loadu_pd(&mat[k+1][l]);
                    __m256d row2 = _mm256_loadu_pd(&mat[k+2][l]);
                    __m256d row3 = _mm256_loadu_pd(&mat[k+3][l]);

                    // Trasposizione 4x4
                    __m256d t0 = _mm256_unpacklo_pd(row0, row1);
                    __m256d t1 = _mm256_unpackhi_pd(row0, row1);
                    __m256d t2 = _mm256_unpacklo_pd(row2, row3);
                    __m256d t3 = _mm256_unpackhi_pd(row2, row3);

                    row0 = _mm256_permute2f128_pd(t0, t2, 0x20);
                    row1 = _mm256_permute2f128_pd(t1, t3, 0x20);
                    row2 = _mm256_permute2f128_pd(t0, t2, 0x31);
                    row3 = _mm256_permute2f128_pd(t1, t3, 0x31);

                    // Memorizza il blocco trasposto
                    _mm256_storeu_pd(&tam[l][k], row0);
                    _mm256_storeu_pd(&tam[l+1][k], row1);
                    _mm256_storeu_pd(&tam[l+2][k], row2);
                    _mm256_storeu_pd(&tam[l+3][k], row3);
                }
            }
        }
    }

    */

}

void matTransposeOMP(int n, float (*mat)[n], float (*tam)[n]) {
#pragma omp parallel for collapse(2)
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
}


int main() {
    srand(time(NULL));

    // FILE* f = fopen(OUTPUT_FILE, "a");
    // if (!f) {
    //     fprintf(stderr, "failed to open file %s\n", OUTPUT_FILE);
    //     return EXIT_FAILURE;
    // }

    // n init
    int n;

    do {
        printf("\nEnter n: ");
        scanf("%d", &n);
    } while (n < 1);
    //fprintf(f, "\nn: %-3d\n\n", n);
    //printf( "\nn: %-3d\n\n", n);

    // matrix init
    float(*mat)[n];
    mat = (float(*)[n])malloc(sizeof(*mat) * n);

    float(*tam)[n];
    tam = (float(*)[n])malloc(sizeof(*tam) * n);

    float(*tamImp)[n];
    tamImp = (float(*)[n])malloc(sizeof(*tamImp) * n);

    if (!mat || !tam) {
        //fprintf(stderr, "failed to allocate mat and/or tam\n");
        printf( "failed to allocate mat and/or tam\n");
        return EXIT_FAILURE;
    }

    // for (int r = 0; r < n; ++r)
    //     for (int c = 0; c < n; ++c) mat[r][c] = rand() % 1000;  // just for debugging purposes

    init_mat(n, mat); //init it

    //print it
    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         printf("[%f]", mat[i][j]);
    //     }
    //     printf("\n");
    // }

    // fprintf(f, "\nMATRIX\n");
    // print_mat(f, n, mat);

    // double checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    struct timeval start, end;
    float checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;


    // check simmetry
    //clock_t t = clock();
    gettimeofday(&start, NULL);
    int res = checkSym(n, mat);
    //t = clock() - t;
    gettimeofday(&end, NULL);
    //checkSymTime = ((double)t) / CLOCKS_PER_SEC;
    checkSymTime = time_diff(&start, &end);

    // check simmetry
    //t = clock();
    gettimeofday(&start, NULL);
    int r = checkSymImp(n, mat);
    //t = clock() - t;
    gettimeofday(&end, NULL);
    //checkSymTimeImp = ((double)t) / CLOCKS_PER_SEC;
    checkSymTimeImp = time_diff(&start, &end);

    if(res == EXIT_SUCCESS) {
        printf( "\n\n è sym\n");
    }else {
        printf("\n\nnon è sym\n");
    }

    if(r == EXIT_SUCCESS) {
        printf( "\n\n è sym\n");
    }else {
        printf("\n\nnon è sym\n");
    }

    // compute transpose
    //t = clock();
    gettimeofday(&start, NULL);
    matTranspose(n, mat, tam);
    //t = clock() - t;
    gettimeofday(&end, NULL);
    //matTransposeTime = ((double)t) / CLOCKS_PER_SEC;
    matTransposeTime = time_diff(&start, &end);

    // printf("\n\n");
    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         printf("[%f]", tam[i][j]);
    //     }
    //     printf("\n");
    // }
    //
    // printf("\n ----- \n");

    //t = clock();
    gettimeofday(&start, NULL);
    matTransposeImp(n, mat, tamImp);
    //t = clock() - t;
    gettimeofday(&end, NULL);
    //matTransposeTimeImp = ((double)t) / CLOCKS_PER_SEC;
    matTransposeTimeImp = time_diff(&start, &end);

    // //print it transp
    // printf("\n\n");
    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         printf("[%f]", tam[i][j]);
    //     }
    //     printf("\n");
    // }

    //tamImp[1][3] = 42;

    check(n, tam, tamImp);

    // fprintf(f, "\nTRANSPOSED\n");
    // print_mat(f, n, tam);

    // fprintf(f, "checkSym    [s]: %-10.10f\n", checkSymTime);
    // fprintf(f, "matTranpose [s]: %-10.10f\n", matTransposeTime);

    printf( "\ncheckSym    [s]: %-10.10f\n", checkSymTime);
    printf( "checkSymImp [s]: %-10.10f\n", checkSymTimeImp);
    printf( "matTranpose [s]: %-10.10f\n", matTransposeTime);
    printf( "matTranposeImp [s]: %-10.10f\n", matTransposeTimeImp);

    printf( "\ncheckSym    [s]: %0.8f\n", checkSymTime);
    printf( "checkSymImp [s]: %0.8f\n", checkSymTimeImp);
    printf( "matTranpose [s]: %0.8f\n", matTransposeTime);
    printf( "matTranposeImp [s]: %0.8f\n", matTransposeTimeImp);

    free(tam);
    free(mat);
    free(tamImp);
    //fclose(f);

    return EXIT_SUCCESS;
}