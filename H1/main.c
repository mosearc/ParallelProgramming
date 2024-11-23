#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <omp.h>

//#define OUTPUT_FILE "output.txt"
#define BLOCK_SIZE_CACHE 64 //trova quello giusto per il cluster - (cahce cluster/float)^0.5 dato che i blocchi son quadrati : (32K/4)^0.5 = 90

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

#pragma isolated_call(checkSymImp)
int checkSymImp(int n, float(*mat)[n]) {  // see cache-oblivous alg or tiling
    int tmp = 1;
    for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
        for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {
            // Trasponi il blocco (i, j)
            for (int k = i; k < i + BLOCK_SIZE_CACHE && k < n; ++k) {
                //#pragma ivdep -> peggiora
                //#pragma omp simd -> qui non ha senso
#pragma unroll
                for (int l = j; l < j + BLOCK_SIZE_CACHE && l < n; ++l) {
#pragma execution_frequency(very_high)
                    // Prefetcha una riga successiva di mat, ad esempio 4 iterazioni in avanti
                    if (l + 4 < n) {
                        __builtin_prefetch(&mat[k][l + 4], 0, 1);  // Prefetch per lettura
                    }
                    if (l + 4 < n) __builtin_prefetch(&mat[l + 4][k], 0, 1);

                    //La distanza tra il dato prefetchato e l'elemento corrente (l + 4 e k + 4 in questo esempio) può variare.
                    //Distanze troppo grandi rischiano di sostituire in cache i dati ancora utili,
                    //mentre distanze troppo corte potrebbero non prefetchare con sufficiente anticipo
                    //in matrici piccole potrebbe addirittua peggiorare

                    if (mat[k][l] != mat[l][k]) tmp = 0;
                }
            }
        }
    }
    return tmp;
}

// #pragma isolated_call(checkSymImp)
// int checkSymImp(int n, float (*mat)[n]) {
//     int ret_val = 1;
//     for (int bi = 0; bi < n; bi += BLOCK_SIZE_CACHE) {
//         for (int bj = 0; bj < n; bj += BLOCK_SIZE_CACHE) {
//             // Controlla solo i blocchi superiori alla diagonale (inclusa la diagonale)
//             for (int i = bi; i < bi + BLOCK_SIZE_CACHE && i < n; ++i) {
//                 // //#pragma ivdep -> peggiora
//                 // //#pragma omp simd -> qui non ha senso
// #pragma unroll
//                 for (int j = (bj == bi ? i + 1 : bj); j < bj + BLOCK_SIZE_CACHE && j < n; ++j) {
// #pragma execution_frequency(very_high)
//                     if (j + 1 < n) __builtin_prefetch(&mat[i][j + 1], 0, 1); // Lettura successiva di mat[i][j+1]
//                     if (j + 1 < n) __builtin_prefetch(&mat[j + 1][i], 0, 1); // Lettura successiva di mat[j+1][i]
//
//                     if (mat[i][j] != mat[j][i]) {
//                         ret_val = 0; // Non è simmetrica
//                     }
//                 }
//             }
//         }
//     }
//     return ret_val;
// }

// int checkSymOMPCancel(int n, float (*mat)[n]) {//inusabile nel cluster
//     int ret = EXIT_SUCCESS;
// #pragma omp parallel for collapse(2)
//     for (int r = 0; r < n; ++r) {
//         for (int c = r + 1 ; c < n; ++c) {
//             if (mat[r][c] != mat[c][r]) {
// #pragma omp atomic read
//                 ret = EXIT_FAILURE;
// #pragma omp cancel for //only for openMP 4.0+
//             }
//         }
//         //#pragma omp cancellation point for //only for openMP 4.0+
//     }
//     return ret;
// }

int checkSymOMP(int n, float (*mat)[n]) { //casino avendo molti thread
//     int ret_val = 1;
// #pragma omp parallel for collapse(2)
//     for (int r = 0; r < n; ++r) {
//         for (int c = 0; c < n; ++c) {
//             if(mat[r][c] != mat[c][r]){
// #pragma omp atomic write
//                 ret_val = 0;
//         }
//     }
// }
//     return ret_val;

    int flag = 1;
    int temp_flag;
#pragma omp parallel for private(temp_flag)
    for (int r = 0; r < n; ++r) {
        {
#pragma omp flush(flag)
#pragma omp atomic read
            temp_flag = flag;
            if(temp_flag) {
#pragma omp simd
                for (int c = 0; c < n; ++c) {
                    if (mat[r][c] != mat[c][r]) {
#pragma omp flush
#pragma omp atomic write
                        flag = 0;
#pragma omp flush(flag)
                    }
                }
            }
        }
    }
    return flag;
}

#pragma option_override(matTranspose, "opt(level, 0)")
void matTransposeHalf(int n, float (*mat)[n], float (*tam)[n]) { //inefficiente la imp per gli accessi distanti, questa sarebbe forte comunque
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            tam[i][j] = mat[j][i];
            tam[j][i] = mat[i][j];
        }
    }
}

#pragma option_override(matTranspose, "opt(level, 0)")
void matTranspose(int n, float (*mat)[n], float (*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
}

#pragma isolated_call(matTransposeImp)
void matTransposeImp(int n, float (*mat)[n], float (*tam)[n]) {  // see cache-oblivous alg or tiling
    for (int i = 0; i < n; i += BLOCK_SIZE_CACHE) {
        for (int j = 0; j < n; j += BLOCK_SIZE_CACHE) {
            // Trasponi il blocco (i, j)
            for (int k = i; k < i + BLOCK_SIZE_CACHE && k < n; ++k) {
#pragma unroll
                for (int l = j; l < j + BLOCK_SIZE_CACHE && l < n; ++l) {
#pragma execution_frequency(very_high)

                    // Prefetcha una riga successiva di mat, ad esempio 4 iterazioni in avanti
                    if (l + 4 < n) {
                        __builtin_prefetch(&mat[k][l + 4], 0, 1);  // Prefetch per lettura
                    }
                    // Prefetcha una posizione di tam per la scrittura
                    if (k + 4 < n) {
                        __builtin_prefetch(&tam[l][k + 4], 1, 1);  // Prefetch per scrittura
                    }

                    //La distanza tra il dato prefetchato e l'elemento corrente (l + 4 e k + 4 in questo esempio) può variare.
                    //Distanze troppo grandi rischiano di sostituire in cache i dati ancora utili,
                    //mentre distanze troppo corte potrebbero non prefetchare con sufficiente anticipo
                    //in matrici piccole potrebbe addirittua peggiorare
                    tam[l][k] = mat[k][l];
                }
            }
        }
    }

}

void matTransposeOMP(int n, float (*mat)[n], float (*tam)[n]) {
#pragma omp parallel for collapse(2)
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            tam[c][r] = mat[r][c];
        }
    }
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
    double start_time, end_time, start_time2, end_time2;
    double checkSymOmpTime, matTransposeOmpTime;


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
#pragma inline
    int r = checkSymImp(n, mat);
    //t = clock() - t;
    gettimeofday(&end, NULL);
    //checkSymTimeImp = ((double)t) / CLOCKS_PER_SEC;
    checkSymTimeImp = time_diff(&start, &end);

    // check simmetry
    //clock_t t = clock();
    start_time = omp_get_wtime();
    int rr = checkSymOMP(n, mat);
    //t = clock() - t;
    end_time = omp_get_wtime();
    //checkSymTime = ((double)t) / CLOCKS_PER_SEC;
    checkSymOmpTime = end_time - start_time;

    if(res) {
        printf( "\n è sym\n");
    }else {
        printf("\nnon è sym\n");
    }

    if(r) {
        printf( "\n è sym\n");
    }else {
        printf("\nnon è sym\n");
    }

    if(rr) {
        printf( "\n è sym\n");
    }else {
        printf("\nnon è sym\n");
    }




    // compute transpose
    //t = clock();
    gettimeofday(&start, NULL);
    matTranspose(n, mat, tam);
    //t = clock() - t;
    gettimeofday(&end, NULL);
    //matTransposeTime = ((double)t) / CLOCKS_PER_SEC;
    matTransposeTime = time_diff(&start, &end);

    check(n, mat, tam);

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
#pragma inline
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

    check(n, mat, tam);


    // compute transpose
    //t = clock();
    start_time2 = omp_get_wtime();
    matTransposeOMP(n, mat, tam);
    //t = clock() - t;
    end_time2 = omp_get_wtime();
    //matTransposeTime = ((double)t) / CLOCKS_PER_SEC;
    matTransposeOmpTime = end_time2 - start_time2;

    check(n, mat, tam);

    // fprintf(f, "\nTRANSPOSED\n");
    // print_mat(f, n, tam);

    // fprintf(f, "checkSym    [s]: %-10.10f\n", checkSymTime);
    // fprintf(f, "matTranpose [s]: %-10.10f\n", matTransposeTime);

    printf( "\ncheckSym    [s]: %-10.10f\n", checkSymTime);
    printf( "checkSymImp [s]: %-10.10f\n", checkSymTimeImp);
    printf( "checkSymOMP [s]: %-10.10f\n", checkSymOmpTime);
    printf( "matTranpose [s]: %-10.10f\n", matTransposeTime);
    printf( "matTranposeImp [s]: %-10.10f\n", matTransposeTimeImp);
    printf( "matTranposeOmpImp [s]: %-10.10f\n", matTransposeOmpTime);


    printf( "\ncheckSym    [s]: %0.8f\n", checkSymTime);
    printf( "checkSymImp [s]: %0.8f\n", checkSymTimeImp);
    printf( "checkSymOMP [s]: %0.8f\n", checkSymOmpTime);
    printf( "matTranpose [s]: %0.8f\n", matTransposeTime);
    printf( "matTranposeImp [s]: %0.8f\n", matTransposeTimeImp);
    printf( "matTranposeOmp [s]: %0.8f\n", matTransposeOmpTime);

    printf( "\ncheckSym    [s]: %f\n", checkSymTime);
    printf( "checkSymImp [s]: %f\n", checkSymTimeImp);
    printf( "checkSymOMP [s]: %f\n", checkSymOmpTime);
    printf( "matTranpose [s]: %f\n", matTransposeTime);
    printf( "matTranposeImp [s]: %f\n", matTransposeTimeImp);
    printf( "matTranposeOmp [s]: %f\n", matTransposeOmpTime);


    free(tam);
    free(mat);
    free(tamImp);
    //fclose(f);

    return EXIT_SUCCESS;
}