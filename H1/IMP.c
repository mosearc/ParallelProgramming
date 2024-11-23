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

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Uso: %s <numero>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);


    srand(time(NULL));

    FILE *file = fopen("output.csv", "a");
    if (file == NULL) {
        perror("Errore nell'apertura del file");
        return 1;
    }

    //fprintf(file, "Threads,dim,checkSym,matTranspose,IMP\n");

    // matrix init
    float(*mat)[n];
    mat = (float(*)[n])malloc(sizeof(*mat) * n);

    float(*tam)[n];
    tam = (float(*)[n])malloc(sizeof(*tam) * n);

    float(*tamImp)[n];
    tamImp = (float(*)[n])malloc(sizeof(*tamImp) * n);

    if (!mat || !tam) {
        printf( "failed to allocate mat and/or tam\n");
        return EXIT_FAILURE;
    }



    init_mat(n, mat); //init it

        //print it
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             printf("[%f]", mat[i][j]);
//         }
//         printf("\n");
//     }


    struct timeval start, end;
    float checkSymTime, matTransposeTime, matTransposeTimeImp, checkSymTimeImp;
    double start_time, end_time, start_time2, end_time2;
    double checkSymOmpTime, matTransposeOmpTime;

    gettimeofday(&start, NULL);
#pragma inline
    int r = checkSymImp(n, mat);
    gettimeofday(&end, NULL);
    checkSymTimeImp = time_diff(&start, &end);

//    if(r) {
//        printf( "\n è sym\n");
//    }else {
//        printf("\nnon è sym\n");
//    }

    gettimeofday(&start, NULL);
#pragma inline
    matTransposeImp(n, mat, tamImp);
    gettimeofday(&end, NULL);
    matTransposeTimeImp = time_diff(&start, &end);

//         for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             printf("[%f]", tam[i][j]);
//         }
//         printf("\n");
//     }


    //check(n, mat, tam);

    printf("\nN: %d\n", n);
    printf( "checkSymImp [s]: %f\n", checkSymTimeImp);
    printf( "matTranposeImp [s]: %f\n", matTransposeTimeImp);

    fprintf(file, "1,%d,%f,%f,IMP\n", n, checkSymTimeImp, matTransposeTimeImp);

    free(tam);
    free(mat);
    free(tamImp);



    return EXIT_SUCCESS;


}
