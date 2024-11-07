#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define OUTPUT_FILE "output.txt"
#define TILE 64

void print_mat(FILE* f, int n, double (*mat)[n]) {
    for (int r = 0; r < n; ++r) {
        fprintf(f, "[ ");
        for (int c = 0; c < n; ++c) fprintf(f, "%3f ", mat[r][c]);
        fprintf(f, "]\n");
    }
    fprintf(f, "\n");
}

void init_mat(int n, double(* mat)[n] ) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat[i][j] = rand() / (double) RAND_MAX;
        }
    }
}

int checkSym(int n, double (*mat)[n]) {
    for (int r = 1; r < n; ++r)
        for (int c = r; c < n; ++c)
            if (mat[r][c] != mat[c][r]) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int checkSymImp(int n, double (*mat)[n]) {  // see cache-oblivous alg or tiling
    // strided access
    for (int r = 1; r < n; ++r)
#pragma simd
        for (int c = r; c < n; ++c)
            if (mat[r][c] != mat[c][r]) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int checkSymOMP(int n, double (*mat)[n]) {
#pragma omp parallel for collapse(2)
    for (int r = 1; r < n; ++r)
        for (int c = r; c < n; ++c)
            if (mat[r][c] != mat[c][r]) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

int checkSymOMPtile(int n, double (*mat)[n]) {
#pragma omp parallel for
    for (int i = 0; i < n; i += TILE)
        for (int j = 0; j < n; j += TILE)

            for (int r = i; r < i + TILE && r < n; ++r)
                for (int c = j; c < j + TILE && c < n; ++c)
                    if (mat[r][c] != mat[c][r]) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}

void matTranspose(int n, double (*mat)[n], double (*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
}

void matTransposeImp(int n, double (*mat)[n], double (*tam)[n]) {  // see cache-oblivous alg or tiling
    if (checkSym(n, mat)) {
        // sequential access
        for (int r = 0; r < n; ++r) {
#pragma simd
            for (int c = 0; c < n; ++c) tam[r][c] = mat[r][c];
        }
    } else {
        // strided access
        for (int r = 0; r < n; ++r) {
#pragma simd
            for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
        }
    }
}

void matTransposeOMP(int n, double (*mat)[n], double (*tam)[n]) {
#pragma omp parallel for collapse(2)
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) tam[c][r] = mat[r][c];
}

void matTransposeOMPtile(int n, double (*mat)[n], double (*tam)[n]) {
#pragma omp parallel for
    for (int i = 0; i < n; i += TILE)
        for (int j = 0; j < n; j += TILE)

            for (int r = i; r < i + TILE && r < n; ++r)
                for (int c = j; c < j + TILE && c < n; ++c) tam[c][r] = mat[r][c];
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
    printf( "\nn: %-3d\n\n", n);

    // matrix init
    double(*mat)[n];
    mat = (double(*)[n])malloc(sizeof(*mat) * n);

    double(*tam)[n];
    tam = (double(*)[n])malloc(sizeof(*tam) * n);

    if (!mat || !tam) {
        //fprintf(stderr, "failed to allocate mat and/or tam\n");
        printf( "failed to allocate mat and/or tam\n");
        return EXIT_FAILURE;
    }

    // for (int r = 0; r < n; ++r)
    //     for (int c = 0; c < n; ++c) mat[r][c] = rand() % 1000;  // just for debugging purposes

    init_mat(n, mat); //init it

    //print it
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("[%f]", mat[i][j]);
        }
        printf("\n");
    }

    // fprintf(f, "\nMATRIX\n");
    // print_mat(f, n, mat);

    double checkSymTime, matTransposeTime;

    // check simmetry
    clock_t t = clock();
    int res = checkSym(n, mat);
    t = clock() - t;
    checkSymTime = ((double)t) / CLOCKS_PER_SEC;

    // compute transpose
    t = clock();
    matTranspose(n, mat, tam);
    t = clock() - t;
    matTransposeTime = ((double)t) / CLOCKS_PER_SEC;

    //print it transp
    printf("\n\n");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("[%f]", tam[i][j]);
        }
        printf("\n");
    }

    // fprintf(f, "\nTRANSPOSED\n");
    // print_mat(f, n, tam);

    // fprintf(f, "checkSym    [s]: %-10.10f\n", checkSymTime);
    // fprintf(f, "matTranpose [s]: %-10.10f\n", matTransposeTime);

    printf( "\ncheckSym    [s]: %-10.10f\n", checkSymTime);
    printf( "matTranpose [s]: %-10.10f\n", matTransposeTime);

    free(tam);
    free(mat);
    //fclose(f);

    return EXIT_SUCCESS;
}