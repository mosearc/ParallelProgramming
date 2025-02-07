#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <math.h>

void check(int n, float (*mat)[n], float(*tam)[n]) {
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) if (mat[r][c] != tam[c][r]) {
            printf("Error in matrix Trasposition\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

    printf("Matrix Traspositions success\n");

}

void init_mat(int n, float(* mat)[n] ) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            mat[i][j] = rand() / (float) RAND_MAX;
			//printf("%f ", mat[i][j]);
        }
        //printf("\n");
    }
    //printf("\n");
}

void checkSymMPI(int n, float(* mat)[n], int myrank, int size, int *tot) {
    int ret_val = 1;


    int base_rows_per_proc = n / size;
	int remainder_rows = n % size;
 	int rows_per_proc = base_rows_per_proc + (myrank < remainder_rows ? 1 : 0);

    int start_row = myrank * rows_per_proc;
    int end_row = start_row + rows_per_proc;


    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < n; j++) {
            if (mat[i][j] != mat[j][i]) { //isn't sym
                ret_val = 0;
            }
        }
    }

    MPI_Reduce(&ret_val, tot, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);
    //printf("Rank %d successfully completed with %d\n", myrank, *tot);
}



void MatTransposeMPI(int N, float (*mat)[N], float (*tam)[N], int rank, int size) {

    int rows_per_proc = N / size;
    int extra_rows = N % size;


    int start_row = rank * rows_per_proc + (rank < extra_rows ? rank : extra_rows);
    int num_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);


    float *row_buffer = (float *)malloc(N * sizeof(float));

    if (rank == 0) {

        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < N; j++) {
                tam[j][start_row + i] = mat[start_row + i][j];
            }
        }


        for (int p = 1; p < size; p++) {
            int other_rows_per_proc = N / size;
            int other_extra = (p < extra_rows ? 1 : 0);
            int other_start = p * rows_per_proc + (p < extra_rows ? p : extra_rows);
            int other_num_rows = other_rows_per_proc + other_extra;

            for (int i = 0; i < other_num_rows; i++) {

                MPI_Recv(row_buffer, N, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int j = 0; j < N; j++) {
                    tam[j][other_start + i] = row_buffer[j];
                }
            }
        }
    } else {

        for (int i = 0; i < num_rows; i++) {

            for (int j = 0; j < N; j++) {
                row_buffer[j] = mat[start_row + i][j];
            }

            MPI_Send(row_buffer, N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        }
    }

    free(row_buffer);
}


int main(int argc, char** argv) {

    if (argc != 3) {
        printf("Uso: %s <matrix dim> <#processes>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int processes = atoi(argv[2]);

    srand(time(NULL));

    FILE *file = fopen("output.csv", "a");
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

    int myrank, size;

    double start_time, end_time, start_time_sym, end_time_sym, total_time_transp, total_time_sym;

    int *sym;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



    if (myrank == 0) {
      init_mat(n, mat);
    }



    start_time_sym = MPI_Wtime();
    MPI_Bcast(mat, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    checkSymMPI(n, mat, myrank, size, sym);
    end_time_sym = MPI_Wtime();
    total_time_sym = end_time_sym - start_time_sym;

    if (myrank == 0) {

        if (*sym) {
            printf("is sym\n");
        } else {
            printf("is not sym\n");
        }
    }

    start_time = MPI_Wtime();
    MPI_Bcast(mat, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MatTransposeMPI(n, mat, tam, myrank, size);
    end_time = MPI_Wtime();
    total_time_transp = end_time - start_time;


    if (myrank == 0){
   		check(n, mat, tam); //check the correctness of the transposition
        fprintf(file, "%d,%d,%f,%f,MPI B\n",processes, n, total_time_sym, total_time_transp);
    }


    MPI_Finalize();

    free(tam);
    free(mat);


    return EXIT_SUCCESS;

}
