#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

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

//    int rows_per_proc = n / size;
//    int start_row = myrank * rows_per_proc;
//    int end_row = start_row + rows_per_proc;
//
//    if (myrank == size - 1) {
//        end_row = n;
//    }

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
    // Calcola quante righe gestirà ogni processo
    int rows_per_proc = N / size;
    int extra_rows = N % size;

    // Calcola l'indice di inizio e fine per questo processo
    int start_row = rank * rows_per_proc + (rank < extra_rows ? rank : extra_rows);
    int num_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);

    if (rank == 0) {
        // Il processo 0 elabora prima le proprie righe
        for (int i = start_row; i < start_row + num_rows; i++) {
            for (int j = 0; j < N; j++) {
                tam[j][i] = mat[i][j];
            }
        }

        // Poi riceve le righe dagli altri processi direttamente in tam
        for (int p = 1; p < size; p++) {
            int other_rows_per_proc = N / size;
            int other_extra = (p < extra_rows ? 1 : 0);
            int other_start = p * rows_per_proc + (p < extra_rows ? p : extra_rows);
            int other_num_rows = other_rows_per_proc + other_extra;

            for (int i = other_start; i < other_start + other_num_rows; i++) {
                // Riceve direttamente in tam[j][i] per ogni j
                MPI_Recv(&tam[0][i], N, MPI_FLOAT, p, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    } else {
        // Gli altri processi traspongono e inviano le proprie righe al processo 0
        float *temp_row = (float *)malloc(N * sizeof(float));

        for (int i = start_row; i < start_row + num_rows; i++) {
            // Traspone la riga nel buffer temporaneo
            for (int j = 0; j < N; j++) {
                temp_row[j] = mat[i][j];
            }
            MPI_Send(temp_row, N, MPI_FLOAT, 0, i, MPI_COMM_WORLD);
        }

        free(temp_row);
    }
}

void MatTransposeMPI_DerivedTypes(int N, float (*mat)[N], float (*tam)[N], int rank, int size) {
    // Calculate number of rows per process
    int base_rows_per_proc = N / size;
    int remainder_rows = N % size;
    int rows_per_proc = base_rows_per_proc + (rank < remainder_rows ? 1 : 0);
    int start_row = rank * base_rows_per_proc + MIN(rank, remainder_rows);

    // Create derived datatype for the block of rows this process will handle
    MPI_Datatype block_row_type;
    MPI_Type_contiguous(N * rows_per_proc, MPI_FLOAT, &block_row_type);
    MPI_Type_commit(&block_row_type);

    // Create derived datatype for columns
    MPI_Datatype column_type;
    MPI_Type_vector(N, 1, N, MPI_FLOAT, &column_type);
    MPI_Type_commit(&column_type);

    // Create a resized column type to handle proper striding
    MPI_Datatype resized_column_type;
    MPI_Type_create_resized(column_type, 0, sizeof(float), &resized_column_type);
    MPI_Type_commit(&resized_column_type);

    // Allocate temporary buffer for storing the local block
    float *local_block = (float *)malloc(N * rows_per_proc * sizeof(float));

    // Scatter the matrix rows to all processes
    MPI_Scatter(mat[0], 1, block_row_type,
                local_block, N * rows_per_proc, MPI_FLOAT,
                0, MPI_COMM_WORLD);

    // Transpose local block
    float *transposed_block = (float *)malloc(N * rows_per_proc * sizeof(float));
    for (int i = 0; i < rows_per_proc; i++) {
        for (int j = 0; j < N; j++) {
            transposed_block[j * rows_per_proc + i] = local_block[i * N + j];
        }
    }

    // Calculate send counts and displacements for gathering the transposed data
    int *sendcounts = (int *)malloc(size * sizeof(int));
    int *displs = (int *)malloc(size * sizeof(int));

    for (int i = 0; i < size; i++) {
        sendcounts[i] = base_rows_per_proc + (i < remainder_rows ? 1 : 0);
        displs[i] = i * base_rows_per_proc + min(i, remainder_rows);
    }

    // Gather the transposed blocks back to the root
    MPI_Gatherv(transposed_block, N * rows_per_proc, MPI_FLOAT,
                tam[0], sendcounts, displs, resized_column_type,
                0, MPI_COMM_WORLD);

    // Clean up
    MPI_Type_free(&block_row_type);
    MPI_Type_free(&column_type);
    MPI_Type_free(&resized_column_type);
    free(local_block);
    free(transposed_block);
    free(sendcounts);
    free(displs);
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

    if(size>n){
      	printf("This implementation requires the number of the row at most equal at the number of the number of processors\n");
    	MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (myrank == 0) {
      init_mat(n, mat);
    }


    //printf("Rank: %d\n", myrank);
//    for (int i = 0; i < n; i++) {
//      for (int j = 0; j < n; j++) {
//        printf("%f ", mat[i][j]);
//      }
//      printf("\n");
//    }
//    printf("\n");

    start_time_sym = MPI_Wtime();
    MPI_Bcast(mat, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    checkSymMPI(n, mat, myrank, size, sym);
    end_time_sym = MPI_Wtime();
    total_time_sym = end_time_sym - start_time_sym;

    if (myrank == 0) {
//      printf("sym: %d\n", *sym);
        if (*sym) {
            printf("is sym\n");
        } else {
            printf("is not sym\n");
        }
    }

    start_time = MPI_Wtime();
    MPI_Bcast(mat, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MatTransposeMPI_DerivedTypes(n, mat, tam, myrank, size);
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

