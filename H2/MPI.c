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

void MatTransposeMPI_RowSiz(int N, float (*mat)[N], float (*tam)[N], int rank, int size) {
    // Works only when N equals size and matrix can be evenly divided
    if (N != size) {
        if (rank == 0) {
            printf("This implementation requires N == size\n");
        }
        return;
    }

    // Each process will handle exactly one row of the matrix
    int rows_per_proc = 1;

    // Allocate local arrays
    float* local_matrix = (float*)malloc(N * sizeof(float));  // One row per process
    if (!local_matrix) {
        printf("Memory allocation failed for local_matrix on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Scatter the matrix - each process gets one row
    MPI_Scatter(mat, N, MPI_FLOAT,
                local_matrix, N, MPI_FLOAT,
                0, MPI_COMM_WORLD);

    // Allocate send buffer for transpose operation
    float* send_buffer = (float*)malloc(N * sizeof(float));
    if (!send_buffer) {
        printf("Memory allocation failed for send_buffer on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Copy and reorder data for sending
    for (int i = 0; i < N; i++) {
        send_buffer[i] = local_matrix[i];
    }

    // Allocate receive buffer for transposed data
    float* local_transpose = (float*)malloc(N * sizeof(float));
    if (!local_transpose) {
        printf("Memory allocation failed for local_transpose on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Perform all-to-all communication
    MPI_Alltoall(send_buffer, 1, MPI_FLOAT,
                 local_transpose, 1, MPI_FLOAT,
                 MPI_COMM_WORLD);

    // Gather the transposed matrix
    MPI_Gather(local_transpose, N, MPI_FLOAT,
               tam, N, MPI_FLOAT,
               0, MPI_COMM_WORLD);

    // Cleanup
    free(local_matrix);
    free(local_transpose);
    free(send_buffer);
}

void MatTransposeMPI(int N, float (*mat)[N], float (*tam)[N], int rank, int size) {

    // Calcola quante righe gestirà ogni processo
    int rows_per_proc = N / size;
    int extra_rows = N % size;

    // Calcola l'indice di inizio e fine per questo processo
    int start_row = rank * rows_per_proc + (rank < extra_rows ? rank : extra_rows);
    int num_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);

    // Buffer temporaneo per memorizzare una riga
    float *temp_row = (float *)malloc(N * sizeof(float));

    // Ogni processo elabora le proprie righe
    for (int i = start_row; i < start_row + num_rows; i++) {
        // Copia la riga corrente nel buffer temporaneo
        for (int j = 0; j < N; j++) {
            temp_row[j] = mat[i][j];
        }

        // Invia la riga a tutti gli altri processi
        for (int p = 0; p < size; p++) {
            if (p != rank) {
                MPI_Send(temp_row, N, MPI_FLOAT, p, i, MPI_COMM_WORLD);
            }
        }

        // Il processo corrente copia i suoi elementi nella posizione trasposta
        for (int j = 0; j < N; j++) {
            tam[j][i] = temp_row[j];
        }
    }

    // Ricevi le righe dagli altri processi
    for (int p = 0; p < size; p++) {
        if (p == rank) continue;

        int other_rows_per_proc = N / size;
        int other_extra = (p < extra_rows ? 1 : 0);
        int other_start = p * rows_per_proc + (p < extra_rows ? p : extra_rows);
        int other_num_rows = other_rows_per_proc + other_extra;

        for (int i = other_start; i < other_start + other_num_rows; i++) {
            MPI_Recv(temp_row, N, MPI_FLOAT, p, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < N; j++) {
                tam[j][i] = temp_row[j];
            }
        }
    }

    free(temp_row);

    // Sincronizza tutti i processi prima di continuare
    //MPI_Barrier(MPI_COMM_WORLD);
}

void MatTransposeMPI_opt(int N, float (*mat)[N], float (*tam)[N], int rank, int size) {
    // Calcola quante righe gestirà ogni processo
    int rows_per_proc = N / size;
    int extra_rows = N % size;
    int start_row = rank * rows_per_proc + (rank < extra_rows ? rank : extra_rows);
    int num_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);

//    // Crea il tipo per la colonna
//    MPI_Datatype column_type;
//    MPI_Type_vector(N, 1, N, MPI_FLOAT, &column_type);
//    MPI_Type_commit(&column_type);

    // Crea il tipo per il blocco di righe
    MPI_Datatype row_block_type;
    MPI_Type_contiguous(N * num_rows, MPI_FLOAT, &row_block_type);
    MPI_Type_commit(&row_block_type);

    // Ogni processo invia il suo blocco di righe a tutti gli altri
    for (int p = 0; p < size; p++) {
        if (p != rank) {
            MPI_Send(&mat[start_row][0], 1, row_block_type, p, 0, MPI_COMM_WORLD);
        }
    }

    // Il processo corrente copia direttamente il suo blocco
    for (int i = start_row; i < start_row + num_rows; i++) {
        for (int j = 0; j < N; j++) {
            tam[j][i] = mat[i][j];
        }
    }

    // Ricevi i blocchi dagli altri processi
    for (int p = 0; p < size; p++) {
        if (p == rank) continue;

        int other_rows_per_proc = N / size;
        int other_extra = (p < extra_rows ? 1 : 0);
        int other_start = p * rows_per_proc + (p < extra_rows ? p : extra_rows);
        int other_num_rows = other_rows_per_proc + other_extra;

        // Buffer temporaneo per ricevere il blocco
        float *temp_block = (float *)malloc(N * other_num_rows * sizeof(float));

        // Ricevi il blocco intero
        MPI_Recv(temp_block, N * other_num_rows, MPI_FLOAT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Copia nella matrice trasposta
        for (int i = 0; i < other_num_rows; i++) {
            for (int j = 0; j < N; j++) {
                tam[j][other_start + i] = temp_block[i * N + j];
            }
        }

        free(temp_block);
    }

    // Pulisci i tipi di dati
    MPI_Type_free(&column_type);
    MPI_Type_free(&row_block_type);
}


int main(int argc, char** argv) {

    if (argc != 3) {
        printf("Uso: %s <matrix dim> <#processes>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int processes = atoi(argv[2]);

    srand(time(NULL));

    FILE *file = fopen("outputMPI.csv", "a");
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

    MPI_Bcast(mat, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //printf("Rank: %d\n", myrank);
//    for (int i = 0; i < n; i++) {
//      for (int j = 0; j < n; j++) {
//        printf("%f ", mat[i][j]);
//      }
//      printf("\n");
//    }
//    printf("\n");

    start_time_sym = MPI_Wtime();
    checkSymMPI(n, mat, myrank, size, sym);
    end_time_sym = MPI_Wtime();
    total_time_sym = start_time_sym - end_time_sym;

    if (myrank == 0) {
//      printf("sym: %d\n", *sym);
        if (*sym) {
            printf("is sym\n");
        } else {
            printf("is not sym\n");
        }
    }

    start_time = MPI_Wtime();
    //printf("before mat transpose: %d \n", myrank);
    MatTransposeMPI(n, mat, tam, myrank, size);
    //MatTransposeMPI_RowSiz(n, mat, tam, myrank, size);
    //MatTransposeMPI_opt(n, mat, tam, myrank, size);
	end_time = MPI_Wtime();
    total_time_transp = end_time - start_time;


    if (myrank == 0){
   		check(n, mat, tam); //check the correctness of the transposition
    }

    fprintf(file, "%d,%d,%f,%f,MPI\n",processes, n, total_time_sym, total_time_transp);

    MPI_Finalize();

    free(tam);
    free(mat);


    return EXIT_SUCCESS;

}
