#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include <sys/param.h>

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





void MatTransposeBlockMPI(int N, float (*mat)[N], float (*tam)[N], int rank, int size) {
    // Calculate the grid dimensions for process layout
    int grid_dim = (int)sqrt(size);
    if (grid_dim * grid_dim != size) {
        if (rank == 0) {
            printf("Number of processes must be a perfect square\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Calculate block dimensions
    int block_dim = N / grid_dim;
    if (N % grid_dim != 0) {
        if (rank == 0) {
            printf("Matrix dimension must be divisible by sqrt(number of processes)\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Calculate process grid coordinates
    int row = rank / grid_dim;
    int col = rank % grid_dim;

    // Allocate memory for local block
    float *local_block = (float *)malloc(block_dim * block_dim * sizeof(float));

    if (!local_block) {
        printf("Memory allocation failed on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Distribute blocks from rank 0 to all processes
    if (rank == 0) {
        // Process 0 sends blocks to all other processes
        for (int p = 1; p < size; p++) {
            int p_row = p / grid_dim;
            int p_col = p % grid_dim;
            float *temp_block = (float *)malloc(block_dim * block_dim * sizeof(float));

            // Pack block for process p
            for (int i = 0; i < block_dim; i++) {
                for (int j = 0; j < block_dim; j++) {
                    temp_block[i * block_dim + j] =
                        mat[p_row * block_dim + i][p_col * block_dim + j];
                }
            }

            MPI_Send(temp_block, block_dim * block_dim, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
            free(temp_block);
        }

        // Process 0 keeps its own block
        for (int i = 0; i < block_dim; i++) {
            for (int j = 0; j < block_dim; j++) {
                local_block[i * block_dim + j] = mat[i][j];
            }
        }
    } else {
        // Other processes receive their blocks
        MPI_Recv(local_block, block_dim * block_dim, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank == 0) {
        // Process 0 transposes its block directly into tam
        for (int i = 0; i < block_dim; i++) {
            for (int j = 0; j < block_dim; j++) {
                tam[j][i] = local_block[i * block_dim + j];
            }
        }

        // Receive blocks from other processes directly into tam
        for (int p = 1; p < size; p++) {
            int p_row = p / grid_dim;
            int p_col = p % grid_dim;

            // Calculate the destination in the transposed matrix
            int dest_col = p_row * block_dim;  // Original row becomes column
            int dest_row = p_col * block_dim;  // Original column becomes row

            // Receive block row by row directly into tam
            for (int i = 0; i < block_dim; i++) {
                MPI_Recv(&tam[dest_row + i][dest_col], block_dim, MPI_FLOAT,
                        p, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    } else {
        // Other processes transpose their block and send row by row
        float temp;
        // First transpose the local block in-place
        for (int i = 0; i < block_dim; i++) {
            for (int j = i + 1; j < block_dim; j++) {
                temp = local_block[i * block_dim + j];
                local_block[i * block_dim + j] = local_block[j * block_dim + i];
                local_block[j * block_dim + i] = temp;
            }
        }

        // Send each row of the transposed block
        for (int i = 0; i < block_dim; i++) {
            MPI_Send(&local_block[i * block_dim], block_dim, MPI_FLOAT,
                    0, i, MPI_COMM_WORLD);
        }
    }

    // Clean up
    free(local_block);
}


int main(int argc, char** argv) {

    if (argc != 3) {
        printf("Uso: %s <matrix dim> <#processes>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int processes = atoi(argv[2]);

    srand(time(NULL));

    FILE *file = fopen("outputBlock.csv", "a");
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
    double bcast_time_start, bcast_time_end, bcast_time_tot;

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

    total_time_sym = -42;

    start_time = MPI_Wtime();
    //printf("before mat transpose: %d \n", myrank);
    MatTransposeBlockMPI(n, mat, tam, myrank, size);
	end_time = MPI_Wtime();
    total_time_transp = end_time - start_time;

    //MatTransposeBlock(n, mat, tam, 2);

    if (myrank == 0){
   		check(n, mat, tam); //check the correctness of the transposition
        fprintf(file, "%d,%d,%f,%f,MPI BB\n",processes, n, total_time_sym, total_time_transp);
    }


    MPI_Finalize();

    free(tam);
    free(mat);


    return EXIT_SUCCESS;

}
