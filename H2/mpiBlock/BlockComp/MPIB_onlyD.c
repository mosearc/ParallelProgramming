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


void MatTransposeBlockMPI(int N, float (*mat)[N], float (*tam)[N], int rank, int size) { //BASIC DATATYPE
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

    // Create MPI datatype for a block
    MPI_Datatype block_type;
    MPI_Type_vector(block_dim, block_dim, N, MPI_FLOAT, &block_type);
    MPI_Type_commit(&block_type);

    // Create MPI datatype for a contiguous block (used for local storage)
    MPI_Datatype contig_block_type;
    MPI_Type_contiguous(block_dim * block_dim, MPI_FLOAT, &contig_block_type);
    MPI_Type_commit(&contig_block_type);

    // Allocate memory for local blocks
    float *local_block = (float *)malloc(block_dim * block_dim * sizeof(float));
    float *transposed_block = (float *)malloc(block_dim * block_dim * sizeof(float));

    if (!local_block || !transposed_block) {
        printf("Memory allocation failed on process %d\n", rank);
        MPI_Abort(MPI_COMM_WORLD, 1);
        return;
    }

    // Distribute blocks using MPI datatypes
    if (rank == 0) {
        // Process 0 sends blocks to all other processes
        for (int p = 1; p < size; p++) {
            int p_row = p / grid_dim;
            int p_col = p % grid_dim;
            MPI_Send(&mat[p_row * block_dim][p_col * block_dim], 1, block_type,
                    p, 0, MPI_COMM_WORLD);
        }

        // Process 0 keeps its own block
        for (int i = 0; i < block_dim; i++) {
            for (int j = 0; j < block_dim; j++) {
                local_block[i * block_dim + j] = mat[i][j];
            }
        }
    } else {
        // Other processes receive their blocks using contiguous datatype
        MPI_Recv(local_block, 1, contig_block_type, 0, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Transpose local block using SIMD instructions if available
    //#pragma omp parallel for collapse(2)
    for (int i = 0; i < block_dim; i++) {
        for (int j = 0; j < block_dim; j++) {
            transposed_block[j * block_dim + i] = local_block[i * block_dim + j];
        }
    }

    // Create MPI datatype for transposed block placement
    MPI_Datatype transposed_block_type;
    MPI_Type_vector(block_dim, block_dim, N, MPI_FLOAT, &transposed_block_type);
    MPI_Type_commit(&transposed_block_type);

    // Gather transposed blocks to process 0
    if (rank == 0) {
        // Process 0 places its block
        for (int i = 0; i < block_dim; i++) {
            for (int j = 0; j < block_dim; j++) {
                tam[i][j] = transposed_block[i * block_dim + j];
            }
        }

        // Receive blocks from other processes
        for (int p = 1; p < size; p++) {
            int p_row = p / grid_dim;
            int p_col = p % grid_dim;
            int dest_row = p_col;  // Swapped for transposition
            int dest_col = p_row;

            float *recv_block = (float *)malloc(block_dim * block_dim * sizeof(float));
            MPI_Recv(recv_block, 1, contig_block_type, p, 1,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Place block using MPI datatype
            for (int i = 0; i < block_dim; i++) {
                for (int j = 0; j < block_dim; j++) {
                    tam[dest_row * block_dim + i][dest_col * block_dim + j] =
                        recv_block[i * block_dim + j];
                }
            }
            free(recv_block);
        }
    } else {
        // Other processes send their transposed blocks
        MPI_Send(transposed_block, 1, contig_block_type, 0, 1, MPI_COMM_WORLD);
    }

    // Clean up
    MPI_Type_free(&block_type);
    MPI_Type_free(&contig_block_type);
    MPI_Type_free(&transposed_block_type);
    free(local_block);
    free(transposed_block);
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
        // for (int i = 0; i < n; i++) {
        //     for (int j = 0; j < n; j++) {
        //         printf("%f ", mat[i][j]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");
    }

    //MPI_Bcast(mat, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //printf("Rank: %d\n", myrank);
    // for (int i = 0; i < n; i++) {
    //   for (int j = 0; j < n; j++) {
    //     printf("%f ", mat[i][j]);
    //   }
    //   printf("\n");
    // }
    // printf("\n");

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
    //MatTransposeMPI(n, mat, tam, myrank, size);
    //MatTransposeMPI_RowSiz(n, mat, tam, myrank, size);
    //MatTransposeMPI_opt(n, mat, tam, myrank, size);
    MatTransposeBlockMPI(n, mat, tam, myrank, size);
	end_time = MPI_Wtime();
    total_time_transp = end_time - start_time;

    //MatTransposeBlock(n, mat, tam, 2);

    if (myrank == 0){
        // for (int i = 0; i < n; i++) {
        //     for (int j = 0; j < n; j++) {
        //         printf("%f ", tam[i][j]);
        //     }
        //     printf("\n");
        // }
        //printf("\n");
   		check(n, mat, tam); //check the correctness of the transposition
    }

    fprintf(file, "%d,%d,%f,%f,MPI BOD\n",processes, n, total_time_sym, total_time_transp);

    MPI_Finalize();

    free(tam);
    free(mat);


    return EXIT_SUCCESS;

}
