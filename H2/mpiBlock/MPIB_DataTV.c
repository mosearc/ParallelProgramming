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


 void MatTransposeBlockMPI(int N, float (*mat)[N], float (*tam)[N], int rank, int size) { //DATATYPE AND SIV

     int grid_dim = (int)sqrt(size);
     if (grid_dim * grid_dim != size) {
         if (rank == 0) {
             printf("Number of processes must be a perfect square\n");
         }
         MPI_Abort(MPI_COMM_WORLD, 1);
         return;
     }

     int block_dim = N / grid_dim;
     if (N % grid_dim != 0) {
         if (rank == 0) {
             printf("Matrix dimension must be divisible by sqrt(number of processes)\n");
         }
         MPI_Abort(MPI_COMM_WORLD, 1);
         return;
     }

     int row = rank / grid_dim;
     int col = rank % grid_dim;

     MPI_Datatype source_block_type, source_resized_type;
     MPI_Type_vector(block_dim, block_dim, N, MPI_FLOAT, &source_block_type);
     MPI_Type_create_resized(source_block_type, 0, sizeof(float), &source_resized_type);
     MPI_Type_commit(&source_resized_type);


     MPI_Datatype dest_block_type, dest_resized_type;
     MPI_Type_vector(block_dim, block_dim, N, MPI_FLOAT, &dest_block_type);
     MPI_Type_create_resized(dest_block_type, 0, sizeof(float), &dest_resized_type);
     MPI_Type_commit(&dest_resized_type);


     float *local_block = (float *)malloc(block_dim * block_dim * sizeof(float));
     float *transposed_block = (float *)malloc(block_dim * block_dim * sizeof(float));
     if (!local_block || !transposed_block) {
         printf("Memory allocation failed on process %d\n", rank);
         MPI_Abort(MPI_COMM_WORLD, 1);
         return;
     }

     int *sendcounts = NULL;
     int *displs = NULL;
     if (rank == 0) {
         sendcounts = (int *)malloc(size * sizeof(int));
         displs = (int *)malloc(size * sizeof(int));
         for (int i = 0; i < size; i++) {
             sendcounts[i] = 1;
             int row_i = i / grid_dim;
             int col_i = i % grid_dim;
             displs[i] = row_i * N * block_dim + col_i * block_dim;
         }
     }


     MPI_Scatterv(&mat[0][0], sendcounts, displs, source_resized_type,
                  local_block, block_dim * block_dim, MPI_FLOAT,
                  0, MPI_COMM_WORLD);

     for (int i = 0; i < block_dim; i++) {
         for (int j = 0; j < block_dim; j++) {
             transposed_block[j * block_dim + i] = local_block[i * block_dim + j];
         }
     }

     int *recvcounts = NULL;
     int *rdispls = NULL;
     if (rank == 0) {
         recvcounts = (int *)malloc(size * sizeof(int));
         rdispls = (int *)malloc(size * sizeof(int));
         for (int i = 0; i < size; i++) {
             recvcounts[i] = 1;
             int row_i = i / grid_dim;
             int col_i = i % grid_dim;

             rdispls[i] = col_i * N * block_dim + row_i * block_dim;
         }
     }


     MPI_Gatherv(transposed_block, block_dim * block_dim, MPI_FLOAT,
                 &tam[0][0], recvcounts, rdispls, dest_resized_type,
                 0, MPI_COMM_WORLD);


     MPI_Type_free(&source_block_type);
     MPI_Type_free(&source_resized_type);
     MPI_Type_free(&dest_block_type);
     MPI_Type_free(&dest_resized_type);
     free(local_block);
     free(transposed_block);
     if (rank == 0) {
         free(sendcounts);
         free(displs);
         free(recvcounts);
         free(rdispls);
     }
 }



int main(int argc, char** argv) {

    if (argc != 3) {
        printf("Uso: %s <matrix dim> <#processes>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    int processes = atoi(argv[2]);

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




     total_time_sym = -42;



    start_time = MPI_Wtime();

    MatTransposeBlockMPI(n, mat, tam, myrank, size);

	end_time = MPI_Wtime();
    total_time_transp = end_time - start_time;




    if (myrank == 0){
   		check(n, mat, tam); //check the correctness of the transposition
        fprintf(file, "%d,%d,%f,%f,MPI BDV\n",processes, n, total_time_sym, total_time_transp);

    }


    MPI_Finalize();

    free(tam);
    free(mat);


    return EXIT_SUCCESS;

}
