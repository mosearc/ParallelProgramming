#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include <sys/time.h>


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

		}
	}
}

int checkSymMPI(int n, float(* mat)[n], int myrank, int size) {
	int ret_val = 1;
	int tot = 1;

	int rows_per_proc = n / size;
	int start_row = myrank * rows_per_proc;
	int end_row = start_row + rows_per_proc;

	if (myrank == size - 1) {
		end_row = n;
	}

	for (int i = start_row; i < end_row; i++) {
		for (int j = 0; j < n; j++) {
			if (mat[i][j] != mat[j][i]) { //isn't sym
				ret_val = 0;
			}
		}
	}
	MPI_Reduce(&ret_val, &tot, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);

	return tot;
}



//only for number of processes divisor of dim
void MatTransposeMPI_limited(int n, float(* mat)[n], float (*tam)[n], int myrank, int size) {

	int rows_per_proc = n / size;

	float* local_matrix = (float*)malloc(rows_per_proc * n * sizeof(float));
	float* local_transpose = (float*)malloc(n * rows_per_proc * sizeof(float));

	if (myrank == 0) {
		MPI_Scatter(mat, rows_per_proc * n, MPI_FLOAT, local_matrix, rows_per_proc * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		MPI_Scatter(NULL, 0, MPI_DOUBLE, local_matrix, rows_per_proc * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	MPI_Datatype block_type, resized_type;
	MPI_Type_vector(rows_per_proc, 1, n, MPI_DOUBLE, &block_type);
	MPI_Type_create_resized(block_type, 0, sizeof(double), &resized_type);
	MPI_Type_commit(&resized_type);

	MPI_Alltoall(local_matrix, rows_per_proc, MPI_DOUBLE, local_transpose, 1, resized_type, MPI_COMM_WORLD);

	if (myrank == 0) {
		MPI_Gather(local_transpose, n * rows_per_proc, MPI_DOUBLE, tam, n * rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	} else {
		MPI_Gather(local_transpose, n * rows_per_proc, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	MPI_Type_free(&resized_type);
	free(local_matrix);
	free(local_transpose);
}



void MatTransposeMPI(int N, float(* mat)[N], float (*tam)[N], int rank, int size) {

	int base_rows_per_proc = N / size;
	int remainder_rows = N % size;

	int rows_per_proc = base_rows_per_proc + (rank < remainder_rows ? 1 : 0);

	int* sendcounts = NULL;
	int* displs = NULL;
	if (rank == 0) {
		sendcounts = (int*)malloc(size * sizeof(int));
		displs = (int*)malloc(size * sizeof(int));

		int offset = 0;
		for (int i = 0; i < size; i++) {
			sendcounts[i] = (base_rows_per_proc + (i < remainder_rows ? 1 : 0)) * N;
			displs[i] = offset * N;
			offset += base_rows_per_proc + (i < remainder_rows ? 1 : 0);
		}
	}

	float* local_matrix = (float*)malloc(rows_per_proc * N * sizeof(float));
	if (!local_matrix) {
		printf("Memory allocation failed for local_matrix on process %d\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (rank == 0) {
		MPI_Scatterv(mat, rows_per_proc * N, MPI_FLOAT, local_matrix, rows_per_proc * N, MPI_FLOAT, 0, MPI_COMM_WORLD);
	} else {
		MPI_Scatterv(NULL, NULL, NULL, local_matrix, rows_per_proc * N, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}


	int* recv_counts = (int*)malloc(size * sizeof(int));
	int* send_counts = (int*)malloc(size * sizeof(int));
	int* rdispls = (int*)malloc(size * sizeof(int));
	int* sdispls = (int*)malloc(size * sizeof(int));

	MPI_Allgather(&rows_per_proc, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < size; i++) {
		send_counts[i] = rows_per_proc * recv_counts[i];
	}

	sdispls[0] = 0;
	rdispls[0] = 0;
	for (int i = 1; i < size; i++) {
		sdispls[i] = sdispls[i-1] + send_counts[i-1];
		rdispls[i] = rdispls[i-1] + recv_counts[i-1];
	}

	int total_recv_elements = 0;
	for (int i = 0; i < size; i++) {
		total_recv_elements += recv_counts[i];
	}

	float* local_transpose = (float*)malloc(total_recv_elements * rows_per_proc * sizeof(float));
	if (!local_transpose) {
		printf("Memory allocation failed for local_transpose on process %d\n", rank);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Perform all-to-all communication with variable sizes
	MPI_Alltoallv(local_matrix, send_counts, sdispls, MPI_FLOAT, local_transpose, recv_counts, rdispls, MPI_FLOAT, MPI_COMM_WORLD);

	if (rank == 0) {
		MPI_Gatherv(local_transpose, total_recv_elements * rows_per_proc, MPI_FLOAT, tam, sendcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
		free(sendcounts);
		free(displs);
	} else {
		MPI_Gatherv(local_transpose, total_recv_elements * rows_per_proc, MPI_FLOAT, NULL, NULL, NULL, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	free(local_matrix);
	free(local_transpose);
	free(recv_counts);
	free(send_counts);
	free(rdispls);
	free(sdispls);

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

	init_mat(n, mat);


	int myrank, size;
	double start_time, end_time, start_time_sym, end_time_sym, total_time_transp, total_time_sym;
	int sym;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	start_time_sym = MPI_Wtime();
	sym = checkSymMPI(n, mat, myrank, size);
	end_time_sym = MPI_Wtime();

	start_time = MPI_Wtime();
	MatTransposeMPI(n, mat, tam, myrank, size);
	end_time = MPI_Wtime();

	MPI_Finalize();

	total_time_sym = end_time_sym - start_time_sym;
	total_time_transp = end_time - start_time;

	check(n, mat, tam); //check the correctness of the transposition

	fprintf(file, "%d,%d,%f,%f,MPI\n",size, n, total_time_sym, total_time_transp);

	free(tam);
	free(mat);

    return EXIT_SUCCESS;
}
