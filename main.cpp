#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {

	MPI_Init(NULL, NULL);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	//if (world_rank == 0) { std::cout << "	number of proccess = " << world_size << std::endl; }

	//int RowNum = N / world_size; // N - matrix dimension; RowNum - number of rows per processor 

	MPI_Finalize();
	return 0;
}