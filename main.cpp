#include <iostream>
#include <mpi.h>
#include <map>
#include <vector>
#include "Header.h"


#include <string>

int main() {

	MPI_Init(NULL, NULL);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if (world_rank == 0) { std::cout << "	number of proccess = " << world_size << std::endl; }

	solver s;
	s.initialization(world_size);
	if (world_rank == 0) {
		//s.printA();
	}
	
	s.Jacobi(world_rank);
	s.printAnswer(world_rank);
	    
	MPI_Finalize();
	return 0;
}