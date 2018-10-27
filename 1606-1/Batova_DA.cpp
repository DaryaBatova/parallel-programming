#include <iostream>
#include "mpi.h"
#include <ctime>
#include <limits>

void FillingArray(double array[], int n, int min, int max) //function to fill the array with real numbers in the range [min, max]
{
	srand(time(nullptr));
	for (int i = 0; i < n; i++)
	{
		double ri = (double)rand() / RAND_MAX;
		array[i] = min + (max - min)*ri;
		std::cout << array[i] << " ";
	}
}
int main(int argc, char **argv)
{
	int k = 10; 
	int ProcNum, ProcRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	if (ProcRank == 0) //zero process sets the size of the array (entered by the user via the command line)
	{
		if (argc == 4)
		{
			k = atoi(argv[1]);
		}
		else k = 10;
	}

	MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD); //broadcasts a size of the array from the process with rank "0" to all other processes of the communicator
	double tstart, tfinish;
	double *array;
	double TotalSum = 0.0; //total amount. counted on zero process

	if (ProcRank == 0)
	{
		array = new double[k];
		FillingArray(array, k, atoi(argv[2]), atoi(argv[3]));
		tstart = MPI_Wtime();
	}
	else {
		array = nullptr;
	}

	int rem = k % ProcNum; 

	int *sendcounts = new int[ProcNum]; //integer array (of length group size) specifying the number of elements to send to each processor
	int *displs = new int[ProcNum]; //integer array (of length group size). Entry i specifies the displacement (relative to sendbuf from which to take the outgoing data to process i
	int sum = 0; //sum of counts. used to calculate displacements
	for (int i = 0; i < ProcNum; i++)
	{
		sendcounts[i] = k / ProcNum;
		if (rem > 0)
		{
			sendcounts[i]++;
			rem--;
		}
		displs[i] = sum;
		sum += sendcounts[i];
	}
	
	double *recbuf = new double[sendcounts[ProcRank]]; //address of receive buffer

	double ProcSum = 0.0; //partial amount for each process

	MPI_Scatterv(array, sendcounts, displs, MPI_DOUBLE, recbuf, sendcounts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD); //sendcounts[ProcRank] - number of elements in receive buffer (recbuf)

	for (int i = 0; i < sendcounts[ProcRank]; i++) //considers partial amount
	{
		ProcSum += recbuf[i];
	}

	MPI_Reduce(&ProcSum, &TotalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //reduces partial amount on all processes to a total amount
	tfinish = MPI_Wtime();
	if (ProcRank == 0)
	{
		std::cout << std::endl;
		std::cout << "Parallel version (using functions Scatterv and Reduce): " << std::endl;
		std::cout << "Elapsed time = " << tfinish - tstart << std::endl;
		std::cout << "TotalSum = " << TotalSum << std::endl;
		
	}
	//Linear version of the algorithm
	double tstart1;
	double TotalSum1 = 0.0;
	if (ProcRank == 0)
	{
		tstart1 = MPI_Wtime();
		for (int i = 0; i < k; i++)
			TotalSum1 += array[i];
		std::cout << std::endl;
		std::cout << "Linear version (one process): " << std::endl;
		std::cout << "TotalSum = " << TotalSum1 << std::endl;
		std::cout << "Elapsed time = " << MPI_Wtime() - tstart1 << std::endl;
	}

	MPI_Finalize();
	delete[] recbuf, sendcounts, displs, array;
	return 0;
}