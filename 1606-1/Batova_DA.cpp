#include <iostream>
#include "mpi.h"
#include <ctime>

void FillingArray(double array[], int n)
{
	for (int i = 0; i < n; i++)
	{
		srand(time(nullptr));
		array[i] = rand() % 10;
		std::cout << array[i] << " ";
	}
}
int main(int argc, char **argv)
{
	double *array;
	int k = 100;
	if (argc == 2)
	{
		k = atoi(argv[1]);
	}
	int ProcNum, ProcRank;
	double ProcSum, TotalSum = 0.0;
	array = new double[k];
	double *recbuf = new double[k];
	int sum = 0;
	double tstart, tfinish;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int rem = k % ProcNum;
	if (ProcRank == 0)
	{
		FillingArray(array, k);
		tstart = MPI_Wtime();
	}
	if (rem != k)
	{
		int *sendcounts = new int[ProcNum * sizeof(int)];
		int *displs = new int[ProcNum * sizeof(int)];
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
		MPI_Scatterv(array, sendcounts, displs, MPI_DOUBLE, recbuf, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int i = 0; i < sendcounts[ProcRank]; i++)
		{
			ProcSum += recbuf[i];
		}
		delete[] recbuf, sendcounts, displs;
	}
	else if (rem == k)
	{
		MPI_Bcast(array, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int i = ProcRank; i < ProcRank + 1; i++)
		{
			if (i < k)
			{
				if (rem > 0)
				{
					ProcSum += array[i];
					rem--;
				}
			}
			else ProcSum = 0.0;
		}
	}
	MPI_Reduce(&ProcSum, &TotalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		std::cout << "TotalSum = " << TotalSum << std::endl;
		tfinish = MPI_Wtime();
		std::cout << "Wtime = " << tfinish - tstart << std::endl;
	}
	//для одного процесса
	double tstart1;
	double TotalSum1 = 0.0;
	if (ProcRank == 0)
	{
		tstart1 = MPI_Wtime();
		for (int i = 0; i < k; i++)
			TotalSum1 += array[i];
		std::cout << "TotalSum = " << TotalSum1 << std::endl;
		std::cout << "Wtime for one process = " << MPI_Wtime() - tstart1 << std::endl;
	}

	MPI_Finalize();
	delete[] array;
	return 0;
}