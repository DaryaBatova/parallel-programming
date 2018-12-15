#include "mpi.h"
#include <iostream>
#include <time.h>

//enum Tags {size_merge, send_merge} tags;
void ShellSort(int* array, int n)
{
	double tmp;
	for (int step = n / 2; step > 0; step /= 2)
		for (int i = step; i < n; i++)
			for (int j = i - step; (j >= 0) && (array[j] > array[j + step]); j -= step)
			{
				tmp = array[j];
				array[j] = array[j + step];
				array[j + step] = tmp;
			}
}



int* FillingArray(int n)
{
	srand(time(nullptr));
	int* array = new int[n];
	for (int i = 0; i < n; i++)
		array[i] = rand() % 100;
	return array;
}

void PrintArray(int* array, int n)
{
	for (int i = 0; i < n; i++)
		std::cout << array[i] << " ";
	std::cout << std::endl;
}

void Sendcounts(int* sendcounts, int ProcNum, int size, int rem)
{
	sendcounts = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
	{
		int k = rem - i;
		sendcounts[i] = k > 0 ? size + 1 : size;
	}
}

void Displs(int* displs, int ProcNum, int size, int rem)
{
	displs = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
	{
		int countsSize = i;
		int countsRem = 0;
		if (rem != 0)
		{
			countsRem = i < rem ? i : rem;
			countsSize = i < rem ? 0 : i - rem;
		}
		displs[i] = countsRem*(size + 1) + countsSize*size;
	}
}

bool Comparison(int& value1, int& value2)
{
	if (value1 > value2)
	{
		std::swap(value1, value2);
		return true;
	}
	else
	{
		return false;
	}
}

template <typename T>
void Comparator(T* arr, size_t size)
{
	for (size_t i = 1; i < size; i++)
	{
		if (arr[i] < arr[i - 1])
			std::swap(arr[i], arr[i - 1]);
	}
}

template <typename T>
void EvenSplitter(T* arr, size_t size1, size_t size2)
{
	T* arr1 = arr;
	T* arr2 = arr + size1;
	std::cout << std::endl;
//	ShowArray(tmp, size1, 250);
	std::cout << std::endl;
	ShowArray(arr, size1 + size2, 250);
	std::cout << std::endl;
	int a = 0;
	int b = 0;
	int i = 0;


	while ((a < size1) && (b < size2))
	{
		if (arr1[a] >= arr2[b])
		{
			std::swap(arr1[i], arr2[b]);
			a += 2;
		}
		else
		{
			b += 2;
		}
		i += 2;
	}

	ShowArray(arr, size1 + size2, 250);
	std::cout << std::endl;


	if (a == size1) {
		for (int j = b; j < size2; j += 2, i += 2)
			std::swap(arr1[i], arr2[j]);
	}
	else {
		// части уже стоят в нужном порядке
	}
		
	ShowArray(arr, size1 + size2, 250);

}

template <typename T>
void OddSplitter(T* array1, int size1, T* tmp, int size2)
{
	for (int i = 1; i < size1; i+=2)
		tmp[i] = array1[i];
	int* array2 = array1 + size1;
	int a = 1;
	int b = 1;
	int i = 1;

	while ((a < size1) && (b < size2))
	{
		if (tmp[a] <= array2[b])
		{
			array1[i] = tmp[a];
			a += 2;
		}
		else
		{
			array1[i] = array2[b];
			b += 2;
		}
		i += 2;
	}

	if (a == size1)
		for (int j = b; j < size2; j += 2, i += 2)
			array1[i] = array2[j];
	else
		for (int j = a; j < size1; j += 2, i += 2)
			array1[i] = tmp[j];


}

template<typename T>
void ShowArray(T * pdArray, size_t unSizeArr, int ProcRank)
{
	for (size_t i = 0; i < unSizeArr; ++i)
	{
		std::cout << "ProcRank " << ProcRank << " has a folowing arr[" << i << "] = " << pdArray[i] << std::endl;
	}
}


template <typename T>
void BatcherMerge(T* arr1, T* arr2, T* res, size_t size1, size_t size2)
{
	//ShowArray(arr1, size1, 150);
	//ShowArray(arr2, size1, 150);
	//memcpy(res, arr1, size1);
	for (size_t i = 0; i < size1; i++)
	{
		res[i] = arr1[i];
	}
	for (size_t j = 0; j < size2; j++)
	{
		res[size1 + j] = arr2[j];
	}
	//ShowArray(res, size1, 150);
	//memcpy(res + size1, arr2, size2);
	ShowArray(res, size1 + size2, 150);
	std::cout << std::endl;
	EvenSplitter(res, size1, size2);
	//ShowArray(res, size1 + size2, 150);
	//std::cout << std::endl;
	OddSplitter(res, size1, arr2, size2);
	//ShowArray(res, size1 + size2, 150);
	std::cout << std::endl;
	Comparator(res, size1 + size2);
	//ShowArray(res, size1 + size2, 150);
	std::cout << std::endl;
}

template <typename T>
void reallocate(T* arr, size_t oldsize, size_t newsize)
{
	T* newarr = new T[newsize];
	memcpy(newarr, arr, newsize);
	delete[] arr;
	arr = newarr;
}

double LinearVersionShellSort(int *array, int size)
{
	double StartTime, FinishTime;
	int *newarray = new int[size];
	for (int i = 0; i < size; i++)
		newarray[i] = array[i];
	StartTime = MPI_Wtime();
	ShellSort(newarray, size);
	FinishTime = MPI_Wtime();
	std::cout << "ShellSort on one process: ";
	PrintArray(newarray, size);
	delete[] newarray;
	return FinishTime - StartTime;
}

template <typename T>
int NEW_TREE_Reduce(T *sendbuf, T *recvbuf, int* counts, int sizearr, MPI_Datatype type, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int RankRecv, RankSend;
	std::cout << "ResBuf's size " << sizearr << std::endl;
	T * ResBuf = new T[sizearr];
	int RecvArrSize = 0;
	int SendArrSize = counts[ProcRank];
	int newProcRank = (ProcRank - root + ProcNum) % ProcNum;
	int mask = 1;
	while (mask < ProcNum)
	{
		if ((newProcRank & mask) == 0)
		{
			RankSend = newProcRank | mask;
			if (RankSend < ProcNum)
			{
				RankSend = (RankSend + root) % ProcNum;
				MPI_Recv(&RecvArrSize, 1, MPI_INT, RankSend, 0, comm, &status);
				MPI_Recv(recvbuf, RecvArrSize, type, RankSend, 0, comm, &status);
				std::cout << ProcRank << " recv from " << RankSend << " array with size " << RecvArrSize << std::endl;
				ShowArray(recvbuf, RecvArrSize, ProcRank);
				BatcherMerge(sendbuf, recvbuf, ResBuf, counts[ProcRank], counts[RankSend]);
				SendArrSize = counts[ProcRank];
				reallocate(sendbuf, SendArrSize, SendArrSize + RecvArrSize);
				SendArrSize += RecvArrSize;
				memcpy(sendbuf, ResBuf, SendArrSize);
			}
		}
		else
		{
			RankRecv = newProcRank&(~mask);
			RankRecv = (RankRecv + root) % ProcNum;
			std::cout << ProcRank << " send to " << RankRecv << " array with size " << SendArrSize << std::endl;
			ShowArray(sendbuf, SendArrSize, ProcRank);
			MPI_Send(&SendArrSize, 1, MPI_INT, RankRecv, 0, comm);
			MPI_Send(sendbuf, SendArrSize, type, RankRecv, 0, comm);
			break;
		}
		mask = mask << 1;
	}
	if (ProcRank != root)
	{
		delete[] sendbuf;
		delete[] ResBuf;
	}

	if (ProcRank == root)
	{
		memcpy(recvbuf, ResBuf, sizearr);
		delete[] sendbuf, ResBuf;
	}

	return 0;
}


int main(int argc, char *argv[])
{
	int ProcNum, ProcRank;
	int *array = nullptr; 
	//int *displs = nullptr;
	//int *sendcounts = nullptr;
	double stime, ftime;
	int size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (argc == 2)
		size = atoi(argv[1]);
	else size = 10;

	if (ProcRank == 0)
	{
		array = FillingArray(size);
		PrintArray(array, size);
		//double linTime = LinearVersionShellSort(array, size);
		//std::cout << "Time for linear version: " << linTime << std::endl;
		//Displs(displs, ProcNum, size / ProcNum, size % ProcNum);
	}

	//Sendcounts(sendcounts, ProcNum, size / ProcNum, size % ProcNum);
	int * sendcounts = new int[ProcNum];
	int * displs = new int[ProcNum];
	if (ProcRank == 0)
	{
		int nSendSum = 0;
		int nRecvSum = 0;
		int nRemSizePerProc = size % ProcNum;
		for (size_t i = 0; i < ProcNum; i++)
		{
			sendcounts[i] = size / ProcNum;
			if (nRemSizePerProc > 0)
			{
				sendcounts[i]++;
				nRemSizePerProc--;
			}
			displs[i] = nRecvSum;
			nRecvSum += sendcounts[i];
		}
	}
	MPI_Bcast(sendcounts, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(displs, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);

	
	int *recvbuf = new int[sendcounts[ProcRank]];
	if (ProcRank == 0)
	{
		stime = MPI_Wtime();
	}
	MPI_Scatterv(array, sendcounts, displs, MPI_INT, recvbuf, sendcounts[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
	ShellSort(recvbuf, sendcounts[ProcRank]);
	int * res = new int[size];
	//ShowArray(recvbuf, size, ProcRank);
	NEW_TREE_Reduce(recvbuf, res, sendcounts, size, MPI_INT, 0, MPI_COMM_WORLD);

	
	if (ProcRank == 0)
	{
		ftime = MPI_Wtime();
		PrintArray(res, size);
		std::cout << "Parallel version: " << ftime - stime << std::endl;

	}
	MPI_Finalize();
	return 0;
}
