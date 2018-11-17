#include "mpi.h"
#include <iostream>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <cassert>
#include<cstdlib>


int NEW_TREE_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm);
size_t stype(MPI_Datatype type);
int *createarray(int n);
void arrayinit(int *array, int n, int min, int max);



int main(int argc, char **argv)
{
	int ProcNum, ProcRank;
	int *array;
	int n, root, k;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	
	n = atoi(argv[1]);
	array = new int [n];
	arrayinit(array, n, atoi(argv[2]), atoi(argv[3]));
	root = atoi(argv[4]);
	int * res;
	if (ProcRank == root)
	{
		res = new int[n];
		for (int i = 0; i < n; i++)
		{
			std::cout << array[i] << " ";
		}
		std::cout << std::endl;
	}
	k = atoi(argv[5]);
	double tstart1;
	
	if (ProcRank == root)
	{
		tstart1 = MPI_Wtime();
	}
	if (k == 0)
		MPI_Reduce(array, res, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
	else if (k == 1)
		NEW_TREE_Reduce(array, res, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
	if (ProcRank == root)
	{
		std::cout << std::endl;
		if (k == 0)
			std::cout << "MPI_Reduce: " << std::endl;
		else if (k == 1)
			std::cout << "TREE_Reduce: " << std::endl;
		std::cout << "Elapsed time = " << MPI_Wtime() - tstart1 << std::endl;
		std::cout << "You are in process " << root << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << res[i] << " ";
	}
	
	
	MPI_Finalize();
	return 0;
	
}

int *createarray(int n)
{
	int *array = new int[n];

	return array;
}
void arrayinit(int *array, int n, int min, int max)
{

	for (int i = 0; i < n; i++)
	{
		array[i] = min + rand() % (max - min);
	}
}


size_t stype(MPI_Datatype type)
{
	size_t result;
	switch (type)
	{
	default:
		std::cout << "Incorrect datatype" << std::endl;
		exit(1);
	case MPI_INT:
		result = sizeof(int);
		break;
	case MPI_FLOAT:
		result = sizeof(float);
		break;
	case MPI_DOUBLE:
		result = sizeof(double);
		break;
	}
	return result;
}

void max (void * buf1, void * buf2, void * res, size_t size, MPI_Datatype type)
{
	switch (type)
	{
	case MPI_INT:
	{
		int * buffer = static_cast<int*>(res);
		int * buffer1 = static_cast<int*>(buf1);
		int * buffer2 = static_cast<int*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] > buffer2[i];
		}
		break;
	}
	case MPI_FLOAT:
	{
		float * buffer = static_cast<float*>(res);
		float * buffer1 = static_cast<float*>(buf1);
		float * buffer2 = static_cast<float*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] > buffer2[i];
		}
		break;
	}
	case MPI_DOUBLE:
	{
		double * buffer = static_cast<double*>(res);
		double * buffer1 = static_cast<double*>(buf1);
		double * buffer2 = static_cast<double*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] > buffer2[i];
		}
		break;
	}
	}
}

void min(void * buf1, void * buf2, void * res, size_t size, MPI_Datatype type)
{
	switch (type)
	{
	case MPI_INT:
	{
		int * buffer = static_cast<int*>(res);
		int * buffer1 = static_cast<int*>(buf1);
		int * buffer2 = static_cast<int*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] < buffer2[i];
		}
		break;
	}
	case MPI_FLOAT:
	{
		float * buffer = static_cast<float*>(res);;
		float * buffer1 = static_cast<float*>(buf1);
		float * buffer2 = static_cast<float*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] < buffer2[i];
		}
		break;
	}
	case MPI_DOUBLE:
	{
		double * buffer = static_cast<double*>(res);;
		double * buffer1 = static_cast<double*>(buf1);
		double * buffer2 = static_cast<double*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] < buffer2[i];
		}
		break;
	}
	}
}

void sum(void * buf1, void * buf2, void* res, size_t size, MPI_Datatype type)
{
	switch (type)
	{
	case MPI_INT:
	{
		int * buffer = static_cast<int*>(res);
		int * buffer1 = static_cast<int*>(buf1);
		int * buffer2 = static_cast<int*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] + buffer2[i];
		}
		break;
	}
	case MPI_FLOAT:
	{
		float * buffer = static_cast<float*>(res);
		float * buffer1 = static_cast<float*>(buf1);
		float * buffer2 = static_cast<float*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] + buffer2[i];
		}
		break;	
	}
	case MPI_DOUBLE:
	{
		double * buffer = static_cast<double*>(res);
		double * buffer1 = static_cast<double*>(buf1);
		double * buffer2 = static_cast<double*>(buf2);
		for (size_t i = 0; i < size; i++)
		{
			buffer[i] = buffer1[i] + buffer2[i];
		}
		break;	
	}
	}
}

void Operation(MPI_Datatype type, MPI_Op op, void* buf1, void* buf2, void* res, size_t size)
{
	switch (op)
	{
	case MPI_MAX:
		max(buf1, buf2, res, size, type);
		break;
	case MPI_MIN:
		min(buf1, buf2, res, size, type);
		break;
	case MPI_SUM:
		sum(buf1, buf2, res, size, type);
		break;
	default:
		break;
	}
}

int number_in_old_degree(int num)
{
	int tmp = 1 << 30;
	while (num < tmp)
	{
		tmp >>= 1;
	}
	return tmp;
}

int old_degree(int num)
{
	int tmp = 1 << 30;
	int count = 31;
	while (num < tmp)
	{
		tmp >>= 1;
		count--;
	}
	return count;
}




int NEW_TREE_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int RankSend{ 0 };
	void * res;
	int finish = std::log(ProcNum) / std::log(2);
	if (ProcNum % 2 > 0)
		finish++;

	int newProcRank = (ProcRank < root ? ProcRank + 1 : ProcRank);
	int start = (ProcRank == root ? 0 : old_degree(newProcRank));
	for (int i = finish; i >= start; i--)
	{
		if (ProcRank == root)
		{
			RankSend = std::pow(2, i);
			RankSend <= root ? RankSend-- : RankSend;
		}
		else
		{
			RankSend = std::pow(2, i) + newProcRank;
			if (RankSend <= root)
				RankSend--;
		}
		if (RankSend < ProcNum)
		{
			MPI_Recv(recvbuf, count, type, RankSend, 0, comm, &status);
			Operation(type, op, sendbuf, recvbuf, sendbuf, count);
		}
	}
	if (ProcRank != root)
	{
		int oldnum = number_in_old_degree(newProcRank);
		int RankRecv = newProcRank & (~oldnum);
		if (RankRecv == 0)
			RankRecv = root;
		else if (RankRecv <= root)
			RankRecv--;
		MPI_Send(sendbuf, count, type, RankRecv, 0, comm);
	}

	if (ProcRank == root)
	{
		memcpy(recvbuf, sendbuf, count*stype(type));
	}
	
	return 0;
}