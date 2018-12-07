#include "mpi.h"
#include <iostream>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <cassert>
#include<cstdlib>


int NEW_TREE_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm);
size_t stype(MPI_Datatype type);
void arrayinit(int *array, int n, int min, int max);
int* NodePosition(int ProcNum, int root, int ProcRank);
int ReduceLeftRight(void* sendbuf, void* recvbuf, void* res, int count, MPI_Datatype type, int left, int right, int root, MPI_Op op, MPI_Comm comm);
int NEW_TREE_Reduce_Binary_Tree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm);

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
	switch (k)
	{
	case 0:
		MPI_Reduce(array, res, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
		break;
	case 1:
		NEW_TREE_Reduce(array, res, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
		break;
	case 2:
		NEW_TREE_Reduce_Binary_Tree(array, res, n, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
		break;
	default:
		MPI_Finalize();
		return 0;
	}

	if (ProcRank == root)
	{
		std::cout << std::endl;
		if (k == 0)
			std::cout << "MPI_Reduce: " << std::endl;
		else if (k == 1)
			std::cout << "TREE_Reduce_Binomial_Tree: " << std::endl;
		else if(k == 2)
			std::cout << "TREE_Reduce_Binary_Tree: " << std::endl;
		std::cout << "Elapsed time = " << MPI_Wtime() - tstart1 << std::endl;
		std::cout << "You are in process " << root << std::endl;
		for (int i = 0; i < n; i++)
			std::cout << res[i] << " ";
	}
	
	
	MPI_Finalize();
	return 0;
	
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


int NEW_TREE_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int RankRecv, RankSend;
	void * res = malloc(count * stype(type));
	if (ProcRank != root)
		recvbuf = malloc(count * stype(type));
	memcpy(recvbuf, sendbuf, count*stype(type));
	int newProcRank = (ProcRank - root + ProcNum) % ProcNum;
	int mask = 1;
	while(mask < ProcNum)
	{
		if ((newProcRank & mask) == 0)
		{
			RankSend = newProcRank|mask;
			if (RankSend < ProcNum)
			{
				RankSend = (RankSend + root) % ProcNum;
				MPI_Recv(res, count, type, RankSend, 0, comm, &status);
				Operation(type, op, recvbuf, res, recvbuf, count);
			}
		}
		else
		{
			RankRecv = newProcRank&(~mask);
			RankRecv = (RankRecv + root) % ProcNum;
			MPI_Send(recvbuf, count, type, RankRecv, 0, comm);
			break;
		}
		mask = mask << 1;
	}
	if (ProcRank != root)
	{
		free(recvbuf);
	}

	if (ProcRank == root)
	{
		free(res);
	}

	return 0;
}

int* NodePosition(int ProcNum, int root, int ProcRank)
{
	int* posProcess = new int[3];
	int left = -1;
	int right = -1;
	int parent = -1;
	int size = ProcNum;
	int node = size - 1;
	while(1)
	{
		int rsize = size / 2; //число узлов в правом поддереве
		int lsize = size - rsize - 1;
		left = right = -1;
		if (size > 1)
		{
			left = node - 1;
			if (size > 2)
			{
				right = node - lsize - 1;
			}
		}
		if (ProcRank == node)
		{
			break;
		}
		parent = node;
		if (ProcRank > right)
		{
			node = left;
			size = lsize;
		}
		else {
			node = right;
			size = rsize;
		}
	}
	posProcess[0] = parent;
	posProcess[1] = left;
	posProcess[2] = right;
	return posProcess;

}

int ReduceLeftRight(void* sendbuf, void* recvbuf, void* res, int count, MPI_Datatype type, int left, int right, int root,MPI_Op op, MPI_Comm comm)
{
	MPI_Status status;
	memcpy(recvbuf, sendbuf, count*stype(type));
	MPI_Recv(res, count, type, left, 0, comm, &status);
	Operation(type, op, res, recvbuf, recvbuf, count);
	if (right >= 0)
	{
		MPI_Recv(res, count, type, right, 0, comm, &status);
		Operation(type, op, res, recvbuf, recvbuf, count);
	}
	return 0;
}

int NEW_TREE_Reduce_Binary_Tree(void *sendbuf, void *recvbuf, int count, MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
{
	MPI_Status status;
	int ProcRank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int RankRecv, RankSend;
	int* posProcess = NodePosition(ProcNum, root, ProcRank);
	if ((posProcess[1] < 0) && (posProcess[2] < 0))
	{
		MPI_Send(sendbuf, count, type, posProcess[0], 0, comm);
	}
	else
	{
		void* res = malloc(count * stype(type));
		if (ProcRank == ProcNum - 1)
		{
			if (root == ProcNum - 1)
			{
				ReduceLeftRight(sendbuf, recvbuf, res, count, type, posProcess[1], posProcess[2], root, op, comm);
			}
			else
			{
				recvbuf = malloc(count * stype(type));
				ReduceLeftRight(sendbuf, recvbuf, res, count, type, posProcess[1], posProcess[2], root, op, comm);
				MPI_Send(recvbuf, count, type, root, 0, comm);
				free(recvbuf);
			}

		}
		else
		{
			if (ProcRank == root)
			{
				ReduceLeftRight(sendbuf, recvbuf, res, count, type, posProcess[1], posProcess[2], root, op, comm);
				MPI_Send(recvbuf, count, type, posProcess[0], 0, comm);
			}
			else
			{
				recvbuf = malloc(count * stype(type));
				ReduceLeftRight(sendbuf, recvbuf, res, count, type, posProcess[1], posProcess[2], root, op, comm);
				MPI_Send(recvbuf, count, type, posProcess[0], 0, comm);
				free(recvbuf);
			}

		}
		free(res);
	}
	if ((ProcRank == root) && (root != ProcNum - 1))
	{
		MPI_Recv(recvbuf, count, type, ProcNum - 1, 0, comm, &status);
	}
	return 0;
}
