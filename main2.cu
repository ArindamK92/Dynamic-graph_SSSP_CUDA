/*extern "C" {} *///it will instruct the compiler to expect C linkage for your C functions, not C++ linkage.
#include <thrust/find.h>
#include <thrust/device_vector.h>
#include <thrust/count.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/device_free.h>
#include <stdio.h>
//#include "all_structures.h"
#include "all_structure_cuda.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>


#include<vector>
#include <chrono> 


#define THREADS_PER_BLOCK 1024 //we can change it

using namespace std;
using namespace std::chrono;

__global__ void initialize(int nodes, int src, RT_Vertex* SSSP, int* stencil)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < nodes)
	{
		if (index == src) { SSSP[index].Root = -1; }
		else { SSSP[index].Root = index; }
		SSSP[index].Dist = 0.0;
		stencil[index] = index;
	}
}

//__global__ void create_tree(Colwt2* cuda_adjlist_full_X, int start, RT_Vertex* SSSP, int src, int* d_affectedPointer, int numberofCudaThread)
//{
//	int index = threadIdx.x + blockIdx.x * blockDim.x;
//	int number_CudaThread = numberofCudaThread;
//	int flag = 0;
//
//	if (index < number_CudaThread)
//	{
//		/*printf("source: %d", src);*/
//		int y = cuda_adjlist_full_X[index + start].col;
//		/*printf("y: %d", y);*/
//		double mywt = cuda_adjlist_full_X[index + start].wt;
//		if (mywt == -1) { flag = 1; }//invalid edge
//		if (SSSP[y].Root == -1) { flag = 1; }
//		if (flag == 0)
//		{
//			SSSP[y].Parent = src; //mark the parent
//			SSSP[y].EDGwt = mywt; //mark the edgewt
//			SSSP[y].Root = SSSP[src].Root;
//			SSSP[y].Dist = SSSP[src].Dist + mywt;
//			d_affectedPointer[y] = 1;
//		}
//		/*printf("end if***");*/
//
//	}
//
//}

__global__ void create_tree2(Colwt2* cuda_adjlist_full_X, int* d_colStartPtr_X, RT_Vertex* SSSP, int* d_affectedPointer, int* change_d, int nodes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int flag = 0;
	
	if (index < nodes && d_affectedPointer[index] == 1)
	{
		/*printf("index: %d", index);
		printf("affected ptr value: %d", d_affectedPointer[index]);*/
		int adjNodestartIndex = d_colStartPtr_X[index];
		/*printf("y: %d", y);*/
		for (int k = 0; k < d_colStartPtr_X[index + 1] - d_colStartPtr_X[index]; k++)
		{
			flag = 0;
			int y = cuda_adjlist_full_X[adjNodestartIndex + k].col;
			double mywt = cuda_adjlist_full_X[adjNodestartIndex + k].wt;
			/*printf("y: %d", y);*/
			/*printf("SSSP[y].Root:%d,SSSP[y].Parent:%d, SSSP[y].Root:%d,SSSP[y].Dist:%f", SSSP[y].Root, SSSP[y].Parent, SSSP[y].Root, SSSP[y].Dist);*/
			if (mywt == -1) {
				/*printf("check 1.1");*/
				flag = 1; }//invalid edge
			if (SSSP[y].Root == -1) { /*printf("check 1.2");*/ flag = 1; }
			if (SSSP[y].Root != y) { /*printf("check 1.3");*/ flag = 1; }
			if (flag == 0)
			{
				SSSP[y].Parent = index; //mark the parent
				SSSP[y].EDGwt = mywt; //mark the edgewt
				SSSP[y].Root = SSSP[index].Root;
				SSSP[y].Dist = SSSP[index].Dist + mywt;
				/*printf("mywt: %f", mywt);*/
				d_affectedPointer[y] = 1;
			}
		}
		d_affectedPointer[index] = 0;
		change_d[0] = 1;
		/*printf("end if***");*/

	}

}

struct is_affected
{
	__host__ __device__
		bool operator()(const int x)
	{
		return (x == 1);
	}
};

__global__ void initializeUpdatedDist(double* d_UpdatedDist, RT_Vertex* SSSP, int X_size)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{
		d_UpdatedDist[index] = SSSP[index].Dist;
	}
}

__global__ void initializeEdgedone(int* Edgedone, int totalChange)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < totalChange)
	{
		Edgedone[index] = -1;
	}
}

__global__ void insertDeleteEdge(xEdge_cuda* allChange_cuda, int* Edgedone, RT_Vertex* SSSP, int numS, int X_size, int* d_colStartPtr_X, Colwt2* cuda_adjlist_full_X, double* d_UpdatedDist, double inf, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < numS)
	{
		int node_1 = allChange_cuda[index].node1;
		int node_2 = allChange_cuda[index].node2;
		double edge_weight = allChange_cuda[index].edge_wt;

		if (node_1 > X_size || allChange_cuda[index].node2 > X_size)
		{
			Edgedone[index] = 0; //mark to not add
		}

		if (SSSP[node_1].Root != SSSP[node_2].Root)
		{
			Edgedone[index] = 0; //mark to not add
		}

		if (allChange_cuda[index].inst == 1)  //check x
		{
			//Check if edge exists--then dont insert 
			for (int k = 0; k < d_colStartPtr_X[node_1 + 1] - d_colStartPtr_X[node_1]; k++)
			{
				int myn = cuda_adjlist_full_X[d_colStartPtr_X[node_1] + k].col;
				double mywt = cuda_adjlist_full_X[d_colStartPtr_X[node_1] + k].wt; //check. added recently 01-15-20
				//****need check
				if (myn == node_2 && mywt <= edge_weight && mywt != -1)
				{
					Edgedone[index] = 0;
					break;
				}

			}//end of for
		}

		if (allChange_cuda[index].inst == 1 && Edgedone[index] != 0)
		{
			//We check the distances based on updateddist, to cull some insertion edges
			//In case of conflicts, actual distance remains correct

				//Default is remainder edge
			Edgedone[index] = 2;
			//Check twice once for  n1->n2 and once for n2->n1
			for (int yy = 0; yy < 2; yy++)
			{
				int node1, node2;
				if (yy == 0)
				{
					node1 = node_1;
					node2 = node_2;
				}
				else
				{
					node1 = node_2;
					node2 = node_1;
				}

				//  printf("%d:%f:::%d::%f:::%f \n", node1, UpdatedDist[node1],node2, UpdatedDist[node2], mye.edge_wt);
		  //Check whether node1 is relaxed
				if (d_UpdatedDist[node2] > d_UpdatedDist[node1] + edge_weight)
				{
					//Update Parent and EdgeWt
					SSSP[node2].Parent = node1;
					SSSP[node2].EDGwt = edge_weight;
					d_UpdatedDist[node2] = d_UpdatedDist[node1] + edge_weight;
					SSSP[node2].Update = true;
					/*printf("@@@@node: %d, parent: %d, dist: %f", node2, SSSP[node2].Parent, d_UpdatedDist[node2]);*/
					//Mark Edge to be added--node1 updated
					Edgedone[index] = 1;
					break;
				}

			}//end of for

		}//end of if insert

		//Deletion case
		//in case of deletion we don't update d_UpdatedDist
		if (allChange_cuda[index].inst == 0 && Edgedone[index] != 0)  //if deleted
		{
			Edgedone[index] = 3;
			//Check if edge exists in the tree
				//this will happen if node1 is parentof node or vice-versa
			bool iskeyedge = false;

			// printf("XXX:%d:%d \n",mye.node1, mye.node2 );

					 //Mark edge as deleted
			if (SSSP[node_1].Parent == node_2)
			{
				//printf("YYY:%d:%d \n",mye.node1, mye.node2 );
				SSSP[node_1].EDGwt = inf;
				SSSP[node_1].Update = true;
				iskeyedge = true;
				/*d_UpdatedDist[node_1] = inf;*/ //check. added recently 01-15-20
			}
			else {
				//Mark edge as deleted
				if (SSSP[node_2].Parent == node_1)
				{
					// printf("ZZZ:%d:%d \n",mye.node1, mye.node2 );
					SSSP[node_2].EDGwt = inf;
					SSSP[node_2].Update = true;
					iskeyedge = true;
					/*d_UpdatedDist[node_2] = inf;*/ //check. added recently 01-15-20
				}
			}


			//If  Key Edge Delete from key edges
		   //Set weights to -1;
			if (iskeyedge)
			{

				for (int k = 0; k < d_colStartPtr_X[node_1 + 1] - d_colStartPtr_X[node_1]; k++)
				{
					////TEPS:
					//*te = *te + 1;
					int myn = cuda_adjlist_full_X[d_colStartPtr_X[node_1] + k].col;
					if (myn == node_2)
					{
						cuda_adjlist_full_X[d_colStartPtr_X[node_1] + k].wt = -1; //set wt -1 in adj list of old sssp
						break;
					}

				}//end of for

				for (int k = 0; k < d_colStartPtr_X[node_2 + 1] - d_colStartPtr_X[node_2]; k++)
				{
					////TEPS:
					//*te = *te + 1;
					int myn = cuda_adjlist_full_X[d_colStartPtr_X[node_2] + k].col;
					if (myn == node_1)
					{
						cuda_adjlist_full_X[d_colStartPtr_X[node_2] + k].wt = -1; //set wt -1 in adj list of old sssp
						break;
					}

				}
			}//end of if


			/*else      // check. recently added 24-01-2020. The below part is required for all as we consider full graph
			{*/

			for (int k = 0; k < colStartPtr_R[node_1 + 1] - colStartPtr_R[node_1]; k++)
			{
				int myn = cuda_adjlist_full_R[colStartPtr_R[node_1] + k].col;
				if (myn == node_2)
				{
					cuda_adjlist_full_R[colStartPtr_R[node_1] + k].wt = -1;
					break;
				}

			}//end of for

			for (int k = 0; k < colStartPtr_R[node_2 + 1] - colStartPtr_R[node_2]; k++)
			{
				int myn = cuda_adjlist_full_R[colStartPtr_R[node_2] + k].col;
				if (myn == node_1)
				{
					cuda_adjlist_full_R[colStartPtr_R[node_2] + k].wt = -1;
					break;
				}

			}//end of for

		//}//end of if

		}//end of else if deleted
	}
}


__global__ void checkInsertedEdges(int numS, int* Edgedone, double* d_UpdatedDist, xEdge_cuda* allChange_cuda, RT_Vertex* SSSP, int* change_d)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < numS)
	{

		if (Edgedone[index] == 1) //Edgedone will be 1 when edge is marked to be inserted
		{

			//get the edge
			int node_1 = allChange_cuda[index].node1;
			int node_2 = allChange_cuda[index].node2;
			double edgeWeight = allChange_cuda[index].edge_wt;
			//reset it to 0
			Edgedone[index] = 0;


			int node1, node2;
			if (d_UpdatedDist[node_1] > d_UpdatedDist[node_2])
			{
				node1 = node_1;
				node2 = node_2;
			}
			else
			{
				node1 = node_2;
				node2 = node_1;
			}

			//Check if some other edge was added--mark edge to be added //check x
			if (d_UpdatedDist[node1] > d_UpdatedDist[node2] + edgeWeight)
			{
				Edgedone[index] = 1;
			}

			//Check if correct edge wt was written--mark edge to be added //check x
			if ((SSSP[node1].Parent == node2) && (SSSP[node1].EDGwt > edgeWeight))
			{
				Edgedone[index] = 1;
			}


			if (Edgedone[index] == 1)
			{
				//Update Parent and EdgeWt
				SSSP[node1].Parent = node2;
				SSSP[node1].EDGwt = edgeWeight;
				d_UpdatedDist[node1] = d_UpdatedDist[SSSP[node1].Parent] + SSSP[node1].EDGwt;
				SSSP[node2].Update = true;
				change_d[0] = 1;
			}


		}//end of if
	}
}

__global__ void updateDistance(int X_size, RT_Vertex* SSSP, double* d_UpdatedDist, double inf)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{
		//do not update source node
		int px = SSSP[index].Parent;
		int flag = 0;
		if (SSSP[index].Parent == -1) { flag = 1; }


		if (flag != 1 && index == SSSP[px].Parent)
		{
			printf("DP: %d:%d %d:%d \n", index, SSSP[index].Parent, px, SSSP[px].Parent);
		}

		//For deletion case
		if (flag != 1 && SSSP[index].EDGwt == inf)
		{
			SSSP[index].Dist = inf;
			SSSP[index].Update = true;
			flag = 1;
		}

		//for insertion case
		if (flag != 1 && SSSP[index].Dist > d_UpdatedDist[index])
		{
			SSSP[index].Dist = d_UpdatedDist[index];
			/*printf("In updateDistance:  index: %d, dist:%f\n", index, SSSP[index].Dist);*/
			SSSP[index].Update = true;
		}

	}
}

__global__ void initializeUpdatedDistOldDist(double* d_UpdatedDist, double* d_OldUpdate, RT_Vertex* SSSP, int X_size)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{
		d_UpdatedDist[index] = SSSP[index].Dist; //this will fill up the d_UpdatedDist values for deletion case also
		d_OldUpdate[index] = SSSP[index].Dist;
		/*printf("****Inside initializeUpdatedDistOldDist: %d edge weight: %f parent: %d dist: %f\n", index, SSSP[index].EDGwt, SSSP[index].Parent, SSSP[index].Dist);*/
	}
}


//revised function //check. recently added function. 24-01-2020
__global__ void updateNeighbors(double* d_UpdatedDist, RT_Vertex* SSSP, int X_size, /*int* d_mychange,*/ int* d_colStartPtr_X, Colwt2* cuda_adjlist_full_X, double inf, int* change_d, int its, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{

		//If i is updated--update its neighbors as required
		if (SSSP[index].Update)
		{
			SSSP[index].Update = false;
			for (int j = 0; j < colStartPtr_R[index + 1] - colStartPtr_R[index]; j++)
			{
				int myn = cuda_adjlist_full_R[colStartPtr_R[index] + j].col;
				double mywt = cuda_adjlist_full_R[colStartPtr_R[index] + j].wt;

				if (SSSP[myn].EDGwt < mywt && SSSP[myn].Parent == index) //check if we have taken an edge with lower weight from the changeEdge set. if yes then don't update edgeweight
				{
					mywt = SSSP[myn].EDGwt;
				}
				//check if edge is deleted
				if (mywt < 0) { continue; } //if mywt = -1, that means node was deleted
				if (SSSP[index].Dist == inf)
				{
					/*printf("$$$$: %d edge weight: %f parent: %d dist: %f\n", index, SSSP[index].EDGwt, SSSP[index].Parent, SSSP[index].Dist);
					printf("$$$$myn: %d myn weight: %f parent: %d dist: %f\n", myn, mywt, SSSP[myn].Parent, SSSP[myn].Dist);
					*/
					if (myn == SSSP[index].Parent)
					{
						continue;
					}
					if (SSSP[myn].Parent == index)
					{
						d_UpdatedDist[myn] = inf;
						SSSP[myn].Dist = inf;
						SSSP[myn].Update = true;
						change_d[0] = 1;
						/*printf("&&&&index: %d edge weight: %f parent: %d dist: %f\n", index, SSSP[index].EDGwt, SSSP[index].Parent, SSSP[index].Dist);
						printf("&&&&myn: %d myn weight: %f parent: %d dist: %f\n", myn, mywt, SSSP[myn].Parent, SSSP[myn].Dist);
						*/
						continue;
					}
					else {
						if (SSSP[myn].Dist != inf)
						{
							d_UpdatedDist[index] = d_UpdatedDist[myn] + mywt;
							SSSP[index].Dist = d_UpdatedDist[myn] + mywt;
							SSSP[index].Parent = myn;
							SSSP[index].EDGwt = mywt;
							SSSP[index].Update = true;
							change_d[0] = 1;
							/*printf("++++index: %d edge weight: %f parent: %d dist: %f\n", index, SSSP[index].EDGwt, SSSP[index].Parent, SSSP[index].Dist);
							printf("++++myn: %d myn weight: %f parent: %d dist: %f\n", myn, mywt, SSSP[myn].Parent, SSSP[myn].Dist);*/

							continue;
						}
					}
				}
				if (SSSP[index].Dist != inf)
				{
					/*printf("Not inf: index: %d edge weight: %f parent: %d dist: %f\n", index, SSSP[index].EDGwt, SSSP[index].Parent, SSSP[index].Dist);
					printf("Not infmyn: %d myn weight: %f parent: %d dist: %f\n", myn, mywt, SSSP[myn].Parent, SSSP[myn].Dist);
					*/
					if (SSSP[myn].Dist == inf)
					{
						if (SSSP[index].Parent != myn)
						{
							d_UpdatedDist[myn] = d_UpdatedDist[index] + mywt;
							SSSP[myn].Dist = SSSP[index].Dist + mywt;
							SSSP[myn].EDGwt = mywt;
							SSSP[myn].Update = true;
							SSSP[myn].Parent = index;
							change_d[0] = 1;
							continue;
						}
						else {
							//don't do anything if myn is parent of index node
							continue;
						}

					}
					if (d_UpdatedDist[myn] > d_UpdatedDist[index] + mywt) //update both cases where parent of myn == index or parent of myn != index
					{
						//if (SSSP[myn].EDGwt < mywt && SSSP[myn].Parent == index) //check if we have taken an edge with lower weight from the changeEdge set. if yes then don't update edgeweight
						//{
						//	mywt = SSSP[myn].EDGwt;
						//}
						d_UpdatedDist[myn] = d_UpdatedDist[index] + mywt;
						SSSP[myn].Dist = d_UpdatedDist[index] + mywt;
						SSSP[myn].Update = true;
						SSSP[myn].Parent = index;
						change_d[0] = 1;
						/*printf("Not inf: index: %d edge weight: %f parent: %d dist: %f\n", index, SSSP[index].EDGwt, SSSP[index].Parent, SSSP[index].Dist);
						printf("Not infmyn: %d myn weight: %f parent: %d dist: %f\n", myn, mywt, SSSP[myn].Parent, SSSP[myn].Dist);
						*/
						continue;
					}
					else
					{
						if (SSSP[myn].Parent == index)
						{
							d_UpdatedDist[myn] = d_UpdatedDist[index] + mywt;
							SSSP[myn].Dist = d_UpdatedDist[index] + mywt;
							SSSP[myn].Update = true;
							/*SSSP[myn].Parent = index;*/ //Parent will remain same
							change_d[0] = 1;
							continue;
						}
						if ((d_UpdatedDist[index] > d_UpdatedDist[myn] + mywt) /*&& (SSSP[myn].Parent != index)*/)
						{
							d_UpdatedDist[index] = d_UpdatedDist[myn] + mywt;
							SSSP[index].Dist = d_UpdatedDist[myn] + mywt;
							SSSP[index].Update = true;
							SSSP[index].Parent = myn;
							change_d[0] = 1;
						}
					}
				}
			}
		}
	}
}

__global__ void checkIfDistUpdated(int X_size, double* d_OldUpdate, double* d_UpdatedDist, RT_Vertex* SSSP)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{
		if (d_OldUpdate[index] != d_UpdatedDist[index])
		{
			d_OldUpdate[index] = d_UpdatedDist[index];
			SSSP[index].Update = true;
		}
		else { SSSP[index].Update = false; }
	}
}

__global__ void updateDistanceFinal(int X_size, double* d_UpdatedDist, RT_Vertex* SSSP, double inf)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{
		int flag = 0;
		//do not update parent
		if (SSSP[index].Parent == -1) { flag = 1; }

		if (flag == 0)
		{
			int px = SSSP[index].Parent;
			if (px > -1)
			{
				//printf("XX %d :%d \n", i, px);  
				if (index == SSSP[px].Parent)
				{
					printf("BBP %d %d \n", index, px);
				}
			}
			if (d_UpdatedDist[index] >= inf)
			{
				SSSP[index].Dist = inf;
			}
			else
			{
				SSSP[index].Dist = d_UpdatedDist[SSSP[index].Parent] + SSSP[index].EDGwt;
				//printf("Check 2. index: %d dist: %f, parent dist:%f, edgewt: %f \n", index, SSSP[index].Dist, d_UpdatedDist[SSSP[index].Parent], SSSP[index].EDGwt); //Test 23-01-2020
			}
		}
	}
}

void edge_update(int* totalChange, int* X_size, int* SSSP_size, xEdge_cuda* allChange_cuda, Colwt2* cuda_adjlist_full_X, int* colStartPtr_X, RT_Vertex* SSSP, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R, int* te, int* nodes);
void rest_update(int* X_size, Colwt2* cuda_adjlist_full_X, int* colStartPtr_X, RT_Vertex* SSSP, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R, int* nodes);

int main(int argc, char* argv[]) {

	double startx, endx, starty, endy;
	/*double inf = std::numeric_limits<double>::infinity();*/

	/***** Preprocessing to Graph (GUI) ***********/
	int nodes, edges;
	cudaError_t cudaStatus;

	printf("Enter number of total nodes: ");
	scanf("%d", &nodes);
	printf("Enter number of total edges: ");
	scanf("%d", &edges);
	printf("check 0");



	/*** Read Full Graph ***/
	int* colStartPtr_R;
	cout << "success 1" << endl;
	cudaStatus = cudaMallocManaged((void**)&colStartPtr_R, (nodes + 1) * sizeof(int)); //we take nodes +1 to store the start ptr of the first row 
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		/*goto Error;*/
	}
	int total_adjmatrix_size_R = edges * 2; //e.g.= (0 1 wt1), (1 0 wt1) both are same edge, but both will be there
	Colwt2* cuda_adjlist_full_R;
	cudaStatus = cudaMallocManaged(&cuda_adjlist_full_R, total_adjmatrix_size_R * sizeof(Colwt2));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		/*goto Error;*/
	}
	printf("check 1");

	//use below for direct path
	/*string file1 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\fullGraph.txt";
	char* cstr1 = &file1[0];
	readin_graphU(&R, nodes, cstr1);*/

	//use below code if we use pass file name as argument
	//readin_graphU(&R, nodes, argv[1]);


	//use below code to pass the file name as relative path.
	//**keep the files in the same folder
	string file1 = "./fullGraph.txt";
	char* cstr1 = &file1[0];
	/*readin_graphU(&R, nodes, cstr1);*/
	/*readin_graphU2(&R, nodes, cstr1);*/ //new running
	readin_graphU4(colStartPtr_R, cuda_adjlist_full_R, cstr1, &nodes);
	/*readin_network2(&R, cstr1, -1);
	*/
	cout << "success 2" << endl;
	/*for (int i = 0; i < nodes + 1; i++)
	{
		cout <<"*******Row start ptr: "<< colStartPtr_R[i] << endl;
			for (int y = colStartPtr_R[i]; y< colStartPtr_R[i+1];y++)
			{
				cout <<"node: "<< cuda_adjlist_full_R[y].col <<"weight: "<< cuda_adjlist_full_R[y].wt << endl;
			}
	}*/

	/*** Finished Reading Full graph **/

	/*** Read the input SSSP ***/
	int* colStartPtr_X = (int*)malloc((nodes + 1) * sizeof(int));//we take nodes +1 to store the start ptr of the first row 
	int total_adjmatrix_size_X = (nodes - 1) * 2; //maximum number of edges in SSSP tree = nodes - 1. Each edge will take 2 places in adjacent list
	Colwt2* cuda_adjlist_full_X;
	cudaStatus = cudaMallocManaged(&cuda_adjlist_full_X, total_adjmatrix_size_X * sizeof(Colwt2));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed at SSSP file");
		/*goto Error;*/
	}
	cout << "success 3" << endl;
	//use below for direct path
	/*string file2 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\SSSP.txt";
	char* cstr2 = &file2[0];
	readin_network(&X, cstr2, -1);*/

	//use below code if we use pass file name as argument
	/*readin_network(&X, argv[2], -1);*/

	//use below code to pass the file name as relative path.
	//**keep the files in the same folder
	string file2 = "./SSSP.txt";
	char* cstr2 = &file2[0];
	readin_graphU4(colStartPtr_X, cuda_adjlist_full_X, cstr2, &nodes);
	cout << "success 4" << endl;
	/*** Finished Reading input SSSP **/

	/*** Read the change file ***/
	//There will be a list for inserts and a list for delete
	vector<xEdge> allChange;
	allChange.clear();

	/*** Read set of Changed Edges ***/
   //use below for direct path
   /*string file3 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\changeEdges.txt";
   char* cstr3 = &file3[0];
   readin_changes(cstr3, &allChange);*/

   //use below code if we use pass file name as argument
   /*readin_changes(argv[3], &allChange);*/

   //use below code to pass the file name as relative path.
	string file3 = "./changeEdges.txt";
	char* cstr3 = &file3[0];
	readin_changes(cstr3, &allChange);
	cout << "success 5" << endl;
	//new addition
	xEdge_cuda* allChange_cuda;
	int totalChange = allChange.size();
	cudaStatus = cudaMallocManaged(&allChange_cuda, totalChange * sizeof(xEdge_cuda));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed at changeEdge");
		/*goto Error;*/
	}
	for (int i = 0; i < totalChange; i++)
	{
		allChange_cuda[i].node1 = allChange.at(i).theEdge.node1;
		allChange_cuda[i].node2 = allChange.at(i).theEdge.node2;
		allChange_cuda[i].edge_wt = allChange.at(i).theEdge.edge_wt;
		allChange_cuda[i].inst = allChange.at(i).inst;
	}
	/*** Finished Reading Changed Edges **/

	//Initializing  Rooted Tree
	RT_Vertex* SSSP;
	cudaStatus = cudaMallocManaged(&SSSP, nodes * sizeof(RT_Vertex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed at SSSP structure");
		/*goto Error;*/
	}
	int* stencil; //stencil is used for tracking which node is being affected. 
	/*cudaMallocManaged(&stencil, nodes * sizeof(int));*/
	cudaStatus = cudaMalloc((void**)&stencil, nodes * sizeof(int));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed at SSSP stencil");
		/*goto Error;*/
	}
	/*int* stencil_c = new int[nodes];*/
	int* stencil_c = (int*)calloc(nodes, sizeof(int));
	/*vector<SCC_vertex>SCC;*/
	int graphDirectedUndirectedIndicator = 0; // Should be 1 for SCC, 0 for not SCC. need to modify if we want SCC

	int source;
	printf("Enter source node: ");
	scanf("%d", &source);
	int p;



	if (graphDirectedUndirectedIndicator == 0) {
		int src = source; //the source from which the paths are computed
		initialize << <(nodes / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (nodes, src, SSSP, stencil); //kernet call
		cudaDeviceSynchronize();
		cudaMemcpy(stencil_c, stencil, nodes * sizeof(int), cudaMemcpyDeviceToHost);
		/*for (int i = 0; i < nodes; i++)
		{
			cout << "stencil_c" << stencil_c[i] << endl;
		}*/
		/*for (int i = 0; i < nodes; i++)
		{

			cout <<"dist"<< SSSP->Dist << endl;
			cout <<"wt"<< SSSP->EDGwt << endl;
			cout << "level"<< SSSP->Root << endl;
			cout << "marked"<< SSSP->Parent << endl;
		}*/


		//Code for create_tree:
		//Time calculation
		auto startTime = high_resolution_clock::now();
		int totalAffectedNode; //alias of numberOfAffectedNode

		int* affectedPointer;
		int* d_affectedPointer;
		cudaStatus = cudaMalloc((void**)&d_affectedPointer, nodes * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed at affectedPointer");
			/*goto Error;*/
		}
		affectedPointer = (int*)calloc(nodes, sizeof(int));
		affectedPointer[0] = 1;
		cudaMemcpy(d_affectedPointer, affectedPointer, nodes * sizeof(int), cudaMemcpyHostToDevice);
		/*cudaMallocManaged(&affectedPointer, nodes * sizeof(int));*/

		//new addition
		int* change_d = new int[1];
		int* change = new int[1];
		change[0] = 1;
		cudaMalloc((void**)&change_d, 1 * sizeof(int));
		int* d_colStartPtr_X;
		cudaMalloc((void**)&d_colStartPtr_X, (nodes + 1) * sizeof(int));
		cudaMemcpy(d_colStartPtr_X, colStartPtr_X, (nodes + 1) * sizeof(int), cudaMemcpyHostToDevice);
		while (change[0] == 1)
		{
			change[0] = 0;
			cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);
			create_tree2 << <(nodes / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK>> > (cuda_adjlist_full_X, d_colStartPtr_X, SSSP, d_affectedPointer, change_d, nodes);
			cudaDeviceSynchronize();
			cudaMemcpy(change, change_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);
			/*cout << "change"<< change[0]<<endl;*/

		}


		//totalAffectedNode = 1;
		//int start = 0, end = 0;
		//int* affected_nodes;
		//affected_nodes = (int*)calloc(totalAffectedNode, sizeof(int));
		//affected_nodes[0] = src;
		//cudaStream_t stream1;
		//cudaError_t result;
		//result = cudaStreamCreate(&stream1);
		//try
		//{
		//	while (totalAffectedNode > 0)
		//	{
		//		vector<int> affectedNodeAlias;
		//		for (int i = 0; i < totalAffectedNode; i++)
		//		{
		//			affectedNodeAlias.push_back(affected_nodes[i]);
		//		}
		//		for (int i = 0; i < totalAffectedNode; i++)
		//		{
		//			p = affectedNodeAlias.at(i);
		//			src = p;
		//			start = colStartPtr_X[p];
		//			end = colStartPtr_X[p + 1];
		//			int numberofCudaThread = end - start;

		//			//create_tree method creates the SSSP tree with values stored in cuda_adjlist_full_X.
		//			//This SSSP tree is the input SSSP tree
		//			create_tree << <(numberofCudaThread / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK, 0, stream1 >> > (cuda_adjlist_full_X, start, SSSP, src, d_affectedPointer, numberofCudaThread);
		//		}
		//		cudaStreamSynchronize(stream1);

		//		thrust::device_ptr<int> affectedPointer_alias(d_affectedPointer); // converting raw ptr to device_ptr
		//		cout << "check X.2";
		//		cudaMemcpy(affectedPointer, d_affectedPointer, nodes * sizeof(int), cudaMemcpyDeviceToHost);
		//		cout << "check Y";
		//		thrust::device_vector<int> affectedPointer_vector(affectedPointer_alias, affectedPointer_alias + nodes); //converting device_ptr to device_vector
		//		cout << "check Y.2";
		//		totalAffectedNode = thrust::count(affectedPointer_vector.begin(), affectedPointer_vector.end(), 1); //count the number of affected node
		//		cout << "total affected node" << totalAffectedNode << endl;
		//		affected_nodes = (int*)realloc(affected_nodes, totalAffectedNode * sizeof(int));
		//		/*affectedPointer = thrust::raw_pointer_cast(&affectedPointer_vector[0]);*/
		//		thrust::copy_if(thrust::host, stencil_c, stencil_c + nodes, affectedPointer, affected_nodes, is_affected());
		//		cout << "check Y.3";


		//		//Test code - start
		//		/*cout << "affected nodes " << endl;
		//		for (int i = 0; i < totalAffectedNode; i++)
		//		{
		//			cout << affected_nodes[i] << endl;
		//		}*/
		//		//Test code - end
		//		free(affectedPointer);
		//		cout << "check Y.4";
		//		affectedPointer = (int*)calloc(nodes, sizeof(int));
		//		cudaMemcpy(d_affectedPointer, affectedPointer, nodes * sizeof(int), cudaMemcpyHostToDevice);
		//		cout << "check Y.5";
		//		////clear thrust vector:
		//		//affectedPointer_vector.clear();
		//		//affectedPointer_vector.shrink_to_fit();
		//		/*cudaFree(threadHelpers);*/
		//	}
		//}
		//catch (thrust::system_error e)
		//{
		//	std::cerr << "Error from thrust: " << e.what() << std::endl;
		//	free(affectedPointer);
		//	free(affected_nodes);
		//	cudaFree(d_affectedPointer);
		//}
		free(affectedPointer);
		/*free(affected_nodes);*/
		cudaFree(d_affectedPointer);
		cudaFree(d_colStartPtr_X);

		//Time calculation
		auto stopTime = high_resolution_clock::now();
		// Time calculation
		auto duration = duration_cast<microseconds>(stopTime - startTime);
		cout << "Time taken by create-tree function: "
			<< duration.count() << " microseconds" << endl;

		//test
		/*cout << "input sssp tree" << endl;
		for (int i = 0; i < nodes; i++)
		{
			cout << "node" << i << endl;
			cout << "dist" << SSSP[i].Dist << endl;
			cout << "parent" << SSSP[i].Parent << endl;
			cout << "parent" << SSSP[i].EDGwt << endl;
		}*/
		//test end
		//edge_update function
		//Update the inserted and delted edges in the tree
		int x_size = nodes;
		int SSSP_size = nodes; //considering all nodes are participating in the SSSP
		int te = 0;
		auto startTime1 = high_resolution_clock::now();
		edge_update(&totalChange, &x_size, &SSSP_size, allChange_cuda, cuda_adjlist_full_X, colStartPtr_X, SSSP, cuda_adjlist_full_R, colStartPtr_R, &te, &nodes);
		cout << "after edge_update fn" << endl;
		//Time calculation
		auto stopTime1 = high_resolution_clock::now();
		// Time calculation
		auto duration1 = duration_cast<microseconds>(stopTime1 - startTime1);
		cout << "Time taken by edge_update function: "
			<< duration1.count() << " microseconds" << endl;

		auto startTime2 = high_resolution_clock::now();
		rest_update(&x_size, cuda_adjlist_full_X, colStartPtr_X, SSSP, cuda_adjlist_full_R, colStartPtr_R, &nodes);
		cout << "after rest_update fn" << endl;
		//Time calculation
		auto stopTime2 = high_resolution_clock::now();
		// Time calculation
		auto duration2 = duration_cast<microseconds>(stopTime2 - startTime2);
		cout << "Time taken by rest_update function: "
			<< duration2.count() << " microseconds" << endl;
	}
	else
	{
		//****below code needs modification
		/*SCC.clear();
		readin_SCC(argv[2], &SCC);
		update_SCC(&X, &SCC, &allChange);
		print_network(X);*/
	}
	//Test code start
	/*cout << "SSSP" << endl;
	for (int i = 0; i < nodes; i++)
	{
		cout << "*******" << endl;
		cout << "node" << i << endl << "dist" << SSSP[i].Dist << endl << "parent" << SSSP[i].Parent << endl;
	}*/
	cout << "*******success*******" << endl;
	//Test code end


	cudaFree(colStartPtr_R);
	cudaFree(cuda_adjlist_full_R);
	cudaFree(colStartPtr_X);
	cudaFree(cuda_adjlist_full_X);
	cudaFree(allChange_cuda);
	cudaFree(SSSP);
	cudaFree(stencil);
Error:
	cudaFree(colStartPtr_R);
	return 0;
}



void edge_update(int* totalChange, int* X_size, int* SSSP_size, xEdge_cuda* allChange_cuda, Colwt2* cuda_adjlist_full_X, int* colStartPtr_X, RT_Vertex* SSSP, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R, int* te, int* nodes)
{
	double inf = std::numeric_limits<double>::infinity();
	/*int* Edgedone;*/
	double* UpdatedDist;

	int iter = 0;

	//Mark how the edge is processed
	int* Edgedone;
	cudaMalloc((void**)&Edgedone, (*totalChange) * sizeof(int));
	//initialize Edgedone array with -1
	initializeEdgedone << <((*totalChange) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (Edgedone, *totalChange);
	cudaDeviceSynchronize();

	/*thrust::device_vector<int> Edgedone_ptr(*totalChange);
	thrust::fill(Edgedone_ptr.begin(), Edgedone_ptr.end(), -1);
	int* Edgedone = thrust::raw_pointer_cast(Edgedone_ptr);*/

	//Store the updated distance value
	UpdatedDist = (double*)calloc(*X_size, sizeof(double));
	double* d_UpdatedDist;
	cudaMalloc((void**)&d_UpdatedDist, (*X_size) * sizeof(double));
	cudaMemcpy(d_UpdatedDist, UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyHostToDevice);

	//Initialize with current distance for each node
	initializeUpdatedDist << <((*X_size) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (d_UpdatedDist, SSSP, *X_size);
	cudaDeviceSynchronize();
	/*	cudaMemcpy(UpdatedDist, d_UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyDeviceToHost);*/ //not required


	int numS = *totalChange;
	int* d_colStartPtr_X;
	cudaMalloc((void**)&d_colStartPtr_X, (*nodes + 1) * sizeof(int));
	cudaMemcpy(d_colStartPtr_X, colStartPtr_X, (*nodes + 1) * sizeof(int), cudaMemcpyHostToDevice);

	insertDeleteEdge << < (numS / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (allChange_cuda, Edgedone, SSSP, numS, *X_size, d_colStartPtr_X, cuda_adjlist_full_X, d_UpdatedDist, inf, cuda_adjlist_full_R, colStartPtr_R);
	cudaDeviceSynchronize();


	/*int* Edgedone_c = new int[*totalChange];
	cudaMemcpy(Edgedone_c, Edgedone, *totalChange * sizeof(int), cudaMemcpyDeviceToHost); *///not req.
	/*cudaMemcpy(UpdatedDist, d_UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyDeviceToHost); *///not req.


	//Go over the inserted edges to see if they need to be changed
	int* change_d = new int[1];
	int* change = new int[1];
	change[0] = 1;
	cudaMalloc((void**)&change_d, 1 * sizeof(int));
	/*cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);*/ //recent change
	while (change[0] == 1)
	{
		change[0] = 0;
		cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);
		checkInsertedEdges << < (numS / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (numS, Edgedone, d_UpdatedDist, allChange_cuda, SSSP, change_d);
		cudaDeviceSynchronize();
		cudaMemcpy(change, change_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);
		/*cout << "change"<< change[0]<<endl;*/

	}

	//Update the distances
	 //Initialize with current distance for each node
	updateDistance << < (numS / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (*X_size, SSSP, d_UpdatedDist, inf);
	cudaDeviceSynchronize();


	cudaFree(change_d);
	cudaFree(d_UpdatedDist);
	cudaFree(d_colStartPtr_X);
	free(UpdatedDist);
	return;
}


void rest_update(int* X_size, Colwt2* cuda_adjlist_full_X, int* colStartPtr_X, RT_Vertex* SSSP, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R, int* nodes)
{
	double inf = std::numeric_limits<double>::infinity();


	int its = 0; //number of iterations

	int* change_d = new int[1];
	int* change = new int[1]; //marking whether the connections changed in the iteration
	change[0] = 1;
	cudaMalloc((void**)&change_d, 1 * sizeof(int));
	/*cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);*/ //recent change

	double* UpdatedDist;
	//Store the updated distance value
	UpdatedDist = (double*)calloc(*X_size, sizeof(double));
	double* d_UpdatedDist;
	cudaMalloc((void**)&d_UpdatedDist, (*X_size) * sizeof(double));
	cudaMemcpy(d_UpdatedDist, UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyHostToDevice);


	double* OldUpdate;
	//Store the old updated distance value
	OldUpdate = (double*)calloc(*X_size, sizeof(double));
	double* d_OldUpdate;
	cudaMalloc((void**)&d_OldUpdate, (*X_size) * sizeof(double));
	cudaMemcpy(d_OldUpdate, OldUpdate, (*X_size) * sizeof(double), cudaMemcpyHostToDevice);


	//int* mychange;
	////Store the old updated distance value
	//mychange = (int*)calloc(*X_size, sizeof(int));
	//int* d_mychange;
	//cudaMalloc((void**)&d_mychange, (*X_size) * sizeof(int));
	//cudaMemcpy(d_mychange, mychange, (*X_size) * sizeof(int), cudaMemcpyHostToDevice);

	int* d_colStartPtr_X;
	cudaMalloc((void**)&d_colStartPtr_X, (*nodes + 1) * sizeof(int));
	cudaMemcpy(d_colStartPtr_X, colStartPtr_X, (*nodes + 1) * sizeof(int), cudaMemcpyHostToDevice);


	//Initialize with current distance for each node
	initializeUpdatedDistOldDist << <((*X_size) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (d_UpdatedDist, d_OldUpdate, SSSP, *X_size);
	cudaDeviceSynchronize();


	int iter = 0;
	while (change[0] == 1 && its < 70)
	{
		printf("Iteration:%d \n", its);

		change[0] = 0;
		cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);
		updateNeighbors << <((*X_size) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (d_UpdatedDist, SSSP, *X_size, /*d_mychange,*/ d_colStartPtr_X, cuda_adjlist_full_X, inf, change_d, its, cuda_adjlist_full_R, colStartPtr_R);
		cudaDeviceSynchronize();
		cudaMemcpy(change, change_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);

		//Test code start
		/*cudaMemcpy(UpdatedDist, d_UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyDeviceToHost);
		for (int i = 0; i < *X_size; i++)
		{
			cout << "UpdatedDist: " << UpdatedDist[i] << endl;
		}*/
		//Test code end


	//Check if distance was updated
		checkIfDistUpdated << <((*X_size) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (*X_size, d_OldUpdate, d_UpdatedDist, SSSP);
		cudaDeviceSynchronize();
		its++;
	}//end of while
	printf("Total Iterations to Converge %d \n", its);

	//Update the distances
	//Initialize with current distance for each node
	updateDistanceFinal << <((*X_size) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (*X_size, d_UpdatedDist, SSSP, inf);
	cudaDeviceSynchronize();


	free(UpdatedDist);
	free(OldUpdate);
	/*free(mychange);*/
	cudaFree(change_d);
	cudaFree(d_UpdatedDist);
	cudaFree(d_OldUpdate);
	/*cudaFree(d_mychange);*/
	cudaFree(d_colStartPtr_X);

	return;
}

