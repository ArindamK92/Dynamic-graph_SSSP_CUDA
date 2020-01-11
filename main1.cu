#include <stdio.h>
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<vector>
#include <set>
#include<queue>
#include <chrono> 
#include <thread>
#include<queue>
#include<stack>
#include<list>
#include "all_structures.h"
#include <thrust/find.h>
#include <thrust/device_vector.h>
#include <thrust/count.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#define THREADS_PER_BLOCK 10 //we can change it

using namespace std;
using namespace std::chrono;

__global__ void initialize(int nodes, int src, RT_Vertex* SSSP, int* stencil)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < nodes)
	{
		if (index == src) { SSSP[index].Root = -1; }
		else { SSSP[index].Root = index; }
		SSSP[index].Level = 0;
		SSSP[index].Dist = 0.0;
		stencil[index] = index;
	}
}

__global__ void create_tree(Colwt2* cuda_adjlist_full_X, int start, RT_Vertex* SSSP, int src, int* d_affectedPointer, int numberofCudaThread)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int number_CudaThread = numberofCudaThread;
	
	if (index < number_CudaThread)
	{
		printf("source: %d", src);
		int y = cuda_adjlist_full_X[index + start].col;
		printf("y: %d", y);
		double mywt = cuda_adjlist_full_X[index + start].wt;

		SSSP[y].Parent = src; //mark the parent
		SSSP[y].EDGwt = mywt; //mark the edgewt
		SSSP[y].Level = 1; //mark the Level
		SSSP[y].Root = -1;
		SSSP[y].Dist = SSSP[src].Dist + mywt;
		d_affectedPointer[y] = 1;
		printf("end if***");
		
	}

}

//__global__ void create_tree2(Colwt2* cuda_adjlist_full_X, RT_Vertex* SSSP, int* d_affectedPointer, ThreadHelper* threadHelpers, int totalNumberofCudaThread)
//{
//	int index = threadIdx.x + blockIdx.x * blockDim.x;
//	int number_CudaThread = totalNumberofCudaThread;
//	int i =0;
//	
//
//	if (index < number_CudaThread)
//	{
//		while (index >= threadHelpers[i+1].offset)
//		{
//			i = i + 1;
//		}
//		int src = 0;
//		int start = 0;
//		src = threadHelpers[i].src;
//		start = threadHelpers[i].start;
//		printf("source: %d", src);
//		int y = cuda_adjlist_full_X[index - threadHelpers[i].offset + start].col;
//		printf("y: %d", y);
//		double mywt = cuda_adjlist_full_X[index - threadHelpers[i].offset + start].wt;
//
//		SSSP[y].Parent = src; //mark the parent
//		SSSP[y].EDGwt = mywt; //mark the edgewt
//		SSSP[y].Level = 1; //mark the Level
//		SSSP[y].Root = -1;
//		SSSP[y].Dist = SSSP[src].Dist + mywt;
//		d_affectedPointer[y] = 1;
//		printf("end if***");
//
//	}
//
//}

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

		if (node_1 > X_size || allChange_cuda[index].node2 >  X_size)
		{
			Edgedone[index] = 0; //mark to not add
		}

		if (SSSP[node_1].Root != SSSP[node_2].Root)
		{
			Edgedone[index] = 0; //mark to not add
		}

		if (allChange_cuda[index].inst == 1)
		{
			//Check if edge exists--then dont insert 
			for (int k = 0; k < d_colStartPtr_X[node_1 + 1] - d_colStartPtr_X[node_1]; k++)
			{
				////TEPS:
				//*te = *te + 1;
				int myn = cuda_adjlist_full_X[d_colStartPtr_X[node_1] + k].col;
				if (myn == node_2)
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
					d_UpdatedDist[node2] = d_UpdatedDist[SSSP[node2].Parent] + SSSP[node2].EDGwt;
					SSSP[node2].Update = true;

					//Mark Edge to be added--node1 updated
					Edgedone[index] = 1;
					break;
				}

			}//end of for

		}//end of if insert

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
			}
			else {
				//Mark edge as deleted
				if (SSSP[node_2].Parent == node_1)
				{
					// printf("ZZZ:%d:%d \n",mye.node1, mye.node2 );
					SSSP[node_2].EDGwt = inf;
					SSSP[node_2].Update = true;
					iskeyedge = true;
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
						cuda_adjlist_full_X[d_colStartPtr_X[node_1] + k].wt = -1;
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
						cuda_adjlist_full_X[d_colStartPtr_X[node_2] + k].wt = -1;
						break;
					}

				}
			}//end of if

				//If not Key Edge Delete from remainder edges
				//Set weights to -1;
			else
			{
				//==From remainder edges
				//TBD: Only delete for mye.node1< mye.node2

				for (int k = 0; k < colStartPtr_R[node_1 + 1] - colStartPtr_R[node_1]; k++)
				{
					////TEPS:
					//*te = *te + 1;
					int myn = cuda_adjlist_full_R[colStartPtr_R[node_1] + k].col;
					if (myn == node_2)
					{
						cuda_adjlist_full_R[colStartPtr_R[node_1] + k].wt = -1;
						break;
					}

				}//end of for

				for (int k = 0; k < colStartPtr_R[node_2 + 1] - colStartPtr_R[node_2]; k++)
				{
					////TEPS:
					//*te = *te + 1;
					int myn = cuda_adjlist_full_R[colStartPtr_R[node_2] + k].col;
					if (myn == node_1)
					{
						cuda_adjlist_full_R[colStartPtr_R[node_2] + k].wt = -1;
						break;
					}

				}//end of for

			}//end of if

		}//end of else if deleted
	}
}


__global__ void checkInsertedEdges(int numS, int* Edgedone, double* d_UpdatedDist, xEdge_cuda* allChange_cuda, RT_Vertex* SSSP, int* change_d)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < numS)
	{

		if (Edgedone[index] == 1)
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

			 //Check if some other edge was added--mark edge to be added
			if (d_UpdatedDist[node1] > d_UpdatedDist[node2] + edgeWeight)
			{
				Edgedone[index] = 1;
			}

			//Check if correct edge wt was written--mark edge to be added
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

		if (flag != 1 &&  SSSP[index].EDGwt == inf)
		{
			SSSP[index].Dist = inf;
			SSSP[index].Update = true;
			flag = 1;
		}

		//what is p and why the below code needed??
		/*if (d_UpdatedDist[p] == *maxW)
		{
			SSSP->at(i).Dist = *maxW; SSSP->at(i).Update = true; continue;
		}*/


		if (SSSP[index].Dist > d_UpdatedDist[index])
		{
			SSSP[index].Dist = d_UpdatedDist[index];
			SSSP[index].Update = true;
		}

	}
}

__global__ void initializeUpdatedDistOldDist(double* d_UpdatedDist, double* d_OldUpdate, RT_Vertex* SSSP, int X_size)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{
		d_UpdatedDist[index] = SSSP[index].Dist;
		d_OldUpdate[index] = SSSP[index].Dist;
	}
}

__global__ void updateNeighbors(double* d_UpdatedDist, RT_Vertex* SSSP, int X_size, int* d_mychange, int* d_colStartPtr_X, Colwt2* cuda_adjlist_full_X, double inf, int* change_d, int its, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < X_size)
	{

		//If i is updated--update its neighbors as required
		if (SSSP[index].Update)
		{
			d_mychange[index] = 0;
			int px = SSSP[index].Parent;

			//For its nieghbors in X
			for (int j = 0; j < d_colStartPtr_X[index+1] - d_colStartPtr_X[index]; j++)
			{
				////TEPS:
				//*te = *te + 1;
				int myn = cuda_adjlist_full_X[d_colStartPtr_X[index] + j].col;
				double mywt = cuda_adjlist_full_X[d_colStartPtr_X[index] + j].wt;

				if (mywt < 0) { continue; }

				//Check that the Edge weight matches the parent
				 //NOTE:using atomic structures can reduce this step
				if (myn == px)
				{
					if (SSSP[index].EDGwt == inf) { continue; }

					double mydist = d_UpdatedDist[SSSP[index].Parent] + SSSP[index].EDGwt;

					if ((SSSP[index].EDGwt != mywt) || (d_UpdatedDist[index] > mydist))
					{
						SSSP[index].EDGwt = mywt;
						if (d_UpdatedDist[index] >= inf)
						{
							d_UpdatedDist[index] = inf;
						}
						else
						{
							d_UpdatedDist[index] = mydist;
						}

						d_mychange[index] = 1;
						change_d[0] = 1;
					} //end of if

					continue;
				}//end of if


				//Update for Neigbors
				//If distance INF set all neghbors distance to INF at the first iteration
				if (d_UpdatedDist[index] >= inf && its == 0)
				{
					d_UpdatedDist[myn] = inf;
					SSSP[myn].Update = true;
					//mychange[i]=true;
					change_d[0] = 1;

				}
				else {
					//If Distance of myn is larger--make i its parent
					if (d_UpdatedDist[myn] > d_UpdatedDist[index] + mywt)
					{
						//printf("came heruu e\n");

						SSSP[myn].Parent = index;
						SSSP[myn].EDGwt = mywt;
						d_UpdatedDist[myn] = d_UpdatedDist[SSSP[myn].Parent] + SSSP[myn].EDGwt;

						SSSP[myn].Update = true;
						//mychange[i]=true;
						change_d[0] = 1;
					}
				}//end of else
			}//end of for

			 //++++++++++++++++++++//

			//Check Possible Neighbors in R
			for (int j = 0; j < colStartPtr_R[index + 1] - colStartPtr_R[index]; j++)
			{
				int myn = cuda_adjlist_full_R[colStartPtr_R[index] + j].col;
				double mywt = cuda_adjlist_full_R[colStartPtr_R[index] + j].wt;

				 //check if edge is deleted
				if (mywt < 0) { continue; }

				//Check that the Edge weight && Distance matches the parent
				//NOTE:using atomic structur es can reduce this step
				if (myn == px)
				{
					//printf("hhhXX %d %f\n", myn, mywt);
					if (SSSP[index].EDGwt == inf) { continue; }

					double mydist = d_UpdatedDist[SSSP[index].Parent] + SSSP[index].EDGwt;

					if ((SSSP[index].EDGwt != mywt) || (d_UpdatedDist[index] > mydist)) {
						SSSP[index].EDGwt = mywt;
						if (d_UpdatedDist[index] >= inf)
						{
							d_UpdatedDist[index] = inf;
						}
						else
						{
							d_UpdatedDist[index] = mydist;
						}

						d_mychange[index] = 1;
						change_d[0] = 1;
					} //end of if
					continue;
				}//end of if


				//Update for Possible Neighbors
				if (d_UpdatedDist[myn] >= inf) { continue; }

				 //Connect disconnected vertices--after first iter
				if ((d_UpdatedDist[index] >= inf) && (d_UpdatedDist[myn] < inf))
				{

					if (its > 0 && (SSSP[myn].Parent != index))
					{
						SSSP[index].Parent = myn;
						SSSP[index].EDGwt = mywt;
						d_UpdatedDist[index] = d_UpdatedDist[SSSP[index].Parent] + SSSP[index].EDGwt;
					}
					//SSSP->at(myn).Update=true;
					d_mychange[index] = 1;
					change_d[0] = 1;
				 //   continue;
				}

				//If Distance of myn is larger--make i its parent
				if (d_UpdatedDist[myn] > d_UpdatedDist[index] + mywt)
				{
					// printf("came herXXe %d::%d %f %f %f\n", i,myn,UpdatedDist[myn],UpdatedDist[i],mywt);
					SSSP[myn].Parent = index;
					SSSP[myn].EDGwt = mywt;

					if (d_UpdatedDist[myn] >= inf)
					{
						d_UpdatedDist[myn] = inf;
					}
					else
					{
						d_UpdatedDist[myn] = d_UpdatedDist[SSSP[myn].Parent] + SSSP[myn].EDGwt;
					}

					SSSP[myn].Update = true;
					//mychange[i]=true;
					change_d[0] = 1;
				}

			}//end of for

			//if no change occured then update is done
			if (d_mychange[index] != 1)
			{
				SSSP[index].Update = false;
			}
			else
			{
				SSSP[index].Update = true;
			}
		}//end of if Updated
	}//end of for all nodes
}

__global__ void checkIfDistUpdated(int X_size, double* d_OldUpdate, double*  d_UpdatedDist, RT_Vertex* SSSP)
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
			}
		}
	}
}

void edge_update(int* totalChange, int* X_size, int* SSSP_size, xEdge_cuda* allChange_cuda, Colwt2* cuda_adjlist_full_X, int* colStartPtr_X, RT_Vertex* SSSP, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R, int* te, int* nodes);
void rest_update(int* X_size, Colwt2* cuda_adjlist_full_X, int* colStartPtr_X, RT_Vertex* SSSP, Colwt2* cuda_adjlist_full_R, int* colStartPtr_R, int* nodes);


int main() {

	double startx, endx, starty, endy;
	/*double inf = std::numeric_limits<double>::infinity();*/

	/***** Preprocessing to Graph (GUI) ***********/
	int nodes;
	printf("Enter number of total nodes: ");
	scanf("%d", &nodes);

	

	/*** Read Remainder Edges as Graph ***/
	A_Network R;
	for (int i = 0; i < nodes; i++)
	{
		ADJ_Bundle adj_bundle;
		adj_bundle.Row = i;
		R.push_back(adj_bundle);
	}
	string file1 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\fullGraph.txt";
	/*printf("Enter file1(Argv[1]) name: ");
	scanf("%s", &file1);*/
	char* cstr1 = &file1[0];
	readin_graphU(&R, nodes, cstr1);
	

	/*cout << "***R***" << endl;
	for (int i = 0; i < R.size(); i++)
	{
		cout << "row: " << R.at(i).Row << endl;
		for (int j = 0; j < R.at(i).ListW.size(); j++)
		{
			cout << "column: " << R.at(i).ListW.at(j).first << endl;
			cout << "weight: " << R.at(i).ListW.at(j).second << endl;
		}
	}*/

	int* key_R = new int[nodes]; //it stores the node. key is used to find the adj list of a specific node
	int* colStartPtr_R;
	cudaMallocManaged(&colStartPtr_R, (nodes + 1) * sizeof(int)); //we take nodes +1 to store the start ptr of the first row 
	int total_adjmatrix_size_R = 0;
	colStartPtr_R[0] = 0;
	key_R[0] = R.at(0).Row;
	for (int i = 0; i < R.size(); i++)
	{
		key_R[i] = R.at(i).Row;
		int size = R.at(i).ListW.size();
		colStartPtr_R[i + 1] = colStartPtr_R[i] + size; //size of adjacency list per row is stored
		total_adjmatrix_size_R = total_adjmatrix_size_R + size;
	}

	Colwt2* cuda_adjlist_full_R;
	cudaMallocManaged(&cuda_adjlist_full_R, total_adjmatrix_size_R * sizeof(Colwt2));
	for (int i = 0; i < R.size(); i++)
	{
		int colsize = R.at(i).ListW.size();

		for (int j = 0; j < colsize; j++)
		{
			cuda_adjlist_full_R[colStartPtr_R[i] + j].col = R.at(i).ListW.at(j).first;
			cuda_adjlist_full_R[colStartPtr_R[i] + j].wt = R.at(i).ListW.at(j).second;
		}
	}

	/*** Finished Reading CRT Tree **/

	/*** Read the Certificate ***/
	A_Network X;
	for (int i = 0; i < nodes; i++)
	{
		ADJ_Bundle adj_bundle;
		adj_bundle.Row = i;
		X.push_back(adj_bundle);
	}
	string file2 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\SSSP.txt";
	/*printf("Enter file1(Argv[2]) name: ");
	scanf("%s", &file2);*/
	char* cstr2 = &file2[0];
	readin_network(&X, cstr2, -1);
	/*for (int i = 0; i < X.size(); i++)
	{
		cout << "row: " << X.at(i).Row << endl;
		for (int j = 0; j < X.at(i).ListW.size(); j++)
		{
			cout <<"column: "<< X.at(i).ListW.at(j).first << endl;
			cout <<"weight: "<< X.at(i).ListW.at(j).second << endl;
		}
		
	}*/


	int* key_X = new int[nodes]; //it stores the node. key is used to find the adj list of a specific node
	int* colStartPtr_X = new int[nodes + 1]; //we take nodes +1 to store the start ptr of the first row 
	int total_adjmatrix_size_X = 0;
	colStartPtr_X[0] = 0;
	key_X[0] = X.at(0).Row;
	for (int i = 0; i < X.size(); i++)
	{
		key_X[i] = X.at(i).Row;
		int size = X.at(i).ListW.size();
		colStartPtr_X[i + 1] = colStartPtr_X[i] + size; //size of adjacency list per row is stored
		total_adjmatrix_size_X = total_adjmatrix_size_X + size;
	}
	/*for (int i = 0; i < X.size()+1; i++)
	{
		cout << colStartPtr_X[i] << endl;
	}*/
	Colwt2* cuda_adjlist_full_X;
	cudaMallocManaged(&cuda_adjlist_full_X, total_adjmatrix_size_X * sizeof(Colwt2));
	for (int i = 0; i < X.size(); i++)
	{
		int colsize = X.at(i).ListW.size();

		for (int j = 0; j < colsize; j++)
		{
			cuda_adjlist_full_X[colStartPtr_X[i] + j].col = X.at(i).ListW.at(j).first;
			cuda_adjlist_full_X[colStartPtr_X[i] + j].wt = X.at(i).ListW.at(j).second;
			/*cout <<"col: "<< cuda_adjlist_full_X[colStartPtr_X[i] + j].col << " weight: " << cuda_adjlist_full_X[colStartPtr_X[i] + j].wt << endl;*/
		}
	}

	/*** Finished Reading CRT Tree **/

	//List of Changes
	//There will be a list for inserts and a list for delete
	vector<xEdge> allChange;
	allChange.clear();
	 /*** Read set of Changed Edges ***/
	string file3 = "C:\\Users\\khand\\Desktop\\PhD\\CUDA test\\Test\\test 1\\changeEdges.txt";
	/*printf("Enter file1(Argv[3]) name: ");
	scanf("%s", &file3);*/
	char* cstr3 = &file3[0];
	readin_changes(cstr3, &allChange);

	/*for (int i = 0; i < allChange.size(); i++)
	{
		cout <<"inst: "<< allChange.at(i).inst << endl;
		cout << "node1: " << allChange.at(i).theEdge.node1 << " node2: " << allChange.at(i).theEdge.node2 << " weight: " << allChange.at(i).theEdge.edge_wt << endl;
	}*/

	//new addition
	xEdge_cuda* allChange_cuda;
	int totalChange = allChange.size();
	cudaMallocManaged(&allChange_cuda, totalChange * sizeof(xEdge_cuda));
	for (int i = 0; i < totalChange; i++)
	{
		allChange_cuda[i].node1 = allChange.at(i).theEdge.node1;
		allChange_cuda[i].node2 = allChange.at(i).theEdge.node2;
		allChange_cuda[i].edge_wt = allChange.at(i).theEdge.edge_wt;
		allChange_cuda[i].inst = allChange.at(i).inst;
		allChange_cuda[i].insertedToDatastructure = allChange.at(i).insertedToDatastructure;
	}
	/*** Finished Reading Changed Edges **/

	//Initializing  Rooted Tree
	RT_Vertex* SSSP;
	cudaMallocManaged(&SSSP, nodes * sizeof(RT_Vertex));
	int* stencil; //stencil is used for tracking which node is being affected. 
	/*cudaMallocManaged(&stencil, nodes * sizeof(int));*/
	cudaMalloc((void**)&stencil, nodes * sizeof(int));
	int* stencil_c = new int[nodes];
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
		cudaMalloc((void**)&d_affectedPointer, nodes * sizeof(int));
		affectedPointer = (int*)calloc(nodes, sizeof(int));
		cudaMemcpy(d_affectedPointer, affectedPointer, nodes * sizeof(int), cudaMemcpyHostToDevice);
		/*cudaMallocManaged(&affectedPointer, nodes * sizeof(int));*/

		totalAffectedNode = 1;
		int start = 0, end = 0;
		int* affected_nodes;
		affected_nodes = (int*)calloc(totalAffectedNode, sizeof(int));
		affected_nodes[0] = src;
		cudaStream_t stream1;
		cudaError_t result;
		result = cudaStreamCreate(&stream1);
		while (totalAffectedNode > 0)
		{
			vector<int> affectedNodeAlias;
			/*ThreadHelper* threadHelpers;
			cudaMallocManaged(&threadHelpers, (totalAffectedNode+1) * sizeof(ThreadHelper));
			int offset = 0;
			int totalNumberofCudaThread = 0;*/

			for (int i = 0; i < totalAffectedNode; i++)
			{
				affectedNodeAlias.push_back(affected_nodes[i]);
			}
			for (int i = 0; i < totalAffectedNode; i++)
			{
				p = affectedNodeAlias.at(i);
				src = p;
				/*cout << "src: " << src << endl;*/
				start = colStartPtr_X[p];
				end = colStartPtr_X[p+1];
				int numberofCudaThread = end - start;
				/*threadHelpers[i].offset = offset;
				offset = offset + numberofCudaThread;*/

				/*cout << "adj node 4:" << cuda_adjlist_full_X[4].col<<endl;
				for (int i = 0; i < numberofCudaThread; i++)
				{
					cout << "adj node ptr" << i + start << endl;
					cout << "adj node"<<cuda_adjlist_full_X[i + start].col<<endl;
					cout << "adj node"<<cuda_adjlist_full_X[i + start].col<<endl;
				}*/
				create_tree << <(numberofCudaThread / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK, 0, stream1 >> > (cuda_adjlist_full_X, start, SSSP, src, d_affectedPointer, numberofCudaThread);
				
				/*threadHelpers[i].src = p;
				threadHelpers[i].start = start;*/
			}
			/*threadHelpers[totalAffectedNode].offset = offset;
			totalNumberofCudaThread = offset;
			create_tree2<<<(offset / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK>> > (cuda_adjlist_full_X, SSSP, d_affectedPointer, threadHelpers, totalNumberofCudaThread);
			cudaDeviceSynchronize();*/
			cudaStreamSynchronize(stream1);
			thrust::device_ptr<int> affectedPointer_alias(d_affectedPointer); // converting raw ptr to device_ptr
			cudaMemcpy(affectedPointer, d_affectedPointer, nodes * sizeof(int), cudaMemcpyDeviceToHost);
			for (int i = 0; i < nodes; i++)
			{
				cout << "after kernel call:" << affectedPointer[i] << endl;
			}
			/*thrust::device_ptr<int> affectedPointer_alias(affectedPointer);*/
			thrust::device_vector<int> affectedPointer_vector(affectedPointer_alias, affectedPointer_alias + nodes); //converting device_ptr to device_vector
			totalAffectedNode = thrust::count(affectedPointer_vector.begin(), affectedPointer_vector.end(), 1); //count the number of affected node
			cout << "totalAffectedNode: " << totalAffectedNode<<endl;
			affected_nodes = (int*)realloc(affected_nodes, totalAffectedNode * sizeof(int));
			/*affectedPointer = thrust::raw_pointer_cast(&affectedPointer_vector[0]);*/
			thrust::copy_if(thrust::host, stencil_c, stencil_c + nodes, affectedPointer, affected_nodes, is_affected());
			cout << "affected nodes " << endl;
			for (int i = 0; i < totalAffectedNode; i++)
			{
				cout << affected_nodes[i] << endl;
			}
			
			affectedPointer = (int*)calloc(nodes, sizeof(int));
			cudaMemcpy(d_affectedPointer, affectedPointer, nodes * sizeof(int), cudaMemcpyHostToDevice);
			/*cudaFree(threadHelpers);*/
		}
		free(affectedPointer);
		free(affected_nodes);
		cudaFree(d_affectedPointer);

		//Time calculation
		auto stopTime = high_resolution_clock::now();
		// Time calculation
		auto duration = duration_cast<microseconds>(stopTime - startTime);
		cout << "Time taken by create-tree function: "
			<< duration.count() << " microseconds" << endl;

		//test
		cout << "sssp" << endl;
		for (int i = 0; i < nodes; i++)
		{
			cout << "node" << i << endl;
			cout << "dist" << SSSP[i].Dist << endl;
			cout << "parent" << SSSP[i].Parent << endl;
		}
		//test end
		//edge_update function
		//Update the inserted and delted edges in the tree
		int x_size = X.size();
		int SSSP_size = nodes; //considering all nodes are participating in the SSSP
		int te = 0;
		edge_update(&totalChange, &x_size, &SSSP_size, allChange_cuda, cuda_adjlist_full_X, colStartPtr_X, SSSP, cuda_adjlist_full_R, colStartPtr_R, &te, &nodes);
		cout << "after edge_update fn" << endl;
		rest_update(&x_size, cuda_adjlist_full_X, colStartPtr_X, SSSP, cuda_adjlist_full_R, colStartPtr_R, &nodes);
		cout << "after rest_update fn" << endl;
	}
	else
	{
		//****below code needs modification
		/*SCC.clear();
		readin_SCC(argv[2], &SCC);
		update_SCC(&X, &SCC, &allChange);
		print_network(X);*/
	}
	cout << "SSSP" << endl;
	for (int i = 0; i < nodes; i++)
	{
		cout << "*******" << endl;
		cout << "node"<<i<<endl<<"dist" << SSSP[i].Dist << endl<< "parent" << SSSP[i].Parent << endl;
	}


	cudaFree(colStartPtr_R);
	cudaFree(cuda_adjlist_full_R);
	cudaFree(colStartPtr_X);
	cudaFree(cuda_adjlist_full_X);
	cudaFree(allChange_cuda);
	cudaFree(SSSP);
	cudaFree(stencil);

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
	cudaMemcpy(UpdatedDist, d_UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyDeviceToHost); //not required


	
	int numS = *totalChange;
	int* d_colStartPtr_X;
	cudaMalloc((void**)&d_colStartPtr_X, (*nodes + 1) * sizeof(int));
	cudaMemcpy(d_colStartPtr_X, colStartPtr_X, (*nodes + 1) * sizeof(int), cudaMemcpyHostToDevice);

	insertDeleteEdge << < (numS / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (allChange_cuda, Edgedone, SSSP, numS, *X_size, d_colStartPtr_X, cuda_adjlist_full_X, d_UpdatedDist, inf, cuda_adjlist_full_R, colStartPtr_R);
	cudaDeviceSynchronize();
	int* Edgedone_c = new int[*totalChange];
	cudaMemcpy(Edgedone_c, Edgedone, *totalChange * sizeof(int), cudaMemcpyDeviceToHost); //not req.
	cudaMemcpy(UpdatedDist, d_UpdatedDist, (*X_size) * sizeof(double), cudaMemcpyDeviceToHost); //not req.

	//Test start
	for (int i = 0; i < *totalChange; i++)
	{
		cout << Edgedone_c[i] << endl;

	}
	//Test end


	//Go over the inserted edges to see if they need to be changed
	int* change_d = new int[1];
	int* change = new int[1];
	change[0] = 1;
	cudaMalloc((void**)&change_d, 1 * sizeof(int));
	cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);
	while (change[0] == 1)
	{
		change[0] = 0;
		cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice);
		checkInsertedEdges << < (numS / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (numS, Edgedone, d_UpdatedDist, allChange_cuda, SSSP, change_d);
		cudaDeviceSynchronize();
		cudaMemcpy(change, change_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);
		cout << "change"<< change[0]<<endl;
	}
	cout << "check1";
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
	cudaMemcpy(change_d, change, 1 * sizeof(int), cudaMemcpyHostToDevice); 

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


	int* mychange;
	//Store the old updated distance value
	mychange = (int*)calloc(*X_size, sizeof(int));
	int* d_mychange;
	cudaMalloc((void**)&d_mychange, (*X_size) * sizeof(int));
	cudaMemcpy(d_mychange, mychange, (*X_size) * sizeof(int), cudaMemcpyHostToDevice);

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
		updateNeighbors << <((*X_size) / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >> > (d_UpdatedDist, SSSP, *X_size, d_mychange, d_colStartPtr_X, cuda_adjlist_full_X, inf, change_d, its, cuda_adjlist_full_R, colStartPtr_R);
		cudaDeviceSynchronize();
		cudaMemcpy(change, change_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);

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
	free(mychange);
	cudaFree(change_d);
	cudaFree(d_UpdatedDist);
	cudaFree(d_OldUpdate);
	cudaFree(d_mychange);
	cudaFree(d_colStartPtr_X);

	return;
}
