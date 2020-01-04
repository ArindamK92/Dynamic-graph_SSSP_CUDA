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

__global__ void create_tree(Colwt2* cuda_adjlist_full_X, int start, int stop, RT_Vertex* SSSP, int src, int* d_affectedPointer, int numberofCudaThread)
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
		SSSP[y].Dist = mywt;
		d_affectedPointer[y] = 1;
		printf("end if***");
		
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

int main() {

	double startx, endx, starty, endy;
	/***** Preprocessing to Graph (GUI) ***********/
	int nodes;
	printf("Enter number of total nodes: ");
	scanf("%d", &nodes);

	//List of Changes
	//There will be a list for inserts and a list for delete
	vector<xEdge> allChange;
	allChange.clear();

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
				/*cout << "adj node 4:" << cuda_adjlist_full_X[4].col<<endl;
				for (int i = 0; i < numberofCudaThread; i++)
				{
					cout << "adj node ptr" << i + start << endl;
					cout << "adj node"<<cuda_adjlist_full_X[i + start].col<<endl;
					cout << "adj node"<<cuda_adjlist_full_X[i + start].col<<endl;
				}*/
				create_tree << <(numberofCudaThread / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK, 0, stream1 >> > (cuda_adjlist_full_X, start, end, SSSP, src, d_affectedPointer, numberofCudaThread);
			}
			/*for (auto& t : exe_threads) t.join();*/
			cudaStreamSynchronize(stream1);
			thrust::device_ptr<int> affectedPointer_alias(d_affectedPointer); // converting raw ptr to device_ptr
			cudaMemcpy(affectedPointer, d_affectedPointer, nodes * sizeof(int), cudaMemcpyDeviceToHost);
			/*for (int i = 0; i < nodes; i++)
			{
				cout << "after kernel call:" << affectedPointer[i] << endl;
			}*/
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
			free(affectedPointer);
			affectedPointer = (int*)calloc(nodes, sizeof(int));
			cudaMemcpy(d_affectedPointer, affectedPointer, nodes * sizeof(int), cudaMemcpyHostToDevice);
			
		}

		free(affected_nodes);


		//****below code needs modification
		//create_treeS(&X, &R, &SSSP, src, p);


		//double maxV = (double)maxW * X.size();

		////Update the inserted and delted edges in the tree
		//int te = 0;
		//edge_update(&allChange, &X, &SSSP, &R, &maxV, &te, p);

		//rest_update(&X, &SSSP, &R, &maxV, &te, p);
		cudaFree(d_affectedPointer);

	}
	else
	{
		//****below code needs modification
		/*SCC.clear();
		readin_SCC(argv[2], &SCC);
		update_SCC(&X, &SCC, &allChange);
		print_network(X);*/
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

