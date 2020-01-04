//#include <stdio.h>
//#include <iostream>
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//#include<vector>
//#include <set>
//#include<queue>
//#include <chrono> 
//#include <thread>
//#include<queue>
//#include<stack>
//#include<list>
//#include "all_structures.h"
//#include <thrust/find.h>
//#include <thrust/device_vector.h>
//#include <thrust/count.h>
//#include <thrust/copy.h>
//#include <thrust/execution_policy.h>
//#define THREADS_PER_BLOCK 512 //we can change it
//
//using namespace std;
//using namespace std::chrono;
//
//__global__ void initialize(int nodes, int src, int* colStartPtrRow, RT_Vertex* SSSP, int* stencil)
//{
//	int index = threadIdx.x + blockIdx.x * blockDim.x;
//	if (index < nodes)
//	{
//		if (index == src) { SSSP[index].Root = -1; }
//		else { SSSP[index].Root = index; }
//		SSSP[index].Level = 0;
//		SSSP[index].Dist = 0.0;
//		SSSP[index].degree = colStartPtrRow[index + 1] - colStartPtrRow[index];
//		stencil[index] = index; 
//	}
//}
//
//__global__ void create_tree(Colwt2* cuda_adjlist_full_X, int start, int stop, RT_Vertex* SSSP, int src, int* affectedPointer, int numberofCudaThread)
//{
//	int index = threadIdx.x + blockIdx.x * blockDim.x;
//	if (index < numberofCudaThread)
//	{
//		int y = cuda_adjlist_full_X[index + start].col;
//		double mywt = cuda_adjlist_full_X[index + start].wt;
//
//		SSSP[y].Parent = src; //mark the parent
//		SSSP[y].EDGwt = mywt; //mark the edgewt
//		SSSP[y].Level = 1; //mark the Level
//		SSSP[y].Root = -1;
//		SSSP[y].Dist = mywt;
//		affectedPointer[y] = 1;
//	}
//
//}
//
//struct is_affected
//{
//	__host__ __device__
//		bool operator()(const int x)
//	{
//		return (x  == 1);
//	}
//};
//
//int main() {
//	
//	double startx, endx, starty, endy;
//	/***** Preprocessing to Graph (GUI) ***********/
//	int nodes;
//	scanf("Enter number of total nodes: %d", &nodes);
//
//	//List of Changes
//	//There will be a list for inserts and a list for delete
//	vector<xEdge> allChange;
//	allChange.clear();
//
//	/*** Read Remainder Edges as Graph ***/
//	A_Network R;
//	string file1;
//	scanf("Enter file1(Argv[1]) name: %s", &file1);
//	char* cstr1 = &file1[0];
//	readin_graphU(&R, nodes, cstr1);
//
//	//**new addition**/
//	int* key_R = new int[nodes]; //it stores the node. key is used to find the adj list of a specific node
//	int* colStartPtr_R;
//	cudaMallocManaged(&colStartPtr_R, (nodes + 1) * sizeof(int)); //we take nodes +1 to store the start ptr of the first row 
//	int total_adjmatrix_size_R = 0;
//	colStartPtr_R[0] = 0;
//	key_R[0] = R.at(0).Row;
//	for (int i = 0; i < R.size(); i++)
//	{
//		key_R[i] = R.at(i).Row;
//		int size = R.at(i).ListW.size();
//		colStartPtr_R[i + 1] = colStartPtr_R[i] + size; //size of adjacency list per row is stored
//		total_adjmatrix_size_R = total_adjmatrix_size_R + size;
//	}
//
//	Colwt2* cuda_adjlist_full_R;
//	cudaMallocManaged(&cuda_adjlist_full_R, total_adjmatrix_size_R * sizeof(Colwt2));
//	for (int i = 0; i < R.size(); i++)
//	{
//		int colsize = R.at(i).ListW.size();
//
//		for (int j = 0; j < colsize; j++)
//		{
//			cuda_adjlist_full_R[colStartPtr_R[i] + j].col = R.at(i).ListW.at(j).first;
//			cuda_adjlist_full_R[colStartPtr_R[i] + j].wt = R.at(i).ListW.at(j).second;
//		}
//	}
//
//	/*** Finished Reading CRT Tree **/
//
//	/*** Read the Certificate ***/
//	A_Network X;
//	string file2;
//	scanf("Enter file1(Argv[2]) name: %s", &file2);
//	char* cstr2 = &file2[0];
//	readin_network(&X, cstr2, -1);
//
//
//	
//	//**new addition**//
//	int* key_X = new int[nodes]; //it stores the node. key is used to find the adj list of a specific node
//	int* colStartPtr_X; 
//	cudaMallocManaged(&colStartPtr_X, (nodes + 1) * sizeof(int)); //we take nodes +1 to store the start ptr of the first row 
//	int total_adjmatrix_size_X = 0;
//	colStartPtr_X[0] = 0;
//	key_X[0] = X.at(0).Row;
//	for (int i = 0; i < X.size(); i++)
//	{
//		key_X[i] = X.at(i).Row;
//		int size = X.at(i).ListW.size();
//		colStartPtr_X[i + 1] = colStartPtr_X[i] + size; //size of adjacency list per row is stored
//		total_adjmatrix_size_X = total_adjmatrix_size_X + size;
//	}
//
//	Colwt2* cuda_adjlist_full_X;
//	cudaMallocManaged(&cuda_adjlist_full_X, total_adjmatrix_size_X * sizeof(Colwt2));
//	for (int i = 0; i < X.size(); i++)
//	{
//		int colsize = X.at(i).ListW.size();
//
//		for (int j = 0; j < colsize; j++)
//		{
//			cuda_adjlist_full_X[colStartPtr_X[i] + j].col = X.at(i).ListW.at(j).first;
//			cuda_adjlist_full_X[colStartPtr_X[i] + j].wt = X.at(i).ListW.at(j).second;
//		}
//	}
//	/*** Finished Reading CRT Tree **/
//
//	 /*** Read set of Changed Edges ***/
//	string file3;
//	scanf("Enter file1(Argv[3]) name: %s", &file3);
//	char* cstr3 = &file3[0];
//	readin_changes(cstr3, &allChange);
//	
//	//new addition
//	xEdge_cuda* allChange_cuda;
//	int totalChange = allChange.size();
//	cudaMallocManaged(&allChange_cuda, totalChange * sizeof(xEdge_cuda));
//	for (int i = 0; i < totalChange; i++)
//	{
//		allChange_cuda[i].node1 = allChange.at(i).theEdge.node1;
//		allChange_cuda[i].node2 = allChange.at(i).theEdge.node2;
//		allChange_cuda[i].edge_wt = allChange.at(i).theEdge.edge_wt;
//		allChange_cuda[i].inst = allChange.at(i).inst;
//		allChange_cuda[i].insertedToDatastructure = allChange.at(i).insertedToDatastructure;
//	}
//	/*** Finished Reading Changed Edges **/
//
//	//Initializing  Rooted Tree
//	RT_Vertex* SSSP;
//	cudaMallocManaged(&SSSP, nodes * sizeof(RT_Vertex));
//	int* stencil; //stencil is used for tracking which node is being affected. 
//	cudaMallocManaged(&stencil, nodes * sizeof(int));
//	/*vector<SCC_vertex>SCC;*/
//	int graphDirectedUndirectedIndicator = 0; // Should be 1 for SCC, 0 for not SCC. need to modify if we want SCC
//	
//	int source;
//	scanf("Enter source node: %d", &source);
//
//	if (graphDirectedUndirectedIndicator == 0) {
//		int src = source; //the source from which the paths are computed
//		initialize <<<(nodes / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK >>> (nodes, src, colStartPtr_X, SSSP, stencil); //kernet call
//
//		//Code for create_tree:
//		int totalAffectedNode; //alias of numberOfAffectedNode
//
//		int* affectedPointer;
//		cudaMallocManaged(&affectedPointer, nodes * sizeof(int));
//
//		totalAffectedNode = 1;
//		int start = 0, end = 0;
//		int* affected_nodes;
//		affected_nodes = (int*)calloc(totalAffectedNode, sizeof(int));
//		affected_nodes[0] = src;
//		cudaStream_t stream1;
//		cudaError_t result;
//		result = cudaStreamCreate(&stream1);
//		while(totalAffectedNode>0)
//		{ 
//			for (int i = 0; i < totalAffectedNode; i++)
//			{
//				start = colStartPtr_X[affected_nodes[i]];
//				end = colStartPtr_X[affected_nodes[i] + 1];
//				int numberofCudaThread = end - start;
//				create_tree <<<(numberofCudaThread / THREADS_PER_BLOCK) + 1, THREADS_PER_BLOCK, 0, stream1 >>> (cuda_adjlist_full_X, start, stop, SSSP, src, affectedPointer, numberofCudaThread);
//			}
//			/*for (auto& t : exe_threads) t.join();*/
//			cudaStreamSynchronize(stream1);
//			thrust::device_ptr<int> affectedPointer_alias(affectedPointer); // converting raw ptr to device_ptr
//			thrust::device_vector<int> affectedPointer_vector(affectedPointer_alias, affectedPointer_alias + nodes); //converting device_ptr to device_vector
//			totalAffectedNode = thrust::count(affectedPointer_vector.begin(), affectedPointer_vector.end(), 1); //count the number of affected node
//			affected_nodes = (int*)realloc(affected_nodes, totalAffectedNode * sizeof(int));
//			thrust::copy_if(thrust::host, affectedPointer, affectedPointer + nodes, stencil, affected_nodes, is_affected());
//
//		}
//
//		free(affected_nodes);
//
//
//		//****below code needs modification
//		//create_treeS(&X, &R, &SSSP, src, p);
//
//
//		//double maxV = (double)maxW * X.size();
//
//		////Update the inserted and delted edges in the tree
//		//int te = 0;
//		//edge_update(&allChange, &X, &SSSP, &R, &maxV, &te, p);
//
//		//rest_update(&X, &SSSP, &R, &maxV, &te, p);
//		cudaFree(affectedPointer);
//
//	}
//	else
//	{
//		//****below code needs modification
//		/*SCC.clear();
//		readin_SCC(argv[2], &SCC);
//		update_SCC(&X, &SCC, &allChange);
//		print_network(X);*/
//	}
//
//	cudaFree(colStartPtr_R);
//	cudaFree(cuda_adjlist_full_R);
//	cudaFree(colStartPtr_X);
//	cudaFree(cuda_adjlist_full_X);
//	cudaFree(allChange_cuda);
//	cudaFree(SSSP);
//	cudaFree(stencil);
//
//	return 0;
//}
//
