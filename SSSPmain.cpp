#include <stdio.h>
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<vector>
#include <set>
#include <unordered_map>
#include<queue>
#include <chrono> 
#include <thread>
#include<queue>
#include<stack>
#include<list>
#include "GraphStructure.h"
#include "readin_data.h"
#define THREADS_PER_BLOCK 512 //we can change it

using namespace std;
using namespace std::chrono;

__global__ void initialize(int nodes, int src, int* colStartPtrRow, RT_Vertex* SSSP)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if (index < nodes)
	{
		if (index == src) { SSSP[index].Root = -1; }
		else { SSSP[index].Root = index; }
		SSSP[index].Level = 0;
		SSSP[index].Dist = 0.0;
		SSSP[index].degree = colStartPtrRow[index + 1] - colStartPtrRow[index];
	}
}

int main() {
	
	double startx, endx, starty, endy;
	/***** Preprocessing to Graph (GUI) ***********/
	int nodes;
	scanf("Enter number of total nodes: %d", &nodes);

	//List of Changes
	//There will be a list for inserts and a list for delete
	vector<xEdge> allChange;
	allChange.clear();

	/*** Read Remainder Edges as Graph ***/
	A_Network R;
	string file1;
	scanf("Enter file1(Argv[1]) name: %s", &file1);
	char* cstr1 = &file1[0];
	readin_graphU(&R, nodes, cstr1);

	//**new addition**/
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
	string file2;
	scanf("Enter file1(Argv[2]) name: %s", &file2);
	char* cstr2 = &file2[0];
	readin_network(&X, cstr2, -1);


	
	//**new addition**//
	int* key_X = new int[nodes]; //it stores the node. key is used to find the adj list of a specific node
	int* colStartPtr_X; 
	cudaMallocManaged(&colStartPtr_X, (nodes + 1) * sizeof(int)); //we take nodes +1 to store the start ptr of the first row 
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

	Colwt2* cuda_adjlist_full_X;
	cudaMallocManaged(&cuda_adjlist_full_X, total_adjmatrix_size_X * sizeof(Colwt2));
	for (int i = 0; i < X.size(); i++)
	{
		int colsize = X.at(i).ListW.size();

		for (int j = 0; j < colsize; j++)
		{
			cuda_adjlist_full_X[colStartPtr_X[i] + j].col = X.at(i).ListW.at(j).first;
			cuda_adjlist_full_X[colStartPtr_X[i] + j].wt = X.at(i).ListW.at(j).second;
		}
	}
	/*** Finished Reading CRT Tree **/

	 /*** Read set of Changed Edges ***/
	string file3;
	scanf("Enter file1(Argv[3]) name: %s", &file3);
	char* cstr3 = &file3[0];
	readin_changes(cstr3, &allChange);
	
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
	vector<SCC_vertex>SCC;
	int graphDirectedUndirectedIndicator = 0; // Should be 1 for SCC, 0 for not SCC. need to modify if we want SCC
	
	int source;
	scanf("Enter source node: %d", &source);

	if (graphDirectedUndirectedIndicator == 0) {
		int src = source; //the source from which the paths are computed
		initialize << <nodes / THREADS_PER_BLOCK, THREADS_PER_BLOCK >> > (nodes, src, colStartPtr_X, SSSP); //kernet call


		//****below code needs modification
		//create_treeS(&X, &R, &SSSP, src, p);


		//double maxV = (double)maxW * X.size();

		////Update the inserted and delted edges in the tree
		//int te = 0;
		//edge_update(&allChange, &X, &SSSP, &R, &maxV, &te, p);

		//rest_update(&X, &SSSP, &R, &maxV, &te, p);


	}
	else
	{
		//****below code needs modification
		/*SCC.clear();
		readin_SCC(argv[2], &SCC);
		update_SCC(&X, &SCC, &allChange);
		print_network(X);*/
	}

	return 0;
}

