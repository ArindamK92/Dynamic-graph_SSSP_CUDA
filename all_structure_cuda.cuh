#include <stdio.h>
#include <iostream>
//#include<list>
#include<vector> 
#include <fstream> 
#include <sstream>
using namespace std;

/******* Network Structures *********/
struct Colwt2 {
	int col;
	double wt;
};
//Note: Edges are not ordered, unless specified by code
// Node+Weight = -1;0 indicates buffer space
//Structure for Edge
struct Edge
{
	int node1;
	int node2;
	double edge_wt;
};

Edge create(int n1, int n2, double wt)
{
	Edge e;
	e.node1 = n1;
	e.node2 = n2;
	e.edge_wt = wt;

	return e;
}
//========================|
//Structure to indicate whether Edge is to be inserted/deleted
struct xEdge {
	Edge theEdge;
	int inst;
	bool insertedToDatastructure;

	xEdge()
	{
		insertedToDatastructure = false;
	}
	void clear()
	{}
};

struct xEdge_cuda {
	//Edge related parameters
	int node1;
	int node2;
	double edge_wt;
	//End of Edge related parameters
	int inst;

	void clear()
	{}
};
struct ThreadHelper
{
	int src; //source
	int start; //start points to the starting of adjlist of the node in the full adj list
	int offset; //it stores the lenght of adjlist upto last node
};
/*** Pairs ***/
typedef pair<int, int> int_int;  /** /typedef pair of integers */
typedef pair<int, double> int_double; /** /typedef pair of integer and double */
typedef pair<double, int> double_int; /** /typedef pair of integer and double */

//Structure in STATIC Adjacency List---For diagram go to () 
//Rows=global ID of the rows
//For edges connected with Rows
//NListW.first=Column number
//NListW.second=Value of edge
struct ADJ_Bundle
{
	int Row;
	vector <int_double> ListW;

	//Constructor
	ADJ_Bundle() { ListW.resize(0); }

	//Destructor
	void clear()
	{
		while (!ListW.empty()) { ListW.pop_back(); }
	}


};
typedef  vector<ADJ_Bundle> A_Network;



// Data Structure for each vertex in the rooted tree
struct RT_Vertex
{
	int Root;    //root fo the tree
	int Parent; //mark the parent in the tree
	double EDGwt; //mark weight of the edge
	int marked; //whether the vertex and the edge connecting its parent ..
				//..exists(-1); has been deleted(-2); is marked for replacement (+ve value index of changed edge)

	double Dist;  //Distance from root
	bool Update;  //Whether the distance of this edge was updated
};
//The Rooted tree is a vector of structure RT_Vertex;

////functions////
//Assumes the all nodes present
//Node starts from 0
//Total number of vertices=nodes and are consecutively arranged
//reads only the edges in the edge list, does not reverse them to make undirected



//modification needed
void readin_changes(char* myfile, vector<xEdge>* allChange)
{
	xEdge myedge;
	int ins;
	string fileName = myfile;
	string line;
	ifstream fin;

	// by default open mode = ios::in mode 
	fin.open(fileName);
	int i = 0;
	vector<float> inputList;
	//Read line 
	while (fin) {

		// Read a Line from File 
		getline(fin, line);
		istringstream iss(line);
		std::string word;

		while (iss >> word)
		{
			inputList.push_back(stof(word));
			if (i % 4 == 0)
			{
				myedge.theEdge.node1 = (int)inputList.at(i);
				i++;
			}
			else if (i % 4 == 1)
			{
				myedge.theEdge.node2 = (int)inputList.at(i);
				i++;
			}
			else if (i % 4 == 2) {
				myedge.theEdge.edge_wt = inputList.at(i);
				i++;
			}
			else {
				ins = (int)inputList.at(i);
				if (ins == 1)
				{
					myedge.inst = true;
				}
				else
				{
					myedge.inst = false;
				}
				allChange->push_back(myedge);
				i = 0;
				inputList.clear();
			}
		}
	}


	return;
}//end of function


void readin_graphU2(A_Network* X, int nodes, char* myfile)
{
	FILE* graph_file;
	char line[128];

	graph_file = fopen(myfile, "r");
	int l = 0;
	int prev_node = 0;
	int_double dummy;
	dummy.first = 1;
	dummy.second = 0;
	vector <int_double> ListW;
	ADJ_Bundle adjnode;

	while (fgets(line, 128, graph_file) != NULL)
	{
		int n1, n2;
		int wt;
		//Read line
		sscanf(line, "%d %d %d", &n1, &n2, &wt);
		//cout <<n1<<"\n";
		//Number of nodes given in the first line
//        if(l==0)
//        {l++; continue;}
		if (prev_node != n1)
		{
			adjnode.Row = prev_node;
			adjnode.ListW = ListW;
			X->push_back(adjnode);
			/*X->at(prev_node).ListW = ListW;*/
			ListW.clear();
			adjnode.clear();
		}
		prev_node = n1;
		dummy.first = n2;
		dummy.second = (double)wt;
		ListW.push_back(dummy);
	}//end of while
	adjnode.Row = prev_node;
	adjnode.ListW = ListW;
	X->push_back(adjnode);
	fclose(graph_file);

	return;
}

void readin_graphU4(int* colStartPtr, Colwt2* cuda_adjlist_full, char* myfile, int* nodes)
{
	FILE* graph_file;
	char line[128];

	graph_file = fopen(myfile, "r");
	int l = 0;
	int prev_node = 0;
	int_double dummy;
	dummy.first = 1;
	dummy.second = 0;
	vector <int_double> ListW;
	ADJ_Bundle adjnode;
	colStartPtr[0] = 0;
	while (fgets(line, 128, graph_file) != NULL)
	{
		int n1, n2, wt;
		//Read line
		sscanf(line, "%d %d %d", &n1, &n2, &wt);
		if (prev_node != n1)
		{
			colStartPtr[n1] = l;
		}
		prev_node = n1;
		cuda_adjlist_full[l].col = n2;
		cuda_adjlist_full[l].wt = (double)wt;
		l++;

	}//end of while
	colStartPtr[*nodes] = l; //last element stores dummy last value to indicate the end of the ptr
	fclose(graph_file);

	return;
}