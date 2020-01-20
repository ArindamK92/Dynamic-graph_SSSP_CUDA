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
	bool insertedToDatastructure;

	xEdge_cuda()
	{
		insertedToDatastructure = false;
	}
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
	int Level; //Distance from the root
	int marked; //whether the vertex and the edge connecting its parent ..
				//..exists(-1); has been deleted(-2); is marked for replacement (+ve value index of changed edge)
	/*vector<int> PossN;*/ //possible neighbor to check for /**can't use vector in device. need workaround if PossN is necessary

	double Dist;  //Distance from root
	bool Update;  //Whether the distance of this edge was updated
	int degree;
	////Constructor
	//RT_Vertex()
	//{
	//	Root = -1;
	//	Parent = -1;
	//	EDGwt = 0.0;
	//	Level = -1;
	//	marked = -1;

	//	/*PossN.clear();*/
	//	Dist = 0.0;
	//	Update = false;
	//}
	////Destructor
	//void clear()
	//{}
};
//The Rooted tree is a vector of structure RT_Vertex;

////functions////
//Assumes the all nodes present
//Node starts from 0
//Total number of vertices=nodes and are consecutively arranged
//reads only the edges in the edge list, does not reverse them to make undirected
void readin_graphU(A_Network* X, int nodes, char* myfile)
{

	//Initialize for X
	ADJ_Bundle AList; //create current Adj_Bundle
	X->resize(nodes, AList);
	int_double colvals;
	int node1;
	double edge_wt;


	/////my fn
	string fileName = myfile;
	string line;
	ifstream fin;

	// by default open mode = ios::in mode 
	fin.open(fileName);

	// Execute a loop until EOF (End of File) 
	/*vector<float> inputList;*/
	vector<int_double> ListW;
	int current_node = 0;
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
		}
	}
	/*for (float a : inputList)
	{
		cout << a << endl;
	}*/
	for (int i = 0; i < inputList.size() - 3; i++)
	{
		if (i % 3 == 0)
		{
			node1 = (int)inputList.at(i);

		}
		else if (i % 3 == 1)
		{
			colvals.first = (int)inputList.at(i);

		}
		else {
			colvals.second = inputList.at(i);
			X->at(node1).ListW.push_back(colvals);
		}


		//if (current_node != node1)
		//{
		//	{
		//		X->at(current_node).ListW = ListW;
		//		ListW.clear(); //Clear List of Edges;

		//	} //end of if

		//	//Set to next node
		//	current_node = node1;
		//}

		//ListW.push_back(colvals);
	}
}

/**** Reading File to Create Network******/
template <class myNetworkType>
void readin_network(myNetworkType* X, char* file, int XtraN)
{
	int nodes;

	vector<Edge> a;//list of all edges of network
	Edge myedge;
	string fileName = file;
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
		}
	}
	for (int i = 0; i < inputList.size() - 3; i++)
	{

		if (i % 3 == 0)
		{
			myedge.node1 = (int)inputList.at(i);

		}
		else if (i % 3 == 1)
		{
			myedge.node2 = (int)inputList.at(i);

		}
		else {
			myedge.edge_wt = inputList.at(i);
			a.push_back(myedge);

		}
	}

	/*for (int i = 0; i < a.size(); i++)
	{
		cout << "node 1: " << a.at(i).node1 << "node 2: " << a.at(i).node2 << "weight: " << a.at(i).edge_wt << endl;
	}*/
	create_Network(&a, 0, X, XtraN);

	return;
}
void create_Network(vector<Edge>* b, int bf_size, A_Network* X, int Ns)
{
	for (int i = 0; i < b->size(); i++)
	{
		int row = b->at(i).node1;
		int_double colWt;
		colWt.first = b->at(i).node2;
		colWt.second = b->at(i).edge_wt;

		X->at(row).ListW.push_back(colWt);
	}
}



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
		}
	}

	for (int i = 0; i < inputList.size() - 4; i++)
	{
		if (i % 4 == 0)
		{
			myedge.theEdge.node1 = (int)inputList.at(i);
		}
		else if (i % 4 == 1)
		{
			myedge.theEdge.node2 = (int)inputList.at(i);
		}
		else if (i % 4 == 2) {
			myedge.theEdge.edge_wt = inputList.at(i);
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
		}
	}
	return;
}//end of function
