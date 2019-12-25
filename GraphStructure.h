#include <stdio.h>
#include <iostream>
#include "structure_defs.h"


//_______________________
struct Colwt2 {
	int col;
	double wt;
};


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
//========================|

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

struct SCC_vertex
{
	int node;
	int sccID;
	void clear()
	{}
};

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
	//Constructor
	RT_Vertex()
	{
		Root = -1;
		Parent = -1;
		EDGwt = 0.0;
		Level = -1;
		marked = -1;

		/*PossN.clear();*/
		Dist = 0.0;
		Update = false;
	}

	//Destructor
	void clear()
	{}



};
//The Rooted tree is a vector of structure RT_Vertex;