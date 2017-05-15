/* dataFemModel
The dataFemModel class stores all problem data.  It provides a means to store and access
all elements, nodes, material, and boundary condition data.  Memory for
all data is allocated dynamically after the size of the required storage
is known.  Simple arrays of pointers are used for the nodes, elements and
materials.  For example, nodes are dynamically allocated by:
d_nodes = new defiNode*[d_numNodes];
...
for(...)
{
d_nodes[id - 1] = new defiNode(id);
}

Elements are dynamically allocated by:
d_elems = new defiElementQ8*[d_numElems];
...
for(...)
{
d_elems[i] = new defiElementQ8(eID, mat, numEnodes, enodes);
}
*/

#ifndef dataFemModel_h
#define dataFemModel_h

#include "defiNode.h"
#include "defiMaterial.h"
#include "defiElementQ8.h"
#include "defiEssentialBC.h"
#include "defiPointBC.h"
#include "defiNaturalBC.h"

//#include "globAccessItems.h"

using namespace std;

class dataFemModel
{
private:
	defiNode** d_nodes;		           //array of ptrs to node class (defiNode). The array size is equal to the total number of nodes in the FEM model.
	defiMaterial** d_matprops;         //array of ptrs to material class (defiMaterial). The array size is equal to the total number of materials in the FEM model.
	defiElementQ8** d_elems;		   //array of ptrs to element class (defiElementQ8). The array size is equal to the total number of elements in the FEM model.
	defiEssentialBC** d_essentialBCs;  //array of ptrs to essential B.C. class (defiEssentialBC). The array size is equal to the total number of essential BCs in the FEM model.
	defiPointBC** d_pointBCs;	       //array of prts to point force B.C. class (defiPointBC). The array size is equal to the total number of point force BCs in the FEM model.
	defiNaturalBC** d_naturalBCs;	   //array of prts to natural B.C. class (defiNaturalBC). The array size is equal to the total number of natural BCs in the FEM model.
	string d_title;		               //title of the fem analysis (from the input data file)
	int d_numNodes;	                   //total number of nodes of the fem model
	int d_numElems;	                   //total number of elements of the fem model
	int d_numMats;		               //total number of materials of the fem model
	int d_probType;	                   //problem type: =0 for axisymmetric (irrelevant in this course); =1 for plane stress; =2 for plane strain.
	int d_numEssentialBCs;		       //number of essential B.C. of the fem model
	int d_numPointBCs;			       //number of point force B.C.
	int d_numNaturalBCs;		       //number of natural B.C
	int d_neq;					       //total number of system equations (which is equal to number of nodes * 2 - number of essential BCs)

public:
	// Constructors
	dataFemModel();
	//Destructor
	~dataFemModel();

	// Functions
	void readData();	               //read all input data from input file
	void writeData();	               //write the output data, see "fem.out" of benchmark problems for format
	int getNumNodes();				   //return d_numNodes
	int getNumElems();				   //return d_numElems
	int getNumMats();                  //return d_numMats
	int getProbType();				   //return d_probType
	int getNumEssentialBCs();		   //return d_numEssentialBCs
	int getNumPointBCs();			   //return d_numPointBCs
	int getNumNaturalBCs();			   //return d_numNaturalBCs
	int getNumEq();                    //return d_neq
	defiNode *getNode(int node);			     //return d_nodes[node] (node is node id)
	defiElementQ8 *getElem(int elem);			 //return d_elems[elem] (elem is element id)
	defiMaterial *getMat(int mat);               //return d_matprops[mat] (mat is material id)  Passed by pointer
	defiEssentialBC *getEssentialBC(int ebc);  	 //return d_essentialBCs[ebc] (ebc is the id for the essential BC)
	defiPointBC *getPointBC(int pbc);			 //return d_pointBCs[pbc] (pbc is the id for the point force boundary condition)
	defiNaturalBC *getNaturalBC(int nbc);		 //return d_naturalBCs[nbc] (nbc is the id for the natural boundary condition)
	void setNumEquation(int value);	   //set d_neq = value
	void writeResults();	           //output results data, see "fem.out" of benchmark problems for format
	void writePlotFile();	           //write the output data for contour plot program "ConPlot". see "fem.out" of benchmark problems for format.
};


//Functions
//-----------------------------------------------------------------------------
dataFemModel::dataFemModel()  { std::cout << "Creating dataFemModel Class Object" << endl; }
dataFemModel::~dataFemModel() { std::cout << "Deleting dataFemModel Class Object" << endl; }

int dataFemModel::getNumNodes() { return d_numNodes; }	   //return d_numNodes
int dataFemModel::getNumElems() { return d_numElems; }	   //return d_numElems
int dataFemModel::getNumMats()  { return d_numMats;  }     //return d_numMats     
int dataFemModel::getProbType() { return d_probType; }	   //return d_probType
int dataFemModel::getNumEssentialBCs() { return d_numEssentialBCs; }      //return d_numEssentialBCs
int dataFemModel::getNumPointBCs()     { return d_numPointBCs; }		  //return d_numPointBCs
int dataFemModel::getNumNaturalBCs()   { return d_numNaturalBCs; }	      //return d_numNaturalBCs
int dataFemModel::getNumEq()    { return d_neq; }          //return d_neq

void dataFemModel::readData()
{   //Reads the model data and dynamically allocates memory for data storage

	//Read Title
	string temp;
	getline(fin,temp);                            
	if (temp[0] == '/') { d_title = temp; }      //Read the first line of the input file and set it as the title
	else { temp = "No Title"; d_title = temp; }  //If no title exists, record 'No Title'

	//Read Problem Size Info
	fin >> d_numNodes;
	fin >> d_numElems;
	fin >> d_numMats;
	fin >> d_probType;

	//Dynamically Allocate Arrays   (Materials, Elements, Nodes)
	d_matprops = new defiMaterial*[d_numMats];   //Material Storage
	d_elems = new defiElementQ8*[d_numElems];    //Element Storage
	d_nodes = new defiNode*[d_numNodes];         //Node Storage

	//Fill Material Arrays
	for (int i = 0; i < d_numMats; i++)
	{
		int id;		    //material ID
		double e;		//Young's modulus
		double nu;	    //Poisson's ratio
		double k;		//conductivity, irrelevant for us
		double thick;   //thickness
		int probType;   //problemType

		//Load Data
		fin >> id >> e >> nu >> k >> thick >> probType;                    //Load info
		d_matprops[i] = new defiMaterial(id, e, nu, k, thick, probType);   //Push info
		fin.ignore(numeric_limits<streamsize>::max(), '\n');               //Skip to end of line
	}

	//Fill Element Arrays
	for (int i = 0; i < d_numElems; i++)
	{
		int eID;                       //Element ID
		int mat;                       //Material ID
		int numNodes = 8;              //8 nodes for a q8 element

		defiNode** enodes;
		enodes = new defiNode*[numNodes];

		//Load Data
		fin >> eID >> mat;
		for (int i = 0; i < numNodes; i++)
		{
			int temp;
			fin >> temp;

			d_nodes[temp] = new defiNode(temp);

			enodes[i] = d_nodes[temp];
		}

		d_elems[i] = new defiElementQ8(eID, getMat(mat-1), numNodes, enodes);
		//defiElementQ8(int id, defiMaterial* mat, int numNodes, defiNode* enodes[]) 
	}

	//Fill Node Arrays
	for (int i = 0; i < d_numNodes; i++)
	{
		int id;
		double x, y;

		//Load Data
		fin >> id >> x >> y;                     //Load info
		d_nodes[id - 1] = new defiNode(id);      //Push ID
		(*getNode(i)).setCoords(x,y);            //Push coordinates
	}

	//Read Boundary Condition Info
	fin >> d_numEssentialBCs;
	fin >> d_numPointBCs;
	fin >> d_numNaturalBCs;

	//Dynamically Allocate Arrays   (Essential, Point, Natural)
	d_essentialBCs = new defiEssentialBC*[d_numEssentialBCs];
	d_pointBCs = new defiPointBC*[d_numPointBCs];
	d_naturalBCs = new defiNaturalBC*[d_numNaturalBCs];

	//Fill EssentialBC Arrays
	for (int i = 0; i < d_numEssentialBCs; i++)
	{
		int id_bc;       //BC ID  
		int id;          //Node ID
		double value;    //Displacement Value
		string temp;     //Hold Enum

		//Load Data
		fin >> id_bc >> id >> temp >> value;

		//Convert from enum to dof
		DOFType dof;
		if (temp[0] == 'U' && temp[1] == 'X') { dof = UX; }
		else { dof = UY; }

		d_essentialBCs[id_bc - 1] = new defiEssentialBC((getNode(id - 1)), dof, value);	
	}

	//Fill PointBC Arrays
	for (int i = 0; i < d_numPointBCs; i++)
	{
		int id_bc;
		int id;
		double value;
		string temp;

		//Load Data
		fin >> id_bc >> id >> temp >> value;

		//Convert from enum to dof
		PointBCType dof;
		if (temp[0] == 'F' && temp[1] == 'X') { dof = FX; }
		else { dof = FY; }

		d_pointBCs[id_bc - 1] =	new defiPointBC((getNode(id - 1)), dof, value);
	}

	//Fill NaturalBC Arrays
	for (int i = 0; i < d_numNaturalBCs; i++)
	{
		int id_bc;               //BC ID
		int id1, id2, id3;       //Node ID's
		double value;            //Normal Traction

		//Load Data
		fin >> id_bc >> id1 >> id2 >> id3 >> value;

		d_naturalBCs[id_bc - 1] = new defiNaturalBC(getNode(id1 - 1), getNode(id2 - 1), getNode(id3 - 1), value);
	}

	//Set Number of System Equations
	int val;
	val = d_numNodes * 2 - d_numEssentialBCs;
	setNumEquation(val);
};

void dataFemModel::writeData()
{	//Prints the loaded input data to fout for verification

	//Print Program
	fout << "+-------------------------------------------------------------------------+" << endl;
	fout << "|              Regolith: A 2-D Finite Elements Program                    |" << endl;
	fout << "+-------------------------------------------------------------------------+" << endl << endl;

	//Print Title
	fout << "***Problem Title***" << endl;
	fout << d_title << endl << endl;

	//Print Basic Problem Info
	fout << "***Problem Size***" << endl;
	fout << "Nodes = " << d_numNodes;
	fout << "  Elements = " << d_numElems;
	fout << "  Materials = " << d_numMats;
	fout << "  Problem Type = " << d_probType;
	if (d_probType == 1) { fout << "(Plane Stress)"; }
	if (d_probType == 2) { fout << "(Plane Strain)"; }
	fout << endl;

	//Print Material Data
	fout << endl;
	fout << "***Material Data***" << endl;
	for (int i = 0; i < d_numMats; i++)
	{
		(*getMat(i)).printData();                //Obtain Material Object
	}

	//Print Element Data
	fout << " ***Element Data***" << endl;
	fout << "ElemID  " << "Material      " << "Node Connectivity" << endl;
	for (int i = 0; i < d_numElems; i++)
	{
		(*getElem(i)).printData();               //Obtain Element Object
	}

	//Print Node Data
	fout << endl;
	fout << "***Nodal Data***" << endl;
	for (int i = 0; i < d_numNodes; i++)
	{
		(*getNode(i)).printData();               //Obtain Node Object
	}

	//Print EssentialBC Data
	fout << endl;
	fout << "***Essential BC Data***" << endl;
	for (int i = 0; i < d_numEssentialBCs; i++)
	{
		(*getEssentialBC(i)).printData();
	}

	//Print PointBC Data
	fout << endl;
	fout << "***Point BC Data***" << endl;
	for (int i = 0; i < d_numPointBCs; i++)
	{
		(*getPointBC(i)).printData();
	}

	//Print NaturalBC Data
	fout << endl;
	fout << "***Natural BC Data***" << endl;
	for (int i = 0; i < d_numNaturalBCs; i++)
	{
		(*getNaturalBC(i)).printData();
	}
}

void dataFemModel::setNumEquation(int value)	   //set d_neq = value
{
	d_neq = value;
}

defiMaterial* dataFemModel::getMat(int mat)                            
{													       //Passes the pointer to the material object
	return d_matprops[mat];                                //Return d_matprops[mat] (mat is material id)
}

defiNode* dataFemModel::getNode(int node)
{                                                          //Passes the pointer to the node object
	return d_nodes[node];                                  //Return d_nodes[node] (node is the node id)
}

defiElementQ8* dataFemModel::getElem(int elem)
{                                                          //Passes the pointer to the element object
	return d_elems[elem];                                  //Return d_elems[elems] (elems is element id)
}

defiEssentialBC* dataFemModel::getEssentialBC(int ebc)  	
{														   //Passes the pointer to the essentialBC object
	return d_essentialBCs[ebc];                            //Return d_essentialBCs[ebc] (ebc is the id for the essential BC)
}

defiPointBC* dataFemModel::getPointBC(int pbc)			 
{														   //Passes the pointer to the pointBC object
	return d_pointBCs[pbc];                                //Return d_pointBCs[pbc] (pbc is the id for the point force boundary condition)
}

defiNaturalBC* dataFemModel::getNaturalBC(int nbc)		 
{														   //Passes the pointer to the naturalBC object	
	return d_naturalBCs[nbc];                              //Return d_naturalBCs[nbc] (nbc is the id for the natural boundary condition)
}

#endif