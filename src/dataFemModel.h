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

using namespace std;

class dataFemModel
{
private:
	defiNode** d_nodes;		           //Array of ptrs to node class (defiNode). The array size is equal to the total number of nodes in the FEM model.
	defiMaterial** d_matprops;         //Array of ptrs to material class (defiMaterial). The array size is equal to the total number of materials in the FEM model.
	defiElementQ8** d_elems;		   //Array of ptrs to element class (defiElementQ8). The array size is equal to the total number of elements in the FEM model.
	defiEssentialBC** d_essentialBCs;  //Array of ptrs to essential B.C. class (defiEssentialBC). The array size is equal to the total number of essential BCs in the FEM model.
	defiPointBC** d_pointBCs;	       //Array of prts to point force B.C. class (defiPointBC). The array size is equal to the total number of point force BCs in the FEM model.
	defiNaturalBC** d_naturalBCs;	   //Array of prts to natural B.C. class (defiNaturalBC). The array size is equal to the total number of natural BCs in the FEM model.
	string d_title;		               //Title of the fem analysis (from the input data file)
	int d_numNodes = NULL;	           //Total number of nodes of the fem model
	int d_numElems = NULL;	           //Total number of elements of the fem model
	int d_numMats  = NULL;		       //Total number of materials of the fem model
	int d_probType = NULL;	           //Problem type: =0 for axisymmetric (irrelevant in this course); =1 for plane stress; =2 for plane strain.
	int d_numEssentialBCs = NULL;	   //Number of essential B.C. of the fem model
	int d_numPointBCs = NULL;		   //Number of point force B.C.
	int d_numNaturalBCs = NULL;		   //Number of natural B.C
	int d_neq;					       //Total number of system equations (which is equal to number of nodes * 2 - number of essential BCs)

public:
	// Constructors
	dataFemModel();
	//Destructor
	~dataFemModel();

	// Functions
	void readData();	               //Read all input data from input file
	void writeData();	               //Write the output data, see "fem.out" of benchmark problems for format
	int getNumNodes();				   //Return d_numNodes
	int getNumElems();				   //Return d_numElems
	int getNumMats();                  //Return d_numMats
	int getProbType();				   //Return d_probType
	int getNumEssentialBCs();		   //Return d_numEssentialBCs
	int getNumPointBCs();			   //Return d_numPointBCs
	int getNumNaturalBCs();			   //Return d_numNaturalBCs
	int getNumEq();                    //Return d_neq
	defiNode *getNode(int node);			     //Return d_nodes[node] (node is node id)
	defiElementQ8 *getElem(int elem);			 //Return d_elems[elem] (elem is element id)
	defiMaterial *getMat(int mat);               //Return d_matprops[mat] (mat is material id)  Passed by pointer
	defiEssentialBC *getEssentialBC(int ebc);  	 //Return d_essentialBCs[ebc] (ebc is the id for the essential BC)
	defiPointBC *getPointBC(int pbc);			 //Return d_pointBCs[pbc] (pbc is the id for the point force boundary condition)
	defiNaturalBC *getNaturalBC(int nbc);		 //Return d_naturalBCs[nbc] (nbc is the id for the natural boundary condition)
	void setNumEquation(int value);	   //Set d_neq = value
	void writeResults();	           //Output results data, see "fem.out" of benchmark problems for format
	void writePlotFile();	           //Write the output data for contour plot program "ConPlot". see "fem.out" of benchmark problems for format.
};


//Functions
//-----------------------------------------------------------------------------
dataFemModel::dataFemModel()  { cout << "Creating dataFemModel Class Object" << endl; }
dataFemModel::~dataFemModel() { cout << "Deleting dataFemModel Class Object" << endl; }

int dataFemModel::getNumNodes()		   { return d_numNodes; }	      //Return d_numNodes
int dataFemModel::getNumElems()		   { return d_numElems; }	      //Return d_numElems
int dataFemModel::getNumMats()		   { return d_numMats;  }         //Return d_numMats     
int dataFemModel::getProbType()		   { return d_probType; }	      //Return d_probType
int dataFemModel::getNumEssentialBCs() { return d_numEssentialBCs; }  //Return d_numEssentialBCs
int dataFemModel::getNumPointBCs()     { return d_numPointBCs; }	  //Return d_numPointBCs
int dataFemModel::getNumNaturalBCs()   { return d_numNaturalBCs; }	  //Return d_numNaturalBCs
int dataFemModel::getNumEq()		   { return d_neq; }              //Return d_neq

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

	if (d_matprops == NULL) {
		ferr<< "ERROR::MATERIAL_CLASS::MEMORY_ALLOCATION" << endl;
		exit(0);
	}
	if (d_elems == NULL) {
		ferr << "ERROR::ELEMENT_CLASS::MEMORY_ALLOCATION" << endl;
		exit(0);
	}
	if (d_nodes == NULL) {
		ferr << "ERROR::NODE_CLASS::MEMORY_ALLOCATION" << endl;
		exit(0);
	}
	if (d_probType == NULL) {
		ferr << "ERROR::PROBLEM_TYPE::DEFINITION" << endl;
		exit(0);
	}

	//Fill Material Arrays
	for (int i = 0; i < d_numMats; i++)
	{
		int id;		    //Material ID
		double e;		//Young's modulus
		double nu;	    //Poisson's ratio
		double k;		//Conductivity, irrelevant for us
		double thick;   //Thickness
		int probType;   //ProblemType

		//Load Data
		fin >> id >> e >> nu >> k >> thick >> probType;                    //Load info
		d_matprops[i] = new defiMaterial(id, e, nu, k, thick, d_probType);   //Push info
		fin.ignore(numeric_limits<streamsize>::max(), '\n');               //Skip to end of line
	}

	//Fill Element Arrays
	for (int i = 0; i < d_numNodes; i++)         //Assign node id's
	{
		d_nodes[i] = new defiNode(i + 1);         
	}

	for (int i = 0; i < d_numElems; i++)         //Assign element data
	{
		int eID;						         //Element ID
		int mat;                                 //Material ID
		int const numNodes = 8;                  //8 nodes for a q8 element

		defiMaterial* eMat;                      //Link to the global materials
		defiNode* eNodes[numNodes];              //Link to the global nodes


		//Load Data
		fin >> eID >> mat;

		if (i != (eID - 1))
		{
			ferr << "ERROR::ELEMENT_CLASS::ELEMENT_NUMBERING" << endl;
			exit(0);
		}
		if (mat > d_numMats)
		{
			ferr << "ERROR::MATERIAL_CLASS::MATERIAL_NUMBERING" << endl;
			exit(0);
		}

		eMat = d_matprops[mat - 1];              //Link to the material

		for (int j = 0; j < numNodes; j++)
		{
			int nodeID;
			fin >> nodeID;                       //Obtain node id
			eNodes[j] = d_nodes[nodeID - 1];     //Link to Nodes
		}

		d_elems[i] = new defiElementQ8(eID, eMat, numNodes, eNodes);      //Build q8 element from loaded info
	}

	//Fill Node Arrays
	for (int i = 0; i < d_numNodes; i++)
	{
		int id;
		double x, y;

		//Load Data
		fin >> id >> x >> y;                     //Load info
		d_nodes[i]->setCoords(x, y);             //Push coordinates
	}

	//Read Boundary Condition Info
	fin >> d_numEssentialBCs;
	fin >> d_numPointBCs;
	fin >> d_numNaturalBCs;

	//Dynamically Allocate Arrays   (Essential, Point, Natural)
	d_essentialBCs = new defiEssentialBC*[d_numEssentialBCs];
	d_pointBCs = new defiPointBC*[d_numPointBCs];
	d_naturalBCs = new defiNaturalBC*[d_numNaturalBCs];

	if ((d_essentialBCs == NULL) && (d_numEssentialBCs > 0)) {
		ferr << "ERROR::ESSENTIALBC_CLASS::MEMORY_ALLOCATION" << endl;
		exit(0);
	}
	if ((d_pointBCs == NULL) && (d_numPointBCs > 0)) {
		ferr << "ERROR::POINTBC_CLASS::MEMORY_ALLOCATION" << endl;
		exit(0);
	}
	if ((d_naturalBCs == NULL) && (d_numNaturalBCs > 0)) {
		ferr << "ERROR::NATURALBC_CLASS::MEMORY_ALLOCATION" << endl;
		exit(0);
	}

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

		d_essentialBCs[id_bc - 1] = new defiEssentialBC(d_nodes[id-1], dof, value);
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

		d_pointBCs[id_bc - 1] = new defiPointBC(d_nodes[id-1], dof, value);
	}

	//Fill NaturalBC Arrays
	for (int i = 0; i < d_numNaturalBCs; i++)
	{
		int id_bc;               //BC ID
		int id1, id2, id3;       //Node ID's
		double value;            //Normal Traction

		//Load Data
		fin >> id_bc >> id1 >> id2 >> id3 >> value;

		d_naturalBCs[id_bc - 1] = new defiNaturalBC(d_nodes[id1-1], d_nodes[id2 - 1], d_nodes[id3 - 1], value);
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
	fout << "|                Gauss : A 2-D Finite Elements Program                    |" << endl;
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
		d_matprops[i]->printData();
	}

	//Print Element Data
	fout << " ***Element Data***" << endl;
	fout << "ElemID  " << "Material      " << "Node Connectivity" << endl;
	for (int i = 0; i < d_numElems; i++)
	{
		d_elems[i]->printData();
	}

	//Print Node Data
	fout << endl;
	fout << "***Nodal Data***" << endl;
	for (int i = 0; i < d_numNodes; i++)
	{
		d_nodes[i]->printData();
	}

	//Print EssentialBC Data
	fout << endl;
	fout << "***Essential BC Data***" << endl;
	for (int i = 0; i < d_numEssentialBCs; i++)
	{
		d_essentialBCs[i]->printData();
	}

	//Print PointBC Data
	fout << endl;
	fout << "***Point BC Data***" << endl;
	for (int i = 0; i < d_numPointBCs; i++)
	{
		d_pointBCs[i]->printData();
	}

	//Print NaturalBC Data
	fout << endl;
	fout << "***Natural BC Data***" << endl;
	for (int i = 0; i < d_numNaturalBCs; i++)
	{
		d_naturalBCs[i]->printData();
	}
}

void dataFemModel::writeResults()
{	//Records the processed results into the output file

	//Print Number of Equations
	fout << endl << endl << endl;
	fout << " ******************** Results *************************************" << endl << endl;
	fout << "     Total number of equations = " << d_neq << endl << endl;

	//Print Element Data
	fout << " ********************* Element Data *******************************" << endl << endl;
	fout << "                  (Stresses at element center)" << endl;
	fout << " Elem ID  Material   SigXX         SigYY         SigXY         SigZZ" << endl;

	for (int i = 0; i < d_numElems; i++)
	{
		calcStress2D stressGPs[9];                //Holder for stresses at nodal positions
		getElem(i)->getNodalData(stressGPs);      //Retrieve stress data for element

		//Get Element Data
		int id, mat;
		id = i + 1;                               //Element number  
		mat = getElem(i)->getMaterial()->getID(); //Material ID


		//Get Stress Data
		double sigxx, sigyy, sigzz, sigxy;        //Stress values at center
		sigxx = stressGPs[8].getSigXX();
		sigyy = stressGPs[8].getSigYY();
		sigzz = stressGPs[8].getSigZZ();
		sigxy = stressGPs[8].getSigXY();

		//Export Data
		fout << "    " << id << "       " << mat;
		fout.setf(ios::scientific);
		fout.precision(4);
		fout.width(17);
		fout << sigxx;
		fout << "    " << sigyy;
		fout << "    " << sigxy;
		fout << "    " << sigzz;
		fout << endl;
	}
	fout << endl << endl;

	//Print Displacement Data
	fout << " ********************** Nodal Data ********************************" << endl << endl;
	fout << " Node ID        X              Y             UX            UY " << endl;
	for (int i = 0; i < d_numNodes; i++)
	{
		int id;                        //Node id
		double x, y, u, v;             //Nodal positions (x,y), Nodal displacements (u,v)

		//Get Data
		id = i+1;
		x = d_nodes[i]->getX();
		y = d_nodes[i]->getY();
		u = d_nodes[i]->getDof(UX)->getValue();
		v = d_nodes[i]->getDof(UY)->getValue();

		//Export Data
		fout << "    " << id;
		fout.setf(ios::scientific);
		fout.precision(4);
		fout.width(17);
		fout << x;
		fout << "    " << y;
		fout << "    " << u;
		fout << "    " << v;
		fout << endl;
	}

	//Print Nodal Stresses
	fout << endl << " Node ID      SigXX           SigYY         SigXY         SigZZ" << endl;
	for (int i = 0; i < d_numNodes; i++)
	{
		int id;                                    //Node id 
		double sigxx, sigyy, sigzz, sigxy;         //Nodal Stresses

		//Get Data
		id = i + 1;
		sigxx = d_nodes[i]->getStress()->getSigXX();
		sigyy = d_nodes[i]->getStress()->getSigYY();
		sigzz = d_nodes[i]->getStress()->getSigZZ();
		sigxy = d_nodes[i]->getStress()->getSigXY();

		//Export Data
		fout << "    " << id;
		fout.setf(ios::scientific);
		fout.precision(4);
		fout.width(17);
		fout << sigxx;
		fout << "    " << sigyy;
		fout << "    " << sigxy;
		fout << "    " << sigzz;
		fout << endl;
	}

	//Termination Message
	fout << endl << " ***** End of results ***** ";

}

void dataFemModel::writePlotFile()
{	//Records the processed results into the conplot file

	//Header
	fplot << "%1%" << endl;
	fplot << "5" << " " << d_numNodes << " " << d_numElems << endl;
	fplot << "SIGXX" << endl;
	fplot << "SIGYY" << endl;
	fplot << "SIGZZ" << endl;
	fplot << "SIGXY" << endl;
	
	//Print Node Postions and Displacments
	fplot << "SIG1" << endl;
	for (int i = 0; i < d_numNodes; i++)
	{
		//Get Data
		double x, y, u, v;             //Nodal positions (x,y), Nodal displacements (u,v)
		x = d_nodes[i]->getX();
		y = d_nodes[i]->getY();
		u = d_nodes[i]->getDof(UX)->getValue();
		v = d_nodes[i]->getDof(UY)->getValue();

		//Export Data
		fplot.setf(ios::scientific);
		fplot.precision(4);
		fplot.width(12);
		fplot << x;
		fplot << "    " << y;
		fplot << "    " << u;
		fplot << "    " << v;
		fplot << endl;
	}

	//Print Stresses for Nodes in Each Element
	for (int i = 0; i < d_numElems; i++)
	{
		//Data Holders
		int const eNumNodes = 8;                 //Number of nodes in q8 element
		defiNode* eNodes[eNumNodes];             //Holder for nodes
		getElem(i)->getNodes(eNodes);            //Get list of nodes

		//Get Nodes in Element and Rearange
		defiVector Nodes(8);
		for (int j = 0; j < getElem(i)->getNumNodes()/2; j++)
		{
			//Record first 4
			int id = eNodes[2*j]->getID();
			Nodes.setCoeff(j, id);

			//Record last 4
			id = eNodes[2*j+1]->getID();
			Nodes.setCoeff(4 + j, id);
		}

		//Print Nodes in Element
		int num = getElem(i)->getNumNodes();
		fplot << num << " ";
		for (int k = 0; k < getElem(i)->getNumNodes(); k++)
		{
			int temp;
			temp = Nodes.getCoeff(k);
			fplot << temp << " ";
		}
		fplot << endl;

		//Print Nodal Stresses
		for (int q = 0; q < getElem(i)->getNumNodes(); q++)
		{
			int id;                                    //Node id 
			double sigxx, sigyy, sigzz, sigxy;         //Nodal Stresses

			id = Nodes.getCoeff(q) - 1;

			//Get Data
			sigxx = d_nodes[id]->getStress()->getSigXX();
			sigyy = d_nodes[id]->getStress()->getSigYY();
			sigzz = d_nodes[id]->getStress()->getSigZZ();
			sigxy = d_nodes[id]->getStress()->getSigXY();

			//Export Data
			fplot.setf(ios::scientific);
			fplot.precision(4);
			fplot.width(11);
			fplot << sigxx;
			fplot << "    " << sigyy;
			fplot << "    " << sigzz;
			fplot << "    " << sigxy;
			fplot << "    " << sigxx;
			fplot << endl;
		}

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