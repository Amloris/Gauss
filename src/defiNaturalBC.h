/* defiNaturalBC
The defiNaturalBC class stores natural (prescribed traction) boundary condition data.
Data is given for the three nodes on an element edge where traction is
applied. The positive tangential direction is to traverse the boundary
counter-clockwisely (the material is on the left side of the travelling
direction). The positive tangential traction direction is along the
positive tangential direction, while the positive normal traction is
pointing away from the element.
*/

#ifndef defiNaturalBC_h
#define defiNaturalBC_h

#include "defiDof.h"
#include "defiNode.h"
#include "defiVector.h"

#include "globAccessItems.h"

using namespace std;

class defiNaturalBC	              //This is for natural (traction) boundary condition
{
private:
	defiNode* d_node1;	          //For Q8 element, there are three nodes for each side. this is node 1
	defiNode* d_node2;	          //Node2
	defiNode* d_node3;	          //Node3
	double d_value;		          //The constant normal traction value

	defiNaturalBC();	          //Never used constructors
	defiNaturalBC(defiNaturalBC &ebc);

public:
	// Constructors
	defiNaturalBC(defiNode* node1, defiNode* node2, defiNode* node3, double value);
	//Destructor
	~defiNaturalBC();
	// Functions
	void printData() const;		  //Print useful info
	void getNaturalBCF(defiVector &f, defiDof* dofs[]);
	     //Calculate the equivalent nodal forces from given traction. 
	     //Used to form the right hand side load vector
};


//Functions
//-----------------------------------------------------------------------------
defiNaturalBC::defiNaturalBC()  { std::cout << "Creating defiNaturalBC Class Object" << endl; }
defiNaturalBC::~defiNaturalBC() { std::cout << "Deleting defiNaturalBC Class Object" << endl; }
defiNaturalBC::defiNaturalBC(defiNode* node1, defiNode* node2, defiNode* node3, double value) 
{
	d_node1 = node1;
	d_node2 = node2;
	d_node3 = node3;
	d_value = value;
}

void defiNaturalBC::printData() const
{
	//Get Data
	int id1, id2, id3;

	id1 = d_node1->getID();
	id2 = d_node2->getID();
	id3 = d_node3->getID();

	//Print to File
	fout << "    " << id1 << "    " << id2 << "    " << id3 << "    "  << d_value << endl;
}



#endif