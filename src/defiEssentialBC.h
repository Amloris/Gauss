/* defiEssentialBC
This class stores the essential (prescribed displacement) boundary conditions.
*/

#ifndef defiEssentialBC_h
#define defiEssentialBC_h

#include "defiDof.h"
#include "defiNode.h"

#include "globAccessItems.h"

using namespace std;

class defiEssentialBC
{
private:
	defiNode* d_node;	//	ptr to node (class defiNode) associated with an essential BC 
	DOFType d_dof;		//	enum DOFType {UX=0, UY}; d_dof is either UX or UY
	double d_value;		//	d_value = prescribed displ. value

	defiEssentialBC(); 	//never used constructors
	defiEssentialBC(defiEssentialBC &ebc);

public:
	// Constructors
	defiEssentialBC(defiNode* node, DOFType dof, double value);
	//Destructor
	~defiEssentialBC();

	// Functions
	defiNode* getNode() const;	//	return d_node;
	DOFType getDof() const;		//	return d_dof;
	double getValue() const;	//	return d_value;
	void printData() const;	//print useful info, e.g UY and UY
};


//Functions
//-----------------------------------------------------------------------------
defiEssentialBC::defiEssentialBC()  { std::cout << "Creating defiEssentialBC Class Object" << endl; }
defiEssentialBC::~defiEssentialBC() { std::cout << "Deleting defiEssentialBC Class Object" << endl; }
defiEssentialBC::defiEssentialBC(defiNode* node, DOFType dof, double value)
{
	d_node = node;
	d_dof = dof;
	d_value = value;
}

defiNode* defiEssentialBC::getNode() const
{
	return d_node;
}

DOFType defiEssentialBC::getDof() const
{
	return d_dof;
}

double defiEssentialBC::getValue() const
{
	return d_value;
}

void defiEssentialBC::printData() const
{
	//Get Data
	int id;                      //Node ID
	id = (*getNode()).getID();   //Get Node ID

	//Print to File
	fout << "    " << id << "    " << d_dof << "    " << d_value << endl;
}

#endif