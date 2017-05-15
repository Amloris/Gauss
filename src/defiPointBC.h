/* defiPointBC
The defiPointBC class stores point force boundary condition data. The
point force boundary condition is concentrated force at a node, a special
case of natural boundary condition data.
*/

#ifndef defiPointBC_h
#define defiPointBC_h

#include "defiDof.h"
#include "defiNode.h"

#include "globAccessItems.h"

using namespace std;

enum PointBCType { FX = 0, FY };

class defiPointBC	//this is the point force B.C. class
{
private:
	defiNode* d_node;		//ptrs to a node which has point force b.c.
	PointBCType d_pbctype;	//type of point force, i.e. FX or FY
	double d_value;			//value of the point force

	defiPointBC();	//never used constructors
	defiPointBC(defiPointBC &ebc);

public:
	// Constructors
	defiPointBC(defiNode* node, PointBCType dof, double value);
	//Destructor
	~defiPointBC();

	// Functions
	defiNode* getNode() const;			// return d_node;
	PointBCType getPbctype() const;		// return d_pbctype;
	double getValue() const;			// return d_value;
	void printData() const;	// print useful data
};


//Functions
//-----------------------------------------------------------------------------
defiPointBC::defiPointBC()  { std::cout << "Creating defiPointBC Class Object" << endl; }
defiPointBC::~defiPointBC() { std::cout << "Deleting defiPointBC Class Object" << endl; }
defiPointBC::defiPointBC(defiNode* node, PointBCType dof, double value)
{
	d_node = node;
	d_pbctype = dof;
	d_value = value;
}

defiNode* defiPointBC::getNode() const
{
	return d_node;
}

PointBCType defiPointBC::getPbctype() const
{
	return d_pbctype;
}

double defiPointBC::getValue() const
{
	return d_value;
}

void defiPointBC::printData() const
{
	//Get Data
	int id;                      //Node ID
	id = (*getNode()).getID();   //Get Node ID

	//Print to File
	fout << "    " << id << "    " << d_pbctype << "    " << d_value << endl;
}


#endif