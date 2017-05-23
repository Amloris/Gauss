/* defiPointBC
The defiPointBC class stores point force boundary condition data. The
point force boundary condition is concentrated force at a node, a special
case of natural boundary condition data.
*/

#ifndef defiPointBC_h
#define defiPointBC_h

#include "defiDof.h"
#include "defiNode.h"

using namespace std;

enum PointBCType { FX = 0, FY };

class defiPointBC	                   //This is the point force B.C. class
{
private:
	defiNode* d_node;		           //Ptrs to a node which has point force b.c.
	PointBCType d_pbctype;	           //Type of point force, i.e. FX or FY
	double d_value;			           //Value of the point force

	defiPointBC();	                   //Never used constructors
	defiPointBC(defiPointBC &ebc);

public:
	// Constructors
	defiPointBC(defiNode* node, PointBCType dof, double value);
	//Destructor
	~defiPointBC();

	// Functions
	defiNode* getNode() const;		   //Return d_node;
	PointBCType getPbctype() const;	   //Return d_pbctype;
	double getValue() const;		   //Return d_value;
	void printData() const;	           //Print useful data
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
	id = getNode()->getID();     //Get Node ID

	//Print to File
	fout << "    " << id << "    " << d_pbctype << "    " << d_value << endl;
}


#endif