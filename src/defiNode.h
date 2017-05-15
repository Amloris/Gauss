/* defiNode
The defiNode class stores a node id, the coordinates, equation numbers for
each degree of freedom (dof), the solution, and boundary condition data.
In this implementation, we allow two dof's at the node, UX and UY.
The data for these dofs is stored in arrays.  The appropriate value is
accessed using an enumerated type as follows:
enum DOF {UX=0, UY};

For example, we can then return the UX dof value by:
return d_dof[UX];
*/

#ifndef defiNode_h
#define defiNode_h

#include "defiDof.h"
#include "calcStress2D.h"

#include "globAccessItems.h"

using namespace std;

class defiNode
{
private:
	int d_id;	//node id
	double d_x, d_y;	//coordinates x and y of this node
	defiDof* d_dof[2];	//each node contains two dofs. Therefore two pointers are used
	calcStress2D* d_nodeStress;	//store nodal stress (averaged over all neighboring element) [new addition] 

	defiNode();	// Never used constructors
	defiNode(const defiNode& node);

public:
	// Constructors
	defiNode(int id);
	// Destructor
	~defiNode();

	// Functions
	int getID() const;					//return d_id;
	void setCoords(double x, double y);	//set d_x = x, d_y = y;
	double getX() const;				//return d_x;
	double getY() const;				//return d_y;
	defiDof* getDof(DOFType dof) const;	//return d_dof[dof];
	calcStress2D* getStress() const;    //return d_nodeStress;
	void printData() const;	            //print useful info
	void printResults() const;	        //print useful results
};


//Functions
//-----------------------------------------------------------------------------
defiNode::~defiNode() { std::cout << "Deleting defiNode Class Object" << endl; }
defiNode::defiNode(int id)
{
	d_id = id;
}

int defiNode::getID() const
{
	return d_id;
}

void defiNode::setCoords(double x, double y)
{
	d_x = x;
	d_y = y;
}

double defiNode::getX() const
{
	return d_x;
}

double defiNode::getY() const
{
	return d_y;
}

void defiNode::printData() const
{
	fout << "    " << d_id << "    " << d_x << "    " << d_y << endl;
}

#endif