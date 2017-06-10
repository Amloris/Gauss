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

void defiNaturalBC::getNaturalBCF(defiVector &f, defiDof* dofs[])
{	//Calculate the equivalent nodal forces from given traction. 
	//Used to form the right hand side load vector

	//Get Tractions
	double tn, ts;                //Surface tractions (Normal, Shear)
	tn = d_value;                 //Normal surface force
	ts = 0.0;                     //Shear surface force (Zero for our case, but it can be non-zero if we allow it)

	//Initial Values
	f.zero();                     //Zero the force vector
	defiGauss3Pt g3;              //Initialize Gauss object

	//Data Holders
	defiVector x(3), y(3);        //x nodal positions, y nodal positions
	defiVector fx(3), fy(3);      //x forces, y forces
	fx.zero();
	fy.zero();

	x.setCoeff(0, d_node1->getX());         //Get x position data
	x.setCoeff(1, d_node2->getX());
	x.setCoeff(2, d_node3->getX());
	y.setCoeff(0, d_node1->getY());         //Get y position data
	y.setCoeff(1, d_node2->getY());
	y.setCoeff(2, d_node3->getY());

	//Equivalent Nodal Forces
	for (int i = 0; i < g3.getNumPts(); i++)
	{
		double z = g3.getGaussPt(i);        //Zeta
		double weight = g3.getWeight(i);    //Weight
		defiVector dndz(3);                 //Derivatives of shape factor w.r.t. zeta
		defiVector fx_temp(3), fy_temp(3);  //Force component holders

		//Compute Internal Sums
		double tempx, tempy;
		for (int j = 0; j < 3; j++)
		{
			//Fill Shape Factor Derivatives
			double dndz1, dndz2, dndz3;

			dndz1 = z - 0.5;                //Node 1
			dndz2 = -2.0*z;                 //Node 2
			dndz3 = z + 0.5;                //Node 3

			dndz.setCoeff(0, dndz1);
			dndz.setCoeff(1, dndz2);
			dndz.setCoeff(2, dndz3);

			//Compute Scalars
			tempx = (tn*globMultiply(dndz, y) + ts*globMultiply(dndz, x))*weight;
			tempy = (-1.0*tn*globMultiply(dndz, x) + ts*globMultiply(dndz, y))*weight;
		}

		//Fill Temporary Force Component Holders
		double n1, n2, n3;
		n1 = z*(z - 1.0) / 2.0;
		n2 = 1.0 - z*z;
		n3 = z*(z + 1.0) / 2.0;

		fx_temp.setCoeff(0, n1*tempx);
		fx_temp.setCoeff(1, n2*tempx);
		fx_temp.setCoeff(2, n3*tempx);
		fy_temp.setCoeff(0, n1*tempy);
		fy_temp.setCoeff(1, n2*tempy);
		fy_temp.setCoeff(2, n3*tempy);

		//Add to Force Holders
		globAdd(fx, fx_temp);
		globAdd(fy, fy_temp);
	}

	//Assemble the Force Vector
	for (int i = 0; i < fx.getNumRows(); i++)
	{
		f.setCoeff(2 * i, fx.getCoeff(i));
		f.setCoeff(2 * i + 1, fy.getCoeff(i));
	}

	//Store Nodal Dofs
	dofs[0] = d_node1->getDof(UX);                  //We should think about storing the nodes in the NaturalBCs as arrays so that we can use For loops to index through them...
	dofs[1] = d_node1->getDof(UY);                  //Because writing every single command out is boring and inefficient...
	dofs[2] = d_node2->getDof(UX);                  //Not to mention that it kills the versitility of this function since it can only work for Q8 elements since so much is hard coded...
	dofs[3] = d_node2->getDof(UY);
	dofs[4] = d_node3->getDof(UX);
	dofs[5] = d_node3->getDof(UY);	
}

#endif