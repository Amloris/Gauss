/* defiElementQ8
The defiElementQ8 class implements the eight-noded quadrilateral element.
The basic information contained in the defiElementQ8 class is an id,
a pointer to a material, the number of nodes in the element, and an
array of pointers to the nodes.  When the element is constructed,
all necessary data is passed to the constructor through the variable
argument list. It is possible to access the node id's, which is
needed when postprocessing.
*/

#ifndef defiElementQ8_h
#define defiElementQ8_h

#include "calcStress2D.h"

#include "defiGauss3Pt.h"

using namespace std;

class defiElementQ8
{
private:
	int d_id;					                 //Element id
	defiMaterial* d_material;	                 //Pointer to the material class
	int d_numNodes;				                 //Total number of nodes of this element. For Q8, it is 8
	defiNode** d_nodes;			                 //Array of ptrs to the (eight) nodes in this element

	defiElementQ8(); 	                         //Never to be used constructors
	defiElementQ8(const defiElementQ8 &elem);
	calcStress2D d_stressGPs[3][3];                                            
	     //d_stress is an 3 by 3 matrix, each element is a calcStress2D class object. 
	     //d_stress stores the stresses at the 3 by 3 Gaussian points

public:
	// Constructors
	defiElementQ8(int id, defiMaterial* mat, int numNodes, defiNode* enodes[]);	    //Construct the element according to id, mat, numNodes and enodes[]
	// Destructor
	~defiElementQ8();

	// Functions
	int getID() const;						     //Return d_id;
	int getNumNodes() const;				     //Return d_numNodes;
	void getNodes(defiNode* nodes[]) const;	     //Assign each nodes[i] = d_nodes[i];
	defiMaterial* getMaterial() const;		     //Return d_material;
	void printData() const;		                 //Print any useful info
	int getNumDofs() const;	                     //For Q8 element, return 16; 
	//void getElementK(defiMatrix &k, defiDof* dofs[]);
	void getElementK(defiMatrix &k, defiDof* dofs[]);             //My modification of the stiffness matrix function (Dofs[] not implemented)
	     //Calculates the element stiffness matrix.
	     //It stores the 16 by 16 stiffness matrix in "k", dofs[16] hold equation number for each dof for assembling element k to the global matrix 
	     //CalcFemSolver::assembleK will call this function to obtain k, and assemble the global matrix
	void calcStressesGaussPtsCorners();
	     //Calculates the stresses at Gaussian points, then constructs a bilinear best fit (a0+a1*zeta+a2*eta) to represent the stress distriution, then obtain the nodal stresses from the best fit
	void getNodalData(calcStress2D stresses[]);  //Gets nodal stresses in a vector form
	void printResults() const;	                 //Print useful results such as stresses

private:
	void getBMatrix(double n, double z, defiMatrix &b, double &detj);
	     //Calculates the B matrix.  It maps the arbitrary Q8 isoparametric element to unit normal space
		 //by using shape functions evaluated at 2D gaussian points and constructs the B matrix.
};


//Functions
//-----------------------------------------------------------------------------
defiElementQ8::defiElementQ8()  { std::cout << "Creating defiElementQ8 Class Object" << endl; }
defiElementQ8::~defiElementQ8() { std::cout << "Deleting defiElementQ8 Class Object" << endl; }
defiElementQ8::defiElementQ8(int id, defiMaterial* mat, int numNodes, defiNode* enodes[]) 
{
	d_id = id;
	d_material = mat;
	d_numNodes = numNodes;

	d_nodes = new defiNode*[numNodes];
	for (int i = 0; i < numNodes; i++)
	{
		d_nodes[i] = enodes[i];                  //Link to the 8 global nodes
	}
}

int defiElementQ8::getID() const       { return d_id; }              //Return element id
int defiElementQ8::getNumNodes() const { return d_numNodes; }        //Return number of nodes
int defiElementQ8::getNumDofs() const  { return d_numNodes * 2; }    //Return degrees of freedom

defiMaterial* defiElementQ8::getMaterial() const
{
	return d_material;
}

void defiElementQ8::getNodes(defiNode* nodes[]) const
{	//Assign each nodes[i] = d_nodes[i];

	int size = getNumNodes();
	for (int i = 0; i < size; i++)
	{
		nodes[i] = d_nodes[i];
	}
}

void defiElementQ8::printData() const
{	//Print:  elemID    matID     NodeConnectivityData

	fout << "   " << d_id << "       " << d_material->getID() << "      ";
	for (int i = 0; i < d_numNodes; i++)
	{
		fout << "    " << d_nodes[i]->getID();
	}
	fout << endl;
}

void defiElementQ8::getBMatrix(double n, double z, defiMatrix &b, double &detj)
{	//Calculates the B matrix by mapping the arbitrary element to unit coordinate space
	//eta=η=n, zeta=ζ=z.

	//Data Holders
	int num_nodes;
	num_nodes = getNumNodes();

	defiVector zeta(num_nodes);       //partial w.r.t. zeta   //partial N_i w.r.t. partial zeta
	defiVector eta(num_nodes);        //partial w.r.t. eta    //partial N_i w.r.t. partial eta
	defiVector dNdx(num_nodes);       //partial of N w.r.t partial of x
	defiVector dNdy(num_nodes);       //partial of N w.r.t partial of y

	//Compute Partials w.r.t. Inputs
	double z1, z2, z3, z4, z5, z6, z7, z8;
	double n1, n2, n3, n4, n5, n6, n7, n8;
	z1 = (n - 1.0)*(-1.0*n - 2.0*z) / 4.0;
	z2 = (n - 1.0)*(n - 2.0*z) / 4.0;
	z3 = (n + 1.0)*(n + 2.0*z) / 4.0;
	z4 = (n + 1.0)*(-1.0*n + 2.0*z) / 4.0;
	z5 = z*(n - 1.0);
	z6 = 0.5*(1.0 - n*n);
	z7 = -1.0*z*(1.0 + n);
	z8 = -0.5*(1.0 - n*n);
	n1 = -1.0*(1.0 - z)*(-2.0*n - z) / 4.0;
	n2 = (1.0 + z)*(2.0*n - z) / 4.0;
	n3 = (1.0 + z)*(2.0*n + z) / 4.0;
	n4 = -1.0*(1.0 - z)*(-2.0*n + z) / 4.0;
	n5 = -0.5*(1.0 - z*z);
	n6 = -1.0*n*(1.0 + z);
	n7 = 0.5*(1.0 - z*z);
	n8 = n*(z - 1.0);

	//Load into Holders
	/*
	//The shape function layout for the element described in the notes
	zeta.setCoeff(0, z1);
	zeta.setCoeff(1, z2);
	zeta.setCoeff(2, z3);
	zeta.setCoeff(3, z4);
	zeta.setCoeff(4, z5);
	zeta.setCoeff(5, z6);
	zeta.setCoeff(6, z7);
	zeta.setCoeff(7, z8);
	eta.setCoeff(0, n1);
	eta.setCoeff(1, n2);
	eta.setCoeff(2, n3);
	eta.setCoeff(3, n4);
	eta.setCoeff(4, n5);
	eta.setCoeff(5, n6);
	eta.setCoeff(6, n7);
	eta.setCoeff(7, n8);
	*/
	//The shape function layout for CASCA elements
	zeta.setCoeff(0, z1);
	zeta.setCoeff(1, z5);
	zeta.setCoeff(2, z2);
	zeta.setCoeff(3, z6);
	zeta.setCoeff(4, z3);
	zeta.setCoeff(5, z7);
	zeta.setCoeff(6, z4);
	zeta.setCoeff(7, z8);
	eta.setCoeff(0, n1);
	eta.setCoeff(1, n5);
	eta.setCoeff(2, n2);
	eta.setCoeff(3, n6);
	eta.setCoeff(4, n3);
	eta.setCoeff(5, n7);
	eta.setCoeff(6, n4);
	eta.setCoeff(7, n8);

	//Solve for the Jacobian
	double dxdz, dxdn, dydz, dydn;
	dxdz = 0.0; dxdn = 0.0; dydz = 0.0; dydn = 0.0;
	for (int i = 0; i < num_nodes; i++)
	{
		dxdz += zeta.getCoeff(i) * d_nodes[i]->getX();
		dxdn += eta.getCoeff(i) * d_nodes[i]->getX();
		dydz += zeta.getCoeff(i) * d_nodes[i]->getY();
		dydn += eta.getCoeff(i) * d_nodes[i]->getY();
	}

	//Solve for the Determinant of the Jacobian
	detj = dxdz*dydn - dxdn*dydz;
	
	if (fabs(detj) == 0.0)
	{
		ferr << "ERROR::BMATRIX_CALCULATION::ZERO_DETERMINANT" << endl;
		exit(0);
	}

	//Compute Components of B Matrix
	for (int i = 0; i < num_nodes; i++)
	{
		dNdx.setCoeff(i, (dydn*zeta.getCoeff(i) - dydz*eta.getCoeff(i)) / detj);
		dNdy.setCoeff(i, (-dxdn*zeta.getCoeff(i) + dxdz*eta.getCoeff(i)) / detj);
	}

	//Construct the B Matrix
	b.zero();                                              //Set default values to zero
	for (int i = 0; i < num_nodes; i++)
	{
		int index = i * 2;                                 //The starting position in the b_matrix

		//Load in Data
		b.setCoeff(0, index, dNdx.getCoeff(i));
		b.setCoeff(2, index, dNdy.getCoeff(i));
		b.setCoeff(1, index + 1, dNdy.getCoeff(i));
		b.setCoeff(2, index + 1, dNdx.getCoeff(i));
	}
}

void defiElementQ8::getElementK(defiMatrix &k, defiDof* dofs[])
{	//Computes the elemental stiffness matrix and stores it in 16x16 matrix "k".
	//Uses 2D Gaussian integration to compute k.

	//Initial Values
	k.zero();                                              //Set stiffness matrix to zero
	double thickness = d_material->getThick();             //Material thickness
	defiMatrix D(3, 3), B(3, 16), Bt(16, 3);               //Initialize Matrices
	defiGauss3Pt g3;                                       //Initialize Gauss Object

	//Obtain D Matrix
	d_material->getDMatrixStress(D);                       //Only do this once since it is constant for one element

	//2D Gaussian Integration                              //This should be replaced by a function in the Gauss3Pt class at some point
	for (int i = 0; i < g3.getNumPts(); i++)
	{
		for (int j = 0; j < g3.getNumPts(); j++)
		{
			//Setup Variables
			double n = g3.getGaussPt(i);                   //Eta value
			double z = g3.getGaussPt(j);                   //Zeta value
			double weight_i = g3.getWeight(i);
			double weight_j = g3.getWeight(j);
			double detj = 0.0;

			//Matrix Formations
			getBMatrix(n, z, B, detj);                     //Get B Matrix
			defiMatrix k_temp(16, 16), BtD(16, 3);         //Temporary Matrices
			globTranspose(B, Bt);                          //Transpose of the B Matrix
			globMultiply(Bt, D, BtD);                      //Multiply Bt and D together, store the result in BtD
			globMultiply(BtD, B, k_temp);                  //Multiply BtD and B together, store the result in k_temp

			double mult_factor = detj*weight_i*weight_j;   //The matrix multiplaction factor, including weights with the jacobian determinant
			mult_factor *= thickness;                      //Multiply by thickness
			globMultiply(k_temp, mult_factor);             //This is the k component matrix that will be added to the main k matrix

			//Add to Primary K Matrix
			globAdd(k, k_temp);
		}
	}

	//Store Element Dofs
	for (int i = 0; i < d_numNodes; i++)
	{
		dofs[2 * i] = d_nodes[i]->getDof(UX);
		dofs[2 * i + 1] = d_nodes[i]->getDof(UY);
	}
}

void defiElementQ8::getNodalData(calcStress2D stresses[])
{ //Gets nodal stresses in a vector form
	int p = -1;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			p += 1;
			stresses[p] = d_stressGPs[i][j];              //Retrieve gaussian stresses and store in vector form
		}
	}
}


#endif