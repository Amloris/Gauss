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
	int d_id;					//element id
	defiMaterial* d_material;	//pointer to the material class
	int d_numNodes;				//total number of nodes of this element. For Q8, it is 8
	defiNode** d_nodes;			//array of ptrs to the (eight) nodes in this element

	defiElementQ8(); 	//never to be used constructors
	defiElementQ8(const defiElementQ8 &elem);
	//calcStress2D d_stressGPs[3][3];                                                     //Remember to re-enable this later!!!                                
	//d_stress is an 3 by 3 matrix, each element is a calcStress2D class object. 
	//d_stress stores the stresses at the 3 by 3 Gaussian points

public:
	// Constructors
	defiElementQ8(int id, defiMaterial* mat, int numNodes, defiNode* enodes[]);	    //Construct the element according to id, mat, numNodes and enodes[]
	// Destructor
	~defiElementQ8();

	// Functions
	int getID() const;						//return d_id;
	int getNumNodes() const;				//return d_numNodes;
	void getNodes(defiNode* nodes[]) const;	//assign each nodes[i] = d_nodes[i];
	defiMaterial* getMaterial() const;		//return d_material;
	void printData() const;		            //print any useful info
	int getNumDofs() const;	                //for Q8 element, return 16; 
	//void getElementK(defiMatrix &k, defiDof* dofs[]);
	void getElementK(defiMatrix &k, defiVector &x, defiVector &y);   //My modification of the stiffness matrix function (it needs the position coordinates supplied externally)
	//calculate the element stiffness matrix.
	//here k stores the 16 by 16 stiffness matrix, dofs[16] hold equation number for each dof for assembling element k to the global matrix 
	//calcFemSolver::assembleK will call this function to obtain k, and assemble the global matrix
	void calcStressesGaussPtsCorners();
	//this function calculates the stresses at Gaussian points, then constructs a bilinear best fit (a0+a1*zeta+a2*eta) to represent the stress distriution, then obtain the nodal stresses from the best fit
	void getNodalData(calcStress2D stresses[]); //This function gets nodal stresses in a vector form
	void printResults() const;	//print useful results such as stresses

private:
	void getBMatrix(double n, double z, defiVector &x, defiVector &y, defiMatrix &b, double &detj);
	/*This function calculate the B matrix.
	You need to evaluate the shape function derivatives with respect to r and s, the
	natural coordinates, and then global deriv of shape functions
	Calculate the determinant of the Jacobian and evaluate the B matrix.  Note that we write fout the inverse of the
	Jacobian matrix.
	*/
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

int defiElementQ8::getID() const       { return d_id; }              //return element id
int defiElementQ8::getNumNodes() const { return d_numNodes; }        //return number of nodes
int defiElementQ8::getNumDofs() const  { return d_numNodes * 2; }    //return degrees of freedom

defiMaterial* defiElementQ8::getMaterial() const
{
	return d_material;
}

void defiElementQ8::getNodes(defiNode* nodes[]) const
{	//assign each nodes[i] = d_nodes[i];

	int size = getNumNodes();
	for (int i = 0; i < size; i++)
	{
		nodes[i] = d_nodes[i];
	}
}

void defiElementQ8::printData() const
{	//print:  elemID    matID     NodeConnectivityData

	fout << "   " << d_id << "       " << d_material->getID() << "      ";
	for (int i = 0; i < d_numNodes; i++)
	{
		fout << "    " << d_nodes[i]->getID();
	}
	fout << endl;
}

void defiElementQ8::getBMatrix(double n, double z, defiVector &x, defiVector &y, defiMatrix &b, double &detj)
{	//Modified from the original to include the nodal positions.
	//x and y represent the global coordinates of the nodes.
	//eta=η=n, zeta=ζ=z.

	//Data Holders
	int num_nodes;
	num_nodes = getNumNodes();

	defiVector zeta(num_nodes);       //partial w.r.t. zeta   //partial N_i w.r.t. partial zeta
	defiVector eta(num_nodes);        //partial w.r.t. eta    //partial N_i w.r.t. partial eta

	defiMatrix J(2, 2);               //Jacobian
	defiMatrix Jinv(2, 2);            //Inverse Jacobian

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

	//zeta.print();
	//eta.print();

	//Solve for the Jacobian
	double x_z = globMultiply(zeta, x);     //partial x w.r.t. partial zeta
	double y_z = globMultiply(zeta, y);		//partial y w.r.t. partial zeta
	double x_n = globMultiply(eta, x);      //partial x w.r.t. partial eta
	double y_n = globMultiply(eta, y);      //partial y w.r.t. partial zeta

	J.setCoeff(0, 0, x_z);
	J.setCoeff(0, 1, y_z);
	J.setCoeff(1, 0, x_n);
	J.setCoeff(1, 1, y_n);

	//J.print();

	//Solve for Det of Jacobian                            //I should probably make a function to do this at some point
	detj = x_z*y_n - y_z*x_n;

	//Solve for Jacobain Inverse                           //I should probably make a function to do this at some point
	Jinv.setCoeff(0, 0, y_n / detj);
	Jinv.setCoeff(1, 1, x_z / detj);
	Jinv.setCoeff(0, 1, -1.0*y_z / detj);
	Jinv.setCoeff(1, 0, -1.0*x_n / detj);

	//Jinv.print();

	//Compute Components of B Matrix
	defiMatrix Partials(2, 8);                             //Holder for the partial derivatives of the shape functions evaluated at zeta and eta
	for (int i = 0; i < num_nodes; i++)
	{
		Partials.setCoeff(0, i, zeta.getCoeff(i));         //Load zeta partial components
		Partials.setCoeff(1, i, eta.getCoeff(i));          //Load eta partial components
	}
	defiMatrix B_Comp(2, 8);                               //Holder for the computed B matrix components
	globMultiply(Jinv, Partials, B_Comp);                  //Computing the B matrix components
	
	//B_Comp.print();

	//Construct the B Matrix
	b.zero();                                              //Set default values to zero
	for (int i = 0; i < num_nodes; i++)                   
	{
		int index = i * 2;                                 //The starting position in the b_matrix
		double N_x, N_y;
		N_x = B_Comp.getCoeff(0, i);
		N_y = B_Comp.getCoeff(1, i);

		//Load in Data
		b.setCoeff(0, index, N_x);
		b.setCoeff(2, index, N_y);
		b.setCoeff(1, index + 1, N_y);
		b.setCoeff(2, index + 1, N_x);
	}

	//b.print();
}

void defiElementQ8::getElementK(defiMatrix &k, defiVector &x, defiVector &y)
{	//Called for each element. Take an initialized 16x16 matrix "k" as well as
	//the node coordinate vectors.  Uses Gaussian integration to compute k.

	//Initial Values
	k.zero();                                              //Set stiffness matrix to zero
	double thickness = (*d_material).getThick();           //Material thickness
	defiMatrix D(3, 3), B(3, 16), Bt(16, 3);               //Initialize Matrices
	defiGauss3Pt g3;                                       //Initialize Gauss Object

	//Obtain D Matrix
	(*d_material).getDMatrixStress(D);                     //Only do this once since it is constant for one element

	//2D Gaussian Integration                              //This should be replaced by a function in the Gauss3Pt class at some point
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//Setup Variables
			double n = g3.getGaussPt(i);
			double z = g3.getGaussPt(j);
			double weight_i = g3.getWeight(i);
			double weight_j = g3.getWeight(j);
			double detj = 0.0;

			//Matrix Formations
			getBMatrix(n, z, x, y, B, detj);               //Get B Matrix
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
}

#endif