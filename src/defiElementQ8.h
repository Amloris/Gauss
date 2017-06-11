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

void defiElementQ8::calcStressesGaussPtsCorners()
{	//1. Calculates the stresses at the 3x3 gaussian points
	//2. Best fit the stresses at Gaussian points to the bilinear function,  sigma = a + b*ζ + c*η
	//3. Calculate the stresses at the corners by using the bilinear function

	//Initial Values
	int NumDofs = getNumDofs();           //Number of displacements
	defiMatrix D(3, 3);                   //Initialize matrices
	defiVector u(NumDofs);                //Initialize displacement vector
	defiGauss3Pt g3;                      //Initialize Gauss Object

	//Get Displacement Vector
	for (int i = 0; i < d_numNodes; i++)
	{
		u.setCoeff(2 * i, d_nodes[i]->getDof(UX)->getValue());       //x component of node displacement
		u.setCoeff(2 * i + 1, d_nodes[i]->getDof(UY)->getValue());   //y component of node displacement
	}
	
	//Get D Matrix
	d_material->getDMatrixStress(D);

	//Stresses at Gaussian Points
	defiMatrix stress_x(3, 3), stress_y(3, 3), stress_xy(3, 3);
	for (int i = 0; i < g3.getNumPts(); i++)
	{
		for (int j = 0; j < g3.getNumPts(); j++)
		{
			//Setup Variables
			double detj = 0.0;                             //Dummy Value
			double n = g3.getGaussPt(2-i);                 //Eta value
			double z = g3.getGaussPt(j);                   //Zeta value

			//Setup Matrices and Vectors
			defiMatrix B(3, 16);                           //B Matrix
			defiVector Strain(3), Stress(3);               //Stress and Strain

			//Matrix Formations
			getBMatrix(n, z, B, detj);                     //Get B Matrix
			globMultiply(B, u, Strain);                    //Get strain at gauss point
			globMultiply(D, Strain, Stress);               //Get stress at gauss point

			//Gaussian Stress Component Printouts
			/*
			cout << z << "  " << n << endl;
			Stress.print();
			cout << endl;
			*/

			//Store Data in d_stresssGPs
			d_stressGPs[i][j].setSigXX(Stress.getCoeff(0));     //Orientation of storage matches the orientation 
			d_stressGPs[i][j].setSigYY(Stress.getCoeff(1));     //of the isoparametric element in device unit space 
			d_stressGPs[i][j].setSigXY(Stress.getCoeff(2));

			//Store Data in Matrices
			stress_x.setCoeff(i, j, Stress.getCoeff(0));
			stress_y.setCoeff(i, j, Stress.getCoeff(1));
			stress_xy.setCoeff(i, j, Stress.getCoeff(2));
		}
	}
	//Separated Gaussian Stress Components Printout
	/*
	cout << endl << endl;
	stress_x.print();
	stress_y.print();
	stress_xy.print();
	*/


	//Bilinear Fit of Gaussian Stress Components (Failed attempt)
	//----------------------------------------------------------------------------
	//Conducts a bilinear fit of the stress compoenents by using matrix
	//operations. The independent values array stores the summations of the
	//dependent variables ζ and η.  The dependent value arrays store info relating
	//the dependent variable and the independent variables. The constants for the 
	//least squares fit are then solved for using Gauss Jordan elimination.
	//----------------------------------------------------------------------------
	//Warning:
	//     This section works most of the time, but if any of the stress
	//components are zero at the calculated Gaussian points the constants are
	//determined to be zero by the Gauss-Jordan solver.  This section (when it
	//decides to work) agrees with the constants that are computed using the 
	//second method listed below. Use this bilinear least sqaures fit calculator
	//with caution until the problem can be fixed.
	//----------------------------------------------------------------------------
	/*
	//Construct Independent Values Array
	defiMatrix A(3, 3);                                    //Holds coefficients of primary array for Least Sqaures fit later
	A.zero();
	A.setCoeff(0, 0, 9.0);                                 //Set N

	for (int i = 0; i < g3.getNumPts(); i++)
	{
		for (int j = 0; j < g3.getNumPts(); j++)
		{
			//Setup Variables
			double n = g3.getGaussPt(2 - i);               //Eta value
			double z = g3.getGaussPt(j);                   //Zeta value

			//Load the Matrix
			A.addCoeff(0, 1, z);
			A.addCoeff(0, 2, n);
			A.addCoeff(1, 0, z);
			A.addCoeff(1, 1, z*z);
			A.addCoeff(1, 2, z*n);
			A.addCoeff(2, 0, n);
			A.addCoeff(2, 1, z*n);
			A.addCoeff(2, 2, n*n);
		}
	}

	//Construct Dependent Value Arrays
	defiVector x(3), y(3), xy(3);                            //Stores the dependent component arrays 
	defiVector x_coeff(3), y_coeff(3), xy_coeff(3);          //Stores the best fit constants
	x.zero(); y.zero(); xy.zero();
	x_coeff.zero(); y_coeff.zero(); xy_coeff.zero();

	for (int i = 0; i < g3.getNumPts(); i++)
	{
		for (int j = 0; j < g3.getNumPts(); j++)
		{
			//Setup Variables
			double n = g3.getGaussPt(2 - i);               //Eta value
			double z = g3.getGaussPt(j);                   //Zeta value

			double x_val = d_stressGPs[i][j].getSigXX();
			double y_val = d_stressGPs[i][j].getSigYY();
			double xy_val = d_stressGPs[i][j].getSigXY();

			//Load the Matrix
			x.addCoeff(0, x_val);
			x.addCoeff(1, z*x_val);
			x.addCoeff(2, n*x_val);

			y.addCoeff(0, y_val);
			y.addCoeff(1, z*y_val);
			y.addCoeff(2, n*y_val);

			xy.addCoeff(0, xy_val);
			xy.addCoeff(1, z*xy_val);
			xy.addCoeff(2, n*xy_val);
		}
	}

	
	globGaussJordan(&A, &x, &x_coeff);
	globGaussJordan(&A, &y, &xy_coeff);
	globGaussJordan(&A, &xy, &xy_coeff);

	x_coeff.print();
	y_coeff.print();
	xy_coeff.print();
	*/
	


	//Bilinear Fit of Gaussian Stress Components
	//----------------------------------------------------------------------------
	//Conducts a bilinear fit of the stress components by using summations
	//involving the independent variables and dependent variable. Formulas are 
	//then used that compute the curve fit constants by using the summations for
	//a system with two indepentent variables.  It is messy and hardcoded, but it
	//works...
	//----------------------------------------------------------------------------

	//Bilinear Curve Fits
	double x_avg = 0.0;      //Sum for x stress average  
	double y_avg = 0.0;      //Sum for y stress average
	double xy_avg = 0.0;     //Sum for xy stress average

	double z_avg = 0.0;      //Sum for zeta average
	double n_avg = 0.0;      //Sum for eta average

	double z_x = 0.0;        //Sum for zeta*(x stress)    
	double z_y = 0.0;        //Sum for zeta*(y_stress)
	double z_xy = 0.0;       //Sum for zeta*(xy_stress)
	double n_x = 0.0;        //Sum for eta*(x stress) 
	double n_y = 0.0;		 //Sum for eta*(y_stress)
	double n_xy = 0.0;       //Sum for eta*(xy_stress)

	double z2 = 0.0;         //Shape factor location sum (zeta^2)
	double n2 = 0.0;         //Shape factor location sum (eta^2)
	double zn = 0.0;         //Shape factor location sum (zeta*eta)
	for (int i = 0; i < g3.getNumPts(); i++)
	{
		for (int j = 0; j < g3.getNumPts(); j++)
		{
			//Gauss Points
			double n = g3.getGaussPt(2-i);                   //Eta value
			double z = g3.getGaussPt(j);                     //Zeta value

			//Computed Stresses       
			double x = d_stressGPs[i][j].getSigXX();	     //Stress x at gaussian point
			double y = d_stressGPs[i][j].getSigYY();         //Stress y at gaussian point
			double xy = d_stressGPs[i][j].getSigXY();        //Stress xy at gaussian point

			//Add to Sums
			x_avg += x;
			y_avg += y;
			xy_avg += xy;
			z_avg += z;
			n_avg += n;

			z2 += z*z;
			n2 += n*n;
			zn += z*n;

			z_x += z*x;
			n_x += n*x;
			z_y += z*y;
			n_y += n*y;
			z_xy += z*xy;
			n_xy += n*xy;
		}
	}
	x_avg /= 9.0;            //Average x stress
	y_avg /= 9.0;            //Average y stress
	xy_avg /= 9.0;           //Average xy stress
	z_avg /= 9.0;            //Average zeta value
	n_avg /= 9.0;            //Average eta value

	//Least Squares Fit
	double a_x, b_x, c_x;         //Constants for x stress
	double a_y, b_y, c_y;         //Constants for y stress
	double a_xy, b_xy, c_xy;      //Constants for xy stress

	double denom;                 //The denominator value that all constant calculations share
	denom = z2*n2 - zn*zn;

	//Solve Constants for x stress Fit
	b_x = (n2*z_x - zn*n_x) / denom;
	c_x = (z2*n_x - zn*z_x) / denom;
	a_x = x_avg - b_x*z_avg - c_x*n_avg;

	//Solve Constants for y stress Fit
	b_y = (n2*z_y - zn*n_y) / denom;
	c_y = (z2*n_y - zn*z_y) / denom;
	a_y = y_avg - b_y*z_avg - c_y*n_avg;

	//Solve Constants for xy stress Fit
	b_xy = (n2*z_xy - zn*n_xy) / denom;
	c_xy = (z2*n_xy - zn*z_xy) / denom;
	a_xy = xy_avg - b_xy*z_avg - c_xy*n_avg;


	//Curvefit Constants Printout
	/*
	cout << a_x << " " << b_x << " " << c_x << endl;
	cout << a_y << " " << b_y << " " << c_y << endl;
	cout << a_xy << " " << b_xy << " " << c_xy << endl;
	*/

	//Calculate Stresses at Corners and Re-map the Locations
	//-------------------------------------------------------------------------
	//The section below will compute the stress components at the nodal
	//positions of the element.  It will re-map where the components are saved 
	//in d_stressGPs so that when getNodalData is called the returned vector
	//matches how the nodes are stored in the element. The stress components of
	//the center of the element are stored in last position of the returned 
	//vector.
	//-------------------------------------------------------------------------

	defiMatrix Pos(9, 2);                   //Stores the mapping positions of the nodes
	Pos.setCoeff(0, 0, -1.0);               //ζ position of first node
	Pos.setCoeff(0, 1, -1.0);               //η position of first node
	Pos.setCoeff(1, 0, 0.0);                //ζ position of second node
	Pos.setCoeff(1, 1, -1.0);               //η position of second node
	Pos.setCoeff(2, 0, 1.0);                //ζ position of third node
	Pos.setCoeff(2, 1, -1.0);               //η position of third node
	Pos.setCoeff(3, 0, 1.0);                             
	Pos.setCoeff(3, 1, 0.0);                            
	Pos.setCoeff(4, 0, 1.0);
	Pos.setCoeff(4, 1, 1.0);
	Pos.setCoeff(5, 0, 0.0);
	Pos.setCoeff(5, 1, 1.0);
	Pos.setCoeff(6, 0, -1.0);
	Pos.setCoeff(6, 1, 1.0);
	Pos.setCoeff(7, 0, -1.0);
	Pos.setCoeff(7, 1, 0.0);
	Pos.setCoeff(8, 0, 0.0);                //ζ position of center of element 
	Pos.setCoeff(8, 1, 0.0);                //η position of center of element

	for (int i = 0; i < Pos.getNumRows(); i++)
	{
		//Save Location
		int m = i / 3;                      //Computes where we will save the stresses in d_stressGPs
		int n = i % 3;                      //Computes where we will save the stresses in d_stressGPs

		//Natural Coordinate Values
		double zeta = Pos.getCoeff(i, 0);   //ζ
		double eta = Pos.getCoeff(i, 1);    //η

		//Compute Stresses at Remapped Positions
		double sigXX = 0.0;
		double sigYY = 0.0;
		double sigZZ = 0.0;
		double sigXY = 0.0;

		sigXX = a_x + b_x*zeta + c_x*eta;
		sigYY = a_y + b_y*zeta + c_y*eta;
		sigXY = a_xy + b_xy*zeta + c_xy*eta;

		//Compute SigZZ Based on Problem Type
		double nu = d_material->getNu();              //Poisson's Ratio
		int probType = d_material->getprobType();     //Problem type

		if (probType == 1)
		{
			sigZZ = 0.0;                              //Zero for plane stress
		}
		if (probType == 2)
		{
			sigZZ = nu*(sigXX + sigYY);               //Non-zero for plane strain
		}

		//Save to d_stressGPs
		d_stressGPs[m][n].setSigXX(sigXX);
		d_stressGPs[m][n].setSigYY(sigYY);
		d_stressGPs[m][n].setSigZZ(sigZZ);
		d_stressGPs[m][n].setSigXY(sigXY);


		//Stress Printout
		//cout << sigXX << "  " << sigYY << "  " << sigXY << "  " << sigZZ <<endl;
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