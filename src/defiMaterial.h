/* defiMaterial
This class implements an elastic material.
*/

#ifndef defiMaterial_h
#define defiMaterial_h

#include "defiMatrix.h"
#include "globAccessItems.h"

using namespace std;

class defiMaterial
{
private:
	int d_id;		//material ID
	double d_e;		//Young's modulus
	double d_nu;	//Poisson's ratio
	double d_k;		//conductivity, irrelevant for us
	double d_thick; //thickness
	int d_probType; //problem type

public:
	// Constructors
	defiMaterial();
	defiMaterial(int id, double e, double nu, double k, double thick, int probType);
	//Destructor
	~defiMaterial();

	//Functions
	int getID() const;	          //return d_id
	double getE();                //return d_e
	double getNu();			      //return d_nu
	double getK();				  //return d_k
	double getThick();		      //return d_thick
	int getprobType();		      //return d_probType
	//void printData() const;		  //print useful info such as id, Young's modulus etc.
	void printData();		  //print useful info such as id, Young's modulus etc.
	void getDMatrixStress(defiMatrix &d);
	/*this sets the value of material stiffness matrix in d (3 by 3) or (4 by 4)
	For elastic plane stress or plane strain, the stresses are sigxx, sigyy, sigzz, sigxy,
	while the D matrix is a 4x4 matrix. One can also calculate sigxx, sigyy, sigxy
	using 3 by 3 D matrix, while calculate sigzz separately.
	*/
};


//Functions
//-----------------------------------------------------------------------------
defiMaterial::defiMaterial() { std::cout << "Creating defiMaterial Class Object" << endl; }
defiMaterial::~defiMaterial() { std::cout << "Deleting defiMaterial Class Object" << endl; }
defiMaterial::defiMaterial(int id, double e, double nu, double k, double thick, int probType) 
{
	d_id = id;
	d_e = e;
	d_nu = nu;
	d_thick = thick;
	d_k = k;
	d_probType = probType;
}

void defiMaterial::printData()
{
	//Get Data
	int id, probType;		  //material ID, problem type
	double e, nu, k, thick;   //Young's modulus, Poissons Ratio, Cond., thickness

	id = getID();
	e = getE();
	nu = getNu();
	k = getK();
	probType = getprobType();
	thick = getThick();

	//Print to File
	fout << "Elastic Material" << endl;
	fout << "Material ID     = " << id << endl;
	fout << "Young's Modulus = " << e << endl;
	fout << "Poisson's Ratio = " << nu << endl;
	fout << "Problem Type    = " << probType << endl;
	fout << "Thickness       = " << thick << endl;
	fout << endl;
}

void defiMaterial::getDMatrixStress(defiMatrix &d)
{	//Takes an initialized 3x3 matrix, d, and returns it as the Material
	//Stiffness Matrix.

	//Get Data 
	int prob_id;
	int row, col;
	double e, nu;
	prob_id = getprobType();      //Problem ID
	row = d.getNumRows();         //Get Matrix Dimensions
	col = d.getNumCols();         //Get Matrix Dimensions
	e = getE();                   //Get Modulus of Elasticity
	nu = getNu();                 //Get Poissons Ratio

	//Set to Zero
	d.zero();

	//Populate Array
	if (row != 3 && col != 3)
	{
		ferr << "ERROR::GET_D_MATRIX_STRESS::DIMENSIONS_INVALID" << endl;
	}
	else
	{
		if (prob_id == 1)
		{	//Plane Stress
			double temp = e / (1.0 - nu*nu);
			d.setCoeff(0, 0, temp);
			d.setCoeff(1, 1, temp);
			d.setCoeff(0, 1, nu*temp);
			d.setCoeff(1, 0, nu*temp);
			d.setCoeff(2, 2, temp*(1.0 - nu) / 2.0);
		}
		if (prob_id == 2)
		{
			//Plane Strain
			double temp = e / ((1.0 + nu)*(1.0 - 2.0*nu));
			d.setCoeff(0, 0, temp*(1.0 - nu));
			d.setCoeff(1, 1, temp*(1.0 - nu));
			d.setCoeff(0, 1, nu*temp);
			d.setCoeff(1, 0, nu*temp);
			d.setCoeff(2, 2, temp*(1.0 - 2*nu) / 2.0);
		}
	}
}

int defiMaterial::getID() const { return d_id; }
double defiMaterial::getE()     { return d_e; }
double defiMaterial::getNu()    { return d_nu; }
double defiMaterial::getK()     { return d_k; }
double defiMaterial::getThick() { return d_thick; }
int defiMaterial::getprobType() { return d_probType; }

#endif