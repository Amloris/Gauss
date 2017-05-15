/* defiMatrix
The defiMatrix class is for a full matrix. It stores the numbers
of rows and columns and all the coefficients of the matrix.
Matrix and vector multiplications, matrix tranpose and solution of linear system equations
using the Gaussian-Jordan elimination are performed using global functions described
in the globalAccessItems.h header file.
*/

#ifndef defiMatrix_h
#define defiMatrix_h

using namespace std;

class defiMatrix {	//this is a matrix class which can be unsymmetric and full
private:
	double **m_coeff;		                     //Coefficients of the matrix
	int m_nr, m_nc;                              //Number of columns and rows

	defiMatrix();	                             //Never to be used constructors
	defiMatrix(const defiMatrix &matrix);
public:
	//Constructors
	defiMatrix(int nr, int nc);               	 //Construct a nr by nc matrix
	//Destructor
	~defiMatrix();

	// Functions
	void zero();		                         //initialize all coefficients of the matrix to zero
	int getNumRows() const;	                     //return m_nr;
	int getNumCols() const;	                     //return m_nc;
	void print() const;	                         //print useful info;
	void setCoeff(int i, int j, double value);	 //set m_coeff[i][j] = value;
	void addCoeff(int i, int j, double value);	 //add value to m_coeff[i][j]
	double getCoeff(int i, int j);	             //return m_coeff[i][j];
};

//Functions
//-----------------------------------------------------------------------------
defiMatrix::defiMatrix()  { std::cout << "Creating defiMatrix Class Object" << endl; }
defiMatrix::~defiMatrix()
{
	//Deallocate Memory
	for (int i = 0; i < m_nr; i++)
	{
		delete[] m_coeff[i];
	}
	delete[] m_coeff;

	//cout << "Deleting defiMatrix Class Object" << endl; 
}
defiMatrix::defiMatrix(int nr, int nc)
{
	//Load Data
	m_nr = nr;
	m_nc = nc;

	//Dynamically Allocate Memory
	m_coeff = new double *[nr];
	for (int i = 0; i < nr; i++)
	{
		m_coeff[i] = new double[nc];
	}
}

int defiMatrix::getNumRows() const
{
	return m_nr;
}

int defiMatrix::getNumCols() const
{
	return m_nc;
}

void defiMatrix::setCoeff(int i, int j, double value)	 
{
	m_coeff[i][j] = value;             //Set m_coeff[i][j] = value;
}

double defiMatrix::getCoeff(int i, int j)	             
{
	return m_coeff[i][j];              //Return m_coeff[i][j];
}

void defiMatrix::addCoeff(int i, int j, double value)	 
{
	m_coeff[i][j] += value;            //Add value to m_coeff[i][j]
}

void defiMatrix::zero()
{
	for (int i = 0; i < m_nr; i++)
	{
		for (int j = 0; j < m_nc; j++)
		{
			m_coeff[i][j] = 0;
		}
	}
}

void defiMatrix::print() const
{
	cout << "Matrix: " << m_nr << "x" << m_nc << endl;
	for (int i = 0; i < m_nr; i++)
	{
		for (int j = 0; j < m_nc; j++)
		{
			cout << m_coeff[i][j] << "    ";
		}
		cout << endl;
	}
	cout << endl;
}

#endif