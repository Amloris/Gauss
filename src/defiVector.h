/* defiVector
The defiVector class stores a vector (a 1D array).
*/

#ifndef defiVector_h
#define defiVector_h

#include <fstream>

using namespace std;

class defiVector {
private:
	double *m_coeff;		// Coefficients of the vector
	int m_nr;				// Number of rows
							// never to be used copy constructor
	defiVector(const defiVector &vector);
public:
	// Constructors
	defiVector();
	defiVector(int nr);
	//Destructor
	~defiVector();

	// Functions
	void zero();	                        //initialize all coefficient of the vector to be zero
	int getNumRows() const;	                //return m_nr;
	void setCoeff(int i, double value);	    //set m_coeff[i] = value;
	void addCoeff(int i, double value);	    //add value to m_coeff[i];
	double getCoeff(int i);	                //return m_coeff[i];
	void print();	                        //print useful info
};


//Functions
//-----------------------------------------------------------------------------
defiVector::defiVector()  { std::cout << "Creating defiVector Class Object" << endl; }
defiVector::~defiVector() { /*std::cout << "Deleting defiVector Class Object" << endl;*/ }
defiVector::defiVector(int nr)
{
	//Dynamically Allocate Space
	m_nr = nr;
	m_coeff = new double[nr];         
}

int defiVector::getNumRows() const
{
	return m_nr;
}

void defiVector::setCoeff(int i, double value)
{
	m_coeff[i] = value;
}

void defiVector::addCoeff(int i, double value)
{
	m_coeff[i] += value;
}

double defiVector::getCoeff(int i)
{
	return m_coeff[i];
}

void defiVector::zero()
{
	for (int i = 0; i < m_nr; i++)
	{
		m_coeff[i] = 0;
	}
}

void defiVector::print()
{
	cout << "Vector: length = " << m_nr << endl;
	for (int i = 0; i < m_nr; i++)
	{
		cout << m_coeff[i] << "    ";
	}
	cout << endl;
}

#endif