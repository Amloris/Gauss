/* calcFemSolver
The calcFemSolver class manages the solution procedure.  The idea is to keep
the solution independent of the problem data and storage.

The solution procedure includes (1) finding the total number of equations
(=total number of node * 2 - total number of essential BCs), (2) assigining
equation number to each active dof, (3) assemblng the global stiffness
matrix and load vector, (4) solving the equation using the Gaussian-Jordan
elimination, (5) storing of the solution data (displacements) in the
d_value of the defiDof data member, and (6) calculating stresses and other
data required in postprocessing.
*/

#ifndef calcFemSolver_h
#define calcFemSolver_h

using namespace std;

class calcFemSolver
{
private:
	defiMatrix *m_k;	     //Ptr to a defiMatrix object. Use dynamic allocation to form a size of m_neq by m_neq global stiffness matrix    //e.g m_k = new defiMatrix(m_neq, m_neq); 
	defiVector *m_f;	     //Ptr to a defiVector object. It contains the right hand side vector. Use dynamic allocation to form a size of m_neq right hand side vector  //e.g. m_f = new defiVector(m_neq);
	defiVector *m_displ;	 //Ptr to a defiVector object. It contains the solution (solved displacement vector)
	int m_neq;		         //Total number of equations (=total number of nodes * 2 - total number of essential BCs).  Note that the index for eqns runs from 0 to m_neq-1

	int setEquationNumbers(dataFemModel &dat);
	/*this calculates the total number of equations
	which is equal to (number of nodes)*2 - (total number of inactive dofs) */
	void assembleRHS(dataFemModel &dat, defiVector *f);
	/*calculate contributions to the global right hand side vector from
	point force boundary condition and natural boundary condition.
	Contribution from the essential boundary will be calculated after element K matrix is obtained.
	*/
	void assembleK(dataFemModel &dat, defiMatrix *k, defiVector *f);
	/* assemble the global stiffness matrix. Only active dofs contribute to m_k
	*/
	void saveDislpacements(dataFemModel &dat, defiVector *sol);
	/* save the solved displacements (m_displ) into the d_value of each active dof
	*/
	void calcNodalStresses(dataFemModel &dat);
	//calculate nodal stresses
public:
	// Constructors
	calcFemSolver();
	//Destructor
	~calcFemSolver();
	void solveFem(dataFemModel &dat);
	/*
	this function solves the system equations (displacements) and then stresses.
	it first sets the dofs with essential b.c. to inactive, calculate the
	number of equations, then allocate the mamery for the global stiffness matrix and right hand vector,
	assemble RHS vector and global stiffness matrix, then solve the system equation using Gauss elimination
	It then saves the solution (displacements) and calculate the nodal stresses
	*/
};


//Functions
//-----------------------------------------------------------------------------
calcFemSolver::calcFemSolver()  { std::cout << "Creating calcFemSolver Class Object" << endl; }
calcFemSolver::~calcFemSolver() { std::cout << "Deleting calcFemSolver Class Object" << endl; }

void calcFemSolver::solveFem(dataFemModel &dat)
{	//This section will solve the finite element system

	//Testing (Element Stiffness Matrix for first element)
	//-------------------------------------------------------------------------
	//Note: This section is a placeholder to demonstrate that the element 
	//		stiffness matrrix can be calculated.  It will later be streamlined
	//		by using more functions and classes.

	//Get the K Matrix
	defiMatrix k(16, 16);                                  //Empty k matrix
	dat.getElem(0)->getElementK(k);                        //Compute k matrix

	//Print the K Matrix                                   //We have to print it to the file here since my Matrix and Vector class cant print to file
	int row, col;
	row = k.getNumRows();
	col = k.getNumCols();
	fout << endl << "***Stiffness Matrix (Element 1)***" << endl;
	for (int i = 0; i < row; i++)
	{
		fout << "Row " << i << endl ;
		fout.precision(4);
		fout << scientific;
		for (int j = 0; j < col; j++)
		{
			double temp = k.getCoeff(i, j);
			fout << temp << " ";
		}
		fout << endl << endl;
	}

	

}

#endif