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
{	//This section will solve the finite element system.
	//It sets eqn numbers based on constraints, builds the global stiffness
	//matrix and force vector, and solves for displacements using guassian elimination.
	//Finally, it computes the nodal stresses.


	//Get Number of System Equations
	m_neq = dat.getNumEq();                      //Total number of equations
	cout << "m_neq = " << m_neq << endl;

	//Set Contrained Dofs to Inactive
	defiEssentialBC* ebc;
	defiNode* node;
	enum DOFType dof;
	for (int i = 0; i < dat.getNumEssentialBCs(); i++)
	{
		ebc = dat.getEssentialBC(i);
		node = ebc->getNode();
		dof = ebc->getDof();
		node->getDof(dof)->setNotActive();                 //Set to inactive
		node->getDof(dof)->setValue(ebc->getValue());      //Set perscribed displacment
	}

	//Assign Equation Numbers
	int EqnNum = setEquationNumbers(dat);                  //Assign eqn numbers to active dofs
	if (EqnNum != m_neq) 
	{
		ferr << "ERROR::CALCFEMSOLVER_CLASS::NUM_SYSTEM_EQNS" << endl;
		exit(0);
	}

	//Dynamically Allocate Arrays (Global, Force, Displacement)
	m_k = new defiMatrix(m_neq, m_neq);
	m_f = new defiVector(m_neq);
	m_displ = new defiVector(m_neq);











	//Testing (Element Stiffness Matrix for first element)
	//-------------------------------------------------------------------------

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



int calcFemSolver::setEquationNumbers(dataFemModel &dat)
{	//Assigns an equation number to each dof that is active
	//Number of active dofs = (Num of nodes)*2 - (Num of inactive dofs)

	int EqnNum = 0;
	for (int i = 0; i < dat.getNumNodes(); i++)
	{
		if (dat.getNode(i)->getDof(UX)->isActive())
		{
			dat.getNode(i)->getDof(UX)->setEqn(EqnNum);
			EqnNum += 1;
		}
		if (dat.getNode(i)->getDof(UY)->isActive())
		{
			dat.getNode(i)->getDof(UY)->setEqn(EqnNum);
			EqnNum += 1;
		}
	}
	return EqnNum;
}

#endif