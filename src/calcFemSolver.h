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
		//this calculates the total number of equations
		//which is equal to (number of nodes)*2 - (total number of inactive dofs)
	void assembleK(dataFemModel &dat, defiMatrix *k, defiVector *f);
		//assemble the global stiffness matrix. Only active dofs contribute to m_k
	void assembleRHS(dataFemModel &dat, defiVector *f);
		//calculate contributions to the global right hand side vector from
		//point force boundary condition and natural boundary condition.
		//Contribution from the essential boundary will be calculated after element K matrix is obtained.
	void saveDislpacements(dataFemModel &dat, defiVector *sol);
		//save the solved displacements (m_displ) into the d_value of each active dof
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

	//Set Arrays to Zero
	m_k->zero();
	m_f->zero();
	m_displ->zero();

	//Assemble Global Stiffness Matrix
	assembleK(dat, m_k, m_f);

	//Assemble RHS Force Vector
	assembleRHS(dat, m_f);

	//Solve for Displacements	
	globGaussJordan(m_k, m_f, m_displ);

	//Save Displacements
	saveDislpacements(dat, m_displ);

	//Compute Nodal Stresses
	calcNodalStresses(dat);
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

void calcFemSolver::assembleK(dataFemModel &dat, defiMatrix *k, defiVector *f)
{	//Assembles the global stiffness matrix based on active Dofs.
	//Computes forces from perscribed displacements (EssentialBCs).

	//Construct Global Stiffness Matrix
	for (int i = 0; i < dat.getNumElems(); i++)
	{
		//Setup Storage
		int NumDofs;
		NumDofs = dat.getElem(i)->getNumDofs();       //Number of dofs in element

		defiMatrix k_elem(NumDofs, NumDofs);          //Element stiffness storage

		defiDof** eDofs;                              //Element Dof storage
		eDofs = new defiDof*[NumDofs];
		for (int j = 0; j < NumDofs; j++)
		{
			eDofs[j] = new defiDof();
		}

		//Element Stiffness Matrix
		dat.getElem(i)->getElementK(k_elem, eDofs);

		//Print First Element Stiffness Matrix
		if (i == 0)
		{
			int row, col;
			row = k_elem.getNumRows();
			col = k_elem.getNumCols();
			fout << endl << "***Stiffness Matrix (Element 1)***" << endl;
			for (int k = 0; k < row; k++)
			{
				fout << "Row " << k << endl;
				for (int j = 0; j < col; j++)
				{
					double temp = k_elem.getCoeff(k, j);

					fout.setf(ios::scientific);
					fout.precision(4);
					fout.width(13);
					fout << temp;
					if ((j + 1) % 6 == 0)
					{
						fout << endl;
					}
				}
				fout << endl << endl;
			}
		}

		//Assemble into Global Stiffness Matrix
		for (int m = 0; m < NumDofs; m++)
		{
			for (int n = 0; n < NumDofs; n++)
			{
				if (eDofs[m]->isActive() && eDofs[n]->isActive())
				{
					//Get k_elem Coefficient Value
					double temp_val;
					temp_val = k_elem.getCoeff(m, n);

					//Add it to Global Stiffness Matrix
					int row, col;
					row = eDofs[m]->getEqn();
					col = eDofs[n]->getEqn();
					k->addCoeff(row, col, temp_val);
				}
			}
		}

		//Compute Forces from Prescribed Displacments
		defiVector u(NumDofs);                             //Prescribed displacment vector
		defiVector force_ebc(NumDofs);                     //Force vector from EssentialBC
		u.zero();
		for (int q = 0; q < NumDofs; q++)
		{
			if (!eDofs[q]->isActive())                  
			{
				u.setCoeff(q, eDofs[q]->getValue());       //If constrained, load the prescibed displacement
			}
		}
		globMultiply(k_elem, -1.0);                        //Set element stiffness matrix to negative
		globMultiply(k_elem, u, force_ebc);
		//force_ebc.print();                               //Print forces from EssentialBCs

		//Assemble Global Load Vector (Perscribed Displacement)
		for (int p = 0; p < NumDofs; p++)
		{
			if (eDofs[p]->isActive())
			{
				//Get force_ebc Coefficient Value
				double temp_val;
				temp_val = force_ebc.getCoeff(p);

				//Add it to Global Force Vector
				int row;
				row = eDofs[p]->getEqn();
				f->addCoeff(row, temp_val);							
			}
		}
	}

	//Print Portion of Global Stiffness Matrix
	int row, col;
	row = k->getNumRows();
	col = k->getNumCols();
	fout << endl << "***Global Stiffness Matrix (" << row << "x" << col << ")***   Max Displayed (20x20)" << endl;
	
	if (col > 20) { col = 20; }        //Set printout limits
	if (row > 20) { row = 20; }        //Set printout limits

	for (int k = 0; k < row; k++)
	{
		fout << "Row " << k << endl;
		for (int j = 0; j < col; j++)
		{
			double temp = m_k->getCoeff(k, j);

			fout.setf(ios::scientific);
			fout.precision(4);
			fout.width(13);
			fout << temp;
			if ((j + 1) % 6 == 0)
			{
				fout << endl;
			}
		}
		fout << endl;
	}
}

void calcFemSolver::assembleRHS(dataFemModel &dat, defiVector *f)
{	//Computes equivalent nodal forces from NaturalBCs (tractions)
	//and adds them to the Global force vector along with Point forces
	//for all active dofs.

	m_neq = dat.getNumEq();                             //Total number of equations
	defiVector force_nbc(m_neq);                        //Holder for natural boundary condition forces
	defiVector force_pbc(m_neq);                        //Holder for point boundary condition forces

	//Compute Forces from NatualBCs
	force_nbc.zero();
	for (int i = 0; i < dat.getNumNaturalBCs(); i++)
	{
		//Setup Storage
		int NumDofs = 6;                                //Number of dofs for a naturalBC with 3 nodes

		defiVector f_nbc(NumDofs);                      //Storage for computed forces

		defiDof** eDofs;                                //Node Dof storage
		eDofs = new defiDof*[NumDofs];
		for (int j = 0; j < NumDofs; j++)
		{
			eDofs[j] = new defiDof();
		}

		//Calculate Equivalent Forces
		dat.getNaturalBC(i)->getNaturalBCF(f_nbc, eDofs);

		//Assemble into Force_nbc
		for (int j = 0; j < NumDofs; j++)
		{
			if (eDofs[j]->isActive())
			{
				//Get force_nbc Coefficient Value
				double temp_val;
				temp_val = f_nbc.getCoeff(j);

				//Add it to Global Force Vector
				int row;
				row = eDofs[j]->getEqn();

				force_nbc.addCoeff(row, temp_val);
			}
		}	
	}

	//Compute Forces from PointBCs
	force_pbc.zero();
	for (int i = 0; i < dat.getNumPointBCs(); i++)
	{
		//Get Direction
		PointBCType dof;
		dof = dat.getPointBC(i)->getPbctype();         //dof = 0 if FX, dof = 1 if FY
		DOFType index;

		if (dof == 0) { index = UX; }                  //Set to UX
		else          { index = UY; }                  //Set to UY

		//Assemble into Force_pbc
		if (dat.getPointBC(i)->getNode()->getDof(index)->isActive())
		{
			//Get Point Force Value
			double val;
			val = dat.getPointBC(i)->getValue();

			//Add it to Global Force Vector
			int row;
			row = dat.getPointBC(i)->getNode()->getDof(index)->getEqn();
			force_pbc.addCoeff(row, val);
		}
	}

	//Assemble into Global Force Vector
	globAdd(*f, force_nbc);
	globAdd(*f, force_pbc);
}

void calcFemSolver::saveDislpacements(dataFemModel &dat, defiVector *sol)
{	//Loads the solved displacements(m_displ) into the d_value of each active dof

	//Save Displacements
	for (int i = 0; i < dat.getNumNodes(); i++)
	{
		if (dat.getNode(i)->getDof(UX)->isActive())
		{
			//Get Eqaution Number
			int eqnNum = dat.getNode(i)->getDof(UX)->getEqn();

			//Get Value
			double val = sol->getCoeff(eqnNum);

			//Set Value
			dat.getNode(i)->getDof(UX)->setValue(val);
		}
		if (dat.getNode(i)->getDof(UY)->isActive())
		{
			//Get Eqaution Number
			int eqnNum = dat.getNode(i)->getDof(UY)->getEqn();

			//Get Value
			double val = sol->getCoeff(eqnNum);

			//Set Value
			dat.getNode(i)->getDof(UY)->setValue(val);
		}
	}
}

void calcFemSolver::calcNodalStresses(dataFemModel &dat)
{	//1. calculate stresses at Gaussian points
	//2. calculate stresses at each node by averaging all stress values 
	//   at the node over connecting elements

	//Gaussian Stresses at Corners for Each Element
	for (int i = 0; i < dat.getNumElems(); i++)
	{
		dat.getElem(i)->calcStressesGaussPtsCorners();
	}


	//Sum Nodal Stresses for Connected Elements
	int numNodes = dat.getNumNodes();               
	defiVector tally(numNodes);                     //Holds instances of found nodes
	tally.zero();

	for (int i = 0; i < dat.getNumElems(); i++)
	{
		//Data Holders
		int const eNumNodes = 8;                    //Number of nodes in q8 element
		defiNode* eNodes[eNumNodes];                //Holder for nodes
		dat.getElem(i)->getNodes(eNodes);           //Get list of nodes

		calcStress2D stressGPs[9];                  //Holder for stresses at nodal positions
		dat.getElem(i)->getNodalData(stressGPs);    //Retrieve stress data for element

		for (int j = 0; j < dat.getElem(i)->getNumNodes(); j++)
		{
			//Add to Tally
			int id = eNodes[j]->getID();
			tally.addCoeff(id-1, 1.0);

			//Get Stresses
			double sigxx = stressGPs[j].getSigXX();
			double sigyy = stressGPs[j].getSigYY();
			double sigzz = stressGPs[j].getSigZZ();
			double sigxy = stressGPs[j].getSigXY();

			//Update Stresses
			dat.getNode(id - 1)->getStress()->addSigXX(sigxx);
			dat.getNode(id - 1)->getStress()->addSigYY(sigyy);
			dat.getNode(id - 1)->getStress()->addSigZZ(sigzz);
			dat.getNode(id - 1)->getStress()->addSigXY(sigxy);
		}
	}

	//Average Nodal Stresses
	for (int i = 0; i < numNodes; i++)
	{
		//Get Stresses
		double sigxx = dat.getNode(i)->getStress()->getSigXX();
		double sigyy = dat.getNode(i)->getStress()->getSigYY();
		double sigzz = dat.getNode(i)->getStress()->getSigZZ();
		double sigxy = dat.getNode(i)->getStress()->getSigXY();

		//Average
		int n = tally.getCoeff(i);
		dat.getNode(i)->getStress()->setSigXX(sigxx / n);
		dat.getNode(i)->getStress()->setSigYY(sigyy / n);
		dat.getNode(i)->getStress()->setSigZZ(sigzz / n);
		dat.getNode(i)->getStress()->setSigXY(sigxy / n);
	}
}
#endif