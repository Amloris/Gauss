//-----------------------------------------------------------------------------
//  Gauss: A 2-D Finite Elements Program
//-----------------------------------------------------------------------------
//  Program Info
//		Creator: Aaron Robertson
//		Date   : April 2017
//
//  Program Notes:
//      The input files must be located in the data directory.
//      Additionally, the user must manually specify the input file name
//		in "globAcessItems.h"
//-----------------------------------------------------------------------------
//	Capabilities:
//		This program will take in 2D CASCA finite element data, construct the
//      global stiffness matrix, and solve for stress, strain, and displacement
//      that the mesh will experience due to forces and constraints.
//      It uses Q8 quadrilateral elements and linear elastic materials to 
//      evaluate the solution.
//	Future Additions:
//		3D Mesh Solutions
//		Mixed Element Meshes
//      Nonlinear Materials
//-----------------------------------------------------------------------------

#include <iostream>
#include <string>
#include "globAccessItems.h"
#include "dataFemModel.h"
#include "calcFemSolver.h"

#include "defiGauss3Pt.h"

using namespace std;

ifstream fin;		//Define global fin here (allocate memory space)
ofstream fout;		//Define global fout here (allocate memory space)
ofstream fplot;		//Define global fplot here (allocate memory space)
ofstream ferr;      //Define global ferr here (allocate memory space)

int main()
{
	dataFemModel o_modeldata;
	calcFemSolver o_sol;
	
	//Open Files
	globOpenFiles();
	cout << "FEM program started..." << endl;

	//Read Model Data
	o_modeldata.readData();

	//write model data for checking
	o_modeldata.writeData();
	
	//solve the fem model including displacements and stresses
	o_sol.solveFem(o_modeldata);

	//write results
	o_modeldata.writeResults();

	//write plot data
	o_modeldata.writePlotFile();


	cout << "FEM program finsihed." << endl;

	//Close Files
	globCloseFiles();

	getchar();          //Hold for user (so we can read the console printouts)

	return 0;
}

