//-----------------------------------------------------------------------------
//  Regolith: A 2-D Finite Elements Program
//-----------------------------------------------------------------------------
//  Program Info
//		Creator: Aaron Robertson
//		Date   : April 2017
//		Class  : ME 862
//
//  Program Notes:
//      The input file must be located in the same directory as main.cpp.
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
#include "dataFemModel.h"
#include "calcFemSolver.h"
#include "globAccessItems.h"

#include "defiGauss3Pt.h"

using namespace std;

ifstream fin;		//define global fin here (allocate memory space)
ofstream fout;		//define global fout here (allocate memory space)
ofstream fplot;		//define global fplot here (allocate memory space)
ofstream ferr;      //define global ferr here (allocate memory space)

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

	/*
	//Work in Progress
	//-------------------------------------------------------------------------

	//write results
	o_modeldata.writeResults();

	//write plot data
	o_modeldata.writePlotFile();

	*/

	cout << "FEM program finsihed." << endl;

	//Close Files
	globCloseFiles();

	getchar();      //Hold for user (so we can read the console printouts)

	return 0;
}

