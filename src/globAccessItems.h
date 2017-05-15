/* globAcessItems.h
global variables and functions for use anywhere

In general, use of global variables and global functions are discouraged. However, judicious use of
global variables and global functions may be helpful sometimes. In this FEM program, setting the file streams
for the input, output and conture plot files to be global may be justifiable, because
1. file streams are for input/output only, and there is no need to make changes to these file streams
during the life of the program, reducing the risk of accidental changes.
2. making file streams global elimibates the need to passing them as function variables and
makes the program more concise. In a similar way that we use "cout" without the need to pass
the iostream "cout" as a function argument, we can declare a file stream "fout" as a global variable
so that we just use "fout" anywhere to output to the without passing it as a function argument.
3. being able to output information to the output file anywhere in any function is particularly
important during debugging. withour declaring "fout" a global argument, then we would have to pass
"fout" as a function argument for almost any function. this is too burdensome.

--- About using global variables ---
If the program contains only one source file (*.cpp file), to declare a global variable, you simply
put it outside of any function. If the program contains multiple source files (as in this FEM program),
then the global varible declaration in the header file must use the keyword "extern", and
the global variable must be defined only once in one source file.

To understand "extern", we need to first distinguish declaration and definition. To declare a variable
means to inform the compiler about the name and type of a variable, but not allocate memory space yet.
To define a variable means to allocation memory space according to the variable's type. In c and c++,
a variable can be declared multiple times, but can defined only once (the so-called One Definition Rule).

"extern" simply means "it is declared here, but it is to be defined (allocated momory space) somewhere
else in the program (externally to this file) - usually in another source file.

Example:
Declare an ofstream object "fout" as a global variable to be used in "main.cpp", "file1.cpp", and
"file2.cpp", with header file "header1.h"

In header file "header1.h":
...
extern ofstream fout;		//declare fout to be a global variable (object)
...

In source file "main.cpp":
#include "header1.h"
ofstream fout;				//define fout here
void main()
{
...
fout.open("fem.out");	//fout is linked to output file "fem.out"
...
}

In source file "file1.cpp":
#include "header1.h"
//ok to use "fout" because it is declared in hearder1.h and is defined in main.cpp
void func1()
{
...
fout << "ok to use fout to output without passing it as an argument!";
}

In source file "file2.cpp":
#include "header1.h"
//ok to use "fout" because it is declared in hearder1.h and is defined in main.cpp
void func2()
{
...
fout << "ok to use fout here too!";
}

In the above example, use of keyword "extern" in header1.h is necessary, because without it,
"fout" would have been defined 3 times (as main.cpp, file1.cpp, and file2.cpp all include
header1.h) which is not allowed in c++.

*/

#ifndef globAccessItems_h
#define globAccessItems_h

#include "defiMatrix.h"
#include "defiVector.h"

using namespace std;


extern ifstream fin;		//declare global fin - as the file stream for input file
extern ofstream fout;		//declare global fout - as the file stream for out file (result file)
extern ofstream fplot;		//declare global fplot - as the file stream for con file (contour plot file)
extern ofstream ferr;       //declare global ferr - as the file stream for error logging

void globOpenFiles()
{	//open input file, out file (results file), con file (contour plot file), and errorlog.

	//Open Files
	fin.open("sq.inp", ios::in);
	fout.open("fem.out", ios::out | ios::trunc);
	fplot.open("conplot.out", ios::out | ios::trunc);
	ferr.open("Error_Log.txt", ios::out | ios::trunc);

	//Check for Successful Open
	bool success = true;
	if (!ferr.is_open())                    { cout << "ERROR::ERRORLOG::FAILURE_TO_OPEN" << endl; success = false; }
	if (!fin.is_open() && ferr.is_open())   { ferr << "ERROR::FILE_INPUT::FAILURE_TO_OPEN" << endl; success = false; }
	if (!fout.is_open() && ferr.is_open())  { ferr << "ERROR::FILE_OUTPUT::FAILURE_TO_OPEN" << endl; success = false; }
	if (!fplot.is_open() && ferr.is_open()) { ferr << "ERROR::FILE_PLOT::FAILURE_TO_OPEN" << endl; success = false; }
	if (!success)                           { cout << "A problem was encountered opening one or more files" << endl; }
};

void globCloseFiles()
{	//Closes all files used by the program

	//Close Files
	fin.close();
	fout.close();
	fplot.close();
	ferr.close();

	//Check for Successful Close
	if (!fin.is_open() && !fout.is_open() && !fplot.is_open() && !ferr.is_open()) { cout << "All files closed successfully" << endl; }
	else { cout << "A problem was encountered closing one or more files" << endl; }
};


void globGaussJordan(defiMatrix *a, defiVector *b, defiVector *x);
// global function. it solves [a]{x}={b} using Gauss-Jordan elimination. 
// When it is done, {x} is the solution (=[a]^(-1) {b})
// and [a] becomes the inverse of the original [a]

void globTranspose(defiMatrix &a, defiMatrix &at)
{	//transpose matrix a and store the result in matrix at

	//Get Dimensions
	int rows, cols;
	rows = a.getNumRows();
	cols = a.getNumCols();

	//Transpose
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			at.setCoeff(j, i, a.getCoeff(i, j));
		}
	}
};

void globMultiply(defiMatrix &a, defiMatrix &b, defiMatrix &c)
{	//multiply matrix a with the matrix b, store the result in matrix c

	//Get Dimensions
	int a_rows, a_cols, b_rows, b_cols;
	a_rows = a.getNumRows();
	a_cols = a.getNumCols();
	b_rows = b.getNumRows();
	b_cols = b.getNumCols();

	//Multiply Arrays
	if (a_cols != b_rows)
	{
		ferr << "ERROR::MATRIX_MULT::DIMENSION_MISMATCH" << endl;
	}
	else
	{
		//Set to Zero
		c.zero();

		//Populate Arrays
		for (int i = 0; i < a_rows; i++)
		{
			for (int j = 0; j < b_cols; j++)
			{
				for (int k = 0; k < a_cols; k++)
				{
					double temp = a.getCoeff(i, k) * b.getCoeff(k, j);
					c.addCoeff(i, j, temp);
				}
			}
		}
	}
}

void globMultiply(defiMatrix &a, defiVector &x, defiVector &y)
{	//multiply matrix a with vector x, store the result in vector y

	//Get Dimensions
	int a_rows, a_cols, x_len;
	a_rows = a.getNumRows();
	a_cols = a.getNumCols();
	x_len = x.getNumRows();

	if (a_cols != x_len)
	{
		ferr << "ERROR::MATRIX_VEC_MULT::DIMENSION_MISMATCH" << endl;
	}
	else
	{
		//Set to Zero
		y.zero();

		//Populate Vector
		double sum;
		for (int i = 0; i < a_rows; i++)
		{
			sum = 0;
			for (int j = 0; j < a_cols; j++)
			{
				sum += a.getCoeff(i, j)*x.getCoeff(j);
			}
			y.setCoeff(i, sum);
		}
	}
}

void globMultiply(defiMatrix &a, double value)
{	//Multiply a matrix by a constant value and return it by the origonal name

	//Get Dimensions
	int rows, cols;
	rows = a.getNumRows();
	cols = a.getNumCols();

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			a.setCoeff(i, j, a.getCoeff(i, j)*value);
		}
	}
}

double globMultiply(defiVector &x, defiVector &y)
{	//Computes the dot product of two vectors

	//Get Dimensions
	int size_x = x.getNumRows();
	int size_y = y.getNumRows();

	if (size_x != size_y)
	{
		ferr << "ERROR::VEC_VEC_MULT::DIMENSION_MISMATCH" << endl;
		return 0;
	}
	else
	{
		//Perform Calculation
		double sum = 0;
		for (int i = 0; i < size_x; i++)
		{
			double tempx, tempy;
			tempx = x.getCoeff(i);
			tempy = y.getCoeff(i);

			sum += (tempx*tempy);
		}
		return sum;
	}
}

void globAdd(defiMatrix &a, defiMatrix &b)
{	//Add matrix b to a, and return a

	//Get Dimensions
	int a_rows, a_cols, b_rows, b_cols;
	a_rows = a.getNumRows();
	a_cols = a.getNumCols();
	b_rows = b.getNumRows();
	b_cols = b.getNumCols();

	//Add b to a
	for (int i = 0; i < a_rows; i++)
	{
		for (int j = 0; j < a_cols; j++)
		{
			a.setCoeff(i, j, a.getCoeff(i, j) + b.getCoeff(i, j));
		}
	}
}


#endif