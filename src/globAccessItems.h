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

#include <fstream>
#include "defiMatrix.h"
#include "defiVector.h"

using namespace std;


extern ifstream fin;		 //Declare global fin - as the file stream for input file
extern ofstream fout;		 //Declare global fout - as the file stream for out file (result file)
extern ofstream fplot;		 //Declare global fplot - as the file stream for con file (contour plot file)
extern ofstream ferr;        //Declare global ferr - as the file stream for error logging

void globOpenFiles()
{	//Open input file, output file, contour plot file, and the errorlog.

	//Open Files
	fin.open("../data/h60_pe_pt.inp", ios::in);
	fout.open("../data/fem.out", ios::out | ios::trunc);
	fplot.open("../data/conplot.out", ios::out | ios::trunc);
	ferr.open("../data/Error_Log.txt", ios::out | ios::trunc);

	//Check for Successful Open
	bool success = true;
	if (!ferr.is_open())                    { cout << "ERROR::ERRORLOG::FAILURE_TO_OPEN"    << endl; success = false; }
	if (!fin.is_open() && ferr.is_open())   { ferr << "ERROR::FILE_INPUT::FAILURE_TO_OPEN"  << endl; success = false; }
	if (!fout.is_open() && ferr.is_open())  { ferr << "ERROR::FILE_OUTPUT::FAILURE_TO_OPEN" << endl; success = false; }
	if (!fplot.is_open() && ferr.is_open()) { ferr << "ERROR::FILE_PLOT::FAILURE_TO_OPEN"   << endl; success = false; }
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


void globGaussJordan(defiMatrix *a, defiVector *b, defiVector *x)
{	// global function. it solves [a]{x}={b} using Gauss-Jordan elimination. 
	// When it is done, {x} is the solution (=[a]^(-1) {b})
	// and [a] becomes the inverse of the original [a]

	//#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
	int n, i, icol, irow, j, k, l, ll;
	double big, dum, pivinv, hold;
	double tmp1, tmp2;
	/*	dynamially create 3 integer arrays. The integer arrays
	p_ipiv, p_indxr, and p_indxc are used for bookkeeping on the pivoting.
	*/
	//n = m_nr;
	n = a->getNumRows();
	int* p_indxc = new int[n];
	int* p_indxr = new int[n];
	int* p_ipiv = new int[n];

	for (j = 0; j < n; j++) p_ipiv[j] = 0;
	for (i = 0; i < n; i++)
	{ // This is the main loop over the columns to be reduced. 
		big = 0.0;
		for (j = 0; j < n; j++)
		{// This is the outer loop of the search for a pivot element. 
			if (p_ipiv[j] != 1)
			{
				for (k = 0; k < n; k++)
				{
					if (p_ipiv[k] == 0)
					{
						if (fabs(a->getCoeff(j, k)) >= big)
						{
							big = fabs(a->getCoeff(j, k));
							irow = j;
							icol = k;
						}
					}
					else if (p_ipiv[k] > 1)
					{
						cout << "gaussj: Singular Matrix-1" << endl;
					}
				}
			}
		}
		++(p_ipiv[icol]);

		/* We now have the pivot element, so we interchange rows, if needed, to put the pivot
		element on the diagonal. The columns are not physically interchanged, only relabeled:
		p_indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
		p_indxr[i] is the row in which that pivot element was originally located. If p_indxr[i]
		6 = p_indxc[i] there is an implied column interchange. With this form of bookkeeping, the
		solution b's will end up in the correct order, and the inverse matrix will be scrambled
		by columns. */
		if (irow != icol)
		{
			for (l = 0; l < n; l++)
			{
				tmp1 = a->getCoeff(irow, l);
				tmp2 = a->getCoeff(icol, l);
				a->setCoeff(irow, l, tmp2);
				a->setCoeff(icol, l, tmp1);
				//SWAP(a.getCoeff[irow][l], a.getCoeff[icol][l]);
			}
			hold = b->getCoeff(irow);
			b->setCoeff(irow, b->getCoeff(icol));
			b->setCoeff(icol, hold);
		}

		// We are now ready to divide the pivot row by the pivot element, located at irow and icol.
		p_indxr[i] = irow;
		p_indxc[i] = icol;
		if (a->getCoeff(icol, icol) == 0.0)
		{
			cout << "gaussj: Singular Matrix-2" << endl;
		}
		pivinv = 1.0 / a->getCoeff(icol, icol);
		a->setCoeff(icol, icol, 1.0);
		for (l = 0; l < n; l++) { a->setCoeff(icol, l, a->getCoeff(icol, l)*pivinv); }
		b->setCoeff(icol, b->getCoeff(icol)*pivinv);
		for (ll = 0; ll < n; ll++) //Next, we reduce the rows... except for the pivot one, of course.
		{
			if (ll != icol)
			{
				dum = a->getCoeff(ll, icol);
				a->setCoeff(ll, icol, 0.0);
				for (l = 0; l < n; l++) a->setCoeff(ll, l, a->getCoeff(ll, l) - a->getCoeff(icol, l) * dum);
				b->addCoeff(ll, -dum * b->getCoeff(icol));
			}
		}
	}

	/* This is the end of the main loop over columns of the reduction. It only remains to unscram-
	ble the solution in view of the column interchanges. We do this by interchanging pairs of
	columns in the reverse order that the permutation was built up. */
	for (l = n - 1; l >= 0; l--)
	{
		if (p_indxr[l] != p_indxc[l])
		{
			for (k = 0; k < n; k++)
			{
				tmp1 = a->getCoeff(k, p_indxr[l]);
				tmp2 = a->getCoeff(k, p_indxc[l]);
				a->setCoeff(k, p_indxr[l], tmp2);
				a->setCoeff(k, p_indxr[l], tmp1);
				//				SWAP(m_coeff[k][p_indxr[l]], m_coeff[k][p_indxc[l]]);
			}
		}
	} // And we are done.
	delete[] p_indxc; //delete dynamically created arrays
	delete[] p_indxr;
	delete[] p_ipiv;

	for (i = 0; i < n; i++) {	//copy the solution vector to x object
		x->setCoeff(i, b->getCoeff(i));
	}
}

void globTranspose(defiMatrix &a, defiMatrix &at)
{	//Transpose matrix a and store the result in matrix at

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
{	//Multiply matrix a with the matrix b, store the result in matrix c

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
{	//Multiply matrix a with vector x, store the result in vector y

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

void globAdd(defiVector &x, defiVector &y)
{	//Adds vector y to x, and returns x

	//Get Dimensions
	int x_rows, y_rows;
	x_rows = x.getNumRows();
	y_rows = y.getNumRows();

	if (x_rows != y_rows)
	{
		ferr << "ERROR::VEC_ADD::DIMENSION_MISMATCH" << endl;
		exit(0);
	}

	//Add y to x
	for (int i = 0; i < x_rows; i++)
	{
		x.setCoeff(i, x.getCoeff(i) + y.getCoeff(i));
	}
}


#endif