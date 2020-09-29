/*  Matrix
The  Matrix class defines a full matrix. 

It stores the numbers of rows and columns and all the coefficients of the matrix.
Matrix and vector multiplications, matrix tranpose and solution of linear sustem equations
using the Gaussian-Jordan elimination are performed using global functions described
in the globalAccessItems.h header file.
*/

#pragma once

//#include "stdafx.h"

#include "Vector.h"
#include<iomanip>
using namespace std;

class Matrix {	//this is a matrix class which can be unsummetric and full
private:
	double* cdp2_matCoeff;		// Coefficients of the matrix
	int ci_matRow, ci_matCol;        // Number of columns and rows

	Matrix();	//never to be used constructors
	//Matrix(const Matrix& o_matrix);
public:
	// Constructors
	Matrix(int nr, int nc);	//construct a nr by nc matrix
	//Destructor
	~Matrix();

	// Functions
	void zero();		//initialize all coefficients of the matrix to zero
	int i_getNumRows() const;	//	return ci_matRow;
	int i_getNumCols() const;	//	return ci_matCol;
	void print(ofstream &fout) const;	//print useful info;
	void setCoeff(int i, int j, double value);	//set cdp2_matCoeff[i][j] = value;
	void addCoeff(int i, int j, double value);	//add value to cdp2_matCoeff[i][j]
	double d_getCoeff(int i, int j);	//return cdp2_matCoeff[i][j];
};

