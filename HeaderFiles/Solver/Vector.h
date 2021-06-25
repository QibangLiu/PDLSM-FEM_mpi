/* Vector
The  Vector class defines a vector (a 1D arrau).
*/

#pragma once

//#include "stdafx.h"

#include <fstream>
#include<iomanip>
#include<cmath>
#include<iostream>
using namespace std;
class Vector {
public:
	double* cdp_vecCoeff;		// Coefficients of the vector
	int ci_vecRow;				// Number of rows
public:
	// Constructors
	Vector();
	Vector(int nr);
	//Destructor
	~Vector();

	//=======
	void setNumRows(int numRow);
	void setCoeff(double* vecCoeff);

	// Functions
	void zero();	// initialize all coefficient of the vector to be zero
	int i_getNumRows() const;	// return cd_vecRow;
	void setCoeff(int i, double value);	// set cdp_vecCoeff[i] = value;
	void addCoeff(int i, double value);	// add value to cdp_vecCoeff[i];
	double d_getCoeff(int i);	// return cdp_vecCoeff[i];
	double d_getmol()const;
	void print(ofstream &fout);	// print useful info
	void print();
	//void Vecmultiply(double value);//each component times value;
private:
	// never to be used copy constructor
	Vector(const Vector& o_vector);
};

