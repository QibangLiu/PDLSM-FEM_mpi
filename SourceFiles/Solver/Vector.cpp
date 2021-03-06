/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "Vector.h"

 Vector:: Vector()
{
	 cdp_vecCoeff = NULL;
}

Vector:: Vector(int nr)
{
	ci_vecRow = nr;
	cdp_vecCoeff = new double[nr];
}

Vector::~ Vector()
{
	if (cdp_vecCoeff)
	{
		delete[]cdp_vecCoeff;
		cdp_vecCoeff = NULL;
	}
}

void Vector::setNumRows(int numRow)
{
	ci_vecRow = numRow;
}

void Vector::setCoeff(double* vecCoeff)
{
	cdp_vecCoeff = vecCoeff;
}

void Vector::zero()
{
	for (int i = 0; i <ci_vecRow; i++)
	{
		cdp_vecCoeff[i] = 0;
	}
}

int  Vector::i_getNumRows() const
{
	return ci_vecRow;
}

void Vector::setCoeff(int i, double value)
{
	cdp_vecCoeff[i] = value;
}

void  Vector::addCoeff(int i, double value)
{
	cdp_vecCoeff[i] = value + cdp_vecCoeff[i];
}

double  Vector::d_getCoeff(int i)
{
	return cdp_vecCoeff[i];
}

double Vector::d_getmol() const
{
	double sum=0;
	for (int i = 0; i < ci_vecRow; i++)
	{
		sum = sum + (cdp_vecCoeff[i])*(cdp_vecCoeff[i]);
	}
	return sqrt(sum);
}

void  Vector::print(ofstream & fout)
{
	fout << "Row Number: " << ci_vecRow << endl;
	for (int i = 0; i < ci_vecRow; i++)
	{
		fout <<setiosflags(ios::scientific)<< setprecision(12) << setw(15) << cdp_vecCoeff[i] << endl;
	}
}

void Vector::print()
{
	std::cout << "Row Number: " << ci_vecRow << endl;
	for (int i = 0; i < ci_vecRow; i++)
	{
		cout << setiosflags(ios::scientific) << setprecision(12) << setw(15) << cdp_vecCoeff[i] << endl;
	}
}

  /* void Vector::Vecmultiply(double value)
   {
	   for (int i = 0; i < ci_vecRow; i++)
	   {
		   cdp_vecCoeff[i] = cdp_vecCoeff[i] * value;
	   }
   }*/
