/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "Matrix.h"

Matrix:: Matrix(int nr, int nc)
{
	//ROW_MAJOR layout;
	ci_matRow = nr;
	ci_matCol = nc;
	cdp2_matCoeff = new double[nr*nc];
	/*for (int i = 0; i < nr; i++)
	{
		cdp2_matCoeff[i] = new double[nc];
	}*/
}

Matrix::~ Matrix()
{
	   /*for (int i = 0; i < ci_matRow; i++)
	   {
		   delete[] cdp2_matCoeff[i];
		   cdp2_matCoeff[i] = NULL;
	   }*/
	   delete[] cdp2_matCoeff;
	   cdp2_matCoeff = NULL;
}

void Matrix::zero()
{
	for (int i = 0; i <ci_matRow; i++)
	{
		for (int j = 0; j <ci_matCol; j++)
		{
			//cdp2_matCoeff[i][j] = 0;
			cdp2_matCoeff[i * ci_matCol + j] = 0;
		}
	}
}

int Matrix::i_getNumRows() const
{
	return ci_matRow;
}

int Matrix::i_getNumCols() const
{
	return ci_matCol;
}

void  Matrix::print(ofstream & fout) const
{
	fout << setiosflags(ios::scientific)
		   << setprecision(12) ;
	fout << "Matrix Size: " << "(" << ci_matRow << ", " << ci_matCol << ")" << endl;
	for (int i = 0; i <ci_matRow; i++)
	{
		//fout << "Row(" << i << "):" << endl;
		for (int j = 0; j <ci_matCol; j++)
		{
			fout << setw(15) << cdp2_matCoeff[i * ci_matCol + j] << ' ';
		}
		fout << endl;
	}
}

void  Matrix::setCoeff(int i, int j, double value)
{
	cdp2_matCoeff[i * ci_matCol + j] = value;
}

void  Matrix::addCoeff(int i, int j, double value)
{
	cdp2_matCoeff[i * ci_matCol + j] = value + cdp2_matCoeff[i * ci_matCol + j];
}

double  Matrix::d_getCoeff(int i, int j)
{
	return cdp2_matCoeff[i * ci_matCol + j];
}
