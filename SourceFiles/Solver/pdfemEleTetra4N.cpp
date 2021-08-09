/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdfemEleTetra4N.h"
extern calMatrixOperations matoperat;
pdfemEleTetra4N::pdfemEleTetra4N(int id, int* nId, int algoType, dataLev2* p_datLev2):
	pdfemEles(id, 4, nId, algoType, p_datLev2)
{
	ci_eleType = 10;
}

pdfemEleTetra4N::~pdfemEleTetra4N()
{
}

void pdfemEleTetra4N::eleCenter(double xc[], double xN[][3])
{
	xc[0] = 0; xc[1] = 0; xc[2] = 0;
	double N[4] = { 0.25,0.25,0.25,0.25 };
	for (int n = 0; n < ci_numNodes; n++)
	{
		for (int i = 0; i < 3; i++)
		{
			xc[i] = xc[i] + xN[n][i] * N[n];
		}
	}
}

void pdfemEleTetra4N::eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3])
{
	Matrix* B = new Matrix(6, 12);
	Matrix* Bt, * BtD;
	Bt = new Matrix(12, 6);
	BtD = new Matrix(12, 6);
	double V = eleVolume(xN);
	BmatFEM(B, xN);
	matoperat.matTranspose(B, Bt);
	matoperat.matMultiply(Bt, D, BtD);
	matoperat.matMultiply(BtD, B, Ke);
	matoperat.matMultiply(Ke, V, Ke);
	delete B, delete Bt, delete BtD;
	B = NULL, Bt = NULL, BtD = NULL;

}

void pdfemEleTetra4N::eleMassMat(Matrix* Me, double rho, double xN[][3])
{
	double ve, fac;
	ve = eleVolume(xN);
	fac = rho * ve / 20.0;

	Me->zero();
	for (int i = 0; i < 12; i++)
	{
		Me->setCoeff(i, i, fac * 2);
	}
	int col;
	for (int row = 0; row < 12; row++)
	{
		col = row + 3;
		while (col<12)
		{
			Me->setCoeff(row, col, fac);
			Me->setCoeff(col, row, fac);
			col = col + 3;
		}
		
	}
}

void pdfemEleTetra4N::eleEquivNodalForce(Vector* Fe, double t, double xN[][3])
{
}

void pdfemEleTetra4N::BmatFEM(Matrix* B, double xN[][3])
{
	double a[4], b[4], c[4], d[4], mat[3][3];
	////==a0
	//mat[0][0] = xN[1][0], mat[0][1] = xN[1][1], mat[0][2] = xN[1][2];
	//mat[1][0] = xN[2][0], mat[1][1] = xN[2][1], mat[1][2] = xN[2][2];
	//mat[2][0] = xN[3][0], mat[2][1] = xN[3][1], mat[2][2] = xN[3][2];
	//a[0] = detMat33(mat);
	////a1
	//mat[0][0] = xN[0][0], mat[0][1] = xN[0][1], mat[0][2] = xN[0][2];
	//mat[1][0] = xN[2][0], mat[1][1] = xN[2][1], mat[1][2] = xN[2][2];
	//mat[2][0] = xN[3][0], mat[2][1] = xN[3][1], mat[2][2] = xN[3][2];
	//a[1] =- detMat33(mat);
	////a2
	//mat[0][0] = xN[0][0], mat[0][1] = xN[0][1], mat[0][2] = xN[0][2];
	//mat[1][0] = xN[1][0], mat[1][1] = xN[1][1], mat[1][2] = xN[1][2];
	//mat[2][0] = xN[3][0], mat[2][1] = xN[3][1], mat[2][2] = xN[3][2];
	//a[2] = detMat33(mat);
	////a3
	//mat[0][0] = xN[0][0], mat[0][1] = xN[0][1], mat[0][2] = xN[0][2];
	//mat[1][0] = xN[1][0], mat[1][1] = xN[1][1], mat[1][2] = xN[1][2];
	//mat[2][0] = xN[2][0], mat[2][1] = xN[2][1], mat[2][2] = xN[2][2];
	//a[3] = -detMat33(mat);
	//=====b
	//b0
	mat[0][0] = 1, mat[0][1] = xN[1][1], mat[0][2] = xN[1][2];
	mat[1][0] = 1, mat[1][1] = xN[2][1], mat[1][2] = xN[2][2];
	mat[2][0] = 1, mat[2][1] = xN[3][1], mat[2][2] = xN[3][2];
	b[0] = -detMat33(mat);
	//b1
	mat[0][0] = 1, mat[0][1] = xN[0][1], mat[0][2] = xN[0][2];
	mat[1][0] = 1, mat[1][1] = xN[2][1], mat[1][2] = xN[2][2];
	mat[2][0] = 1, mat[2][1] = xN[3][1], mat[2][2] = xN[3][2];
	b[1] = detMat33(mat);
	//b2
	mat[0][0] = 1, mat[0][1] = xN[0][1], mat[0][2] = xN[0][2];
	mat[1][0] = 1, mat[1][1] = xN[1][1], mat[1][2] = xN[1][2];
	mat[2][0] = 1, mat[2][1] = xN[3][1], mat[2][2] = xN[3][2];
	b[2] = -detMat33(mat);
	//b3
	mat[0][0] = 1, mat[0][1] = xN[0][1], mat[0][2] = xN[0][2];
	mat[1][0] = 1, mat[1][1] = xN[1][1], mat[1][2] = xN[1][2];
	mat[2][0] = 1, mat[2][1] = xN[2][1], mat[2][2] = xN[2][2];
	b[3] = detMat33(mat);
	//=====c
	//c0
	mat[0][0] = 1, mat[0][1] = xN[1][0], mat[0][2] = xN[1][2];
	mat[1][0] = 1, mat[1][1] = xN[2][0], mat[1][2] = xN[2][2];
	mat[2][0] = 1, mat[2][1] = xN[3][0], mat[2][2] = xN[3][2];
	c[0] = detMat33(mat);
	//c1
	mat[0][0] = 1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][2];
	mat[1][0] = 1, mat[1][1] = xN[2][0], mat[1][2] = xN[2][2];
	mat[2][0] = 1, mat[2][1] = xN[3][0], mat[2][2] = xN[3][2];
	c[1] =-detMat33(mat);
	//c2
	mat[0][0] = 1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][2];
	mat[1][0] = 1, mat[1][1] = xN[1][0], mat[1][2] = xN[1][2];
	mat[2][0] = 1, mat[2][1] = xN[3][0], mat[2][2] = xN[3][2];
	c[2] = detMat33(mat);
	//c3
	mat[0][0] = 1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][2];
	mat[1][0] = 1, mat[1][1] = xN[1][0], mat[1][2] = xN[1][2];
	mat[2][0] = 1, mat[2][1] = xN[2][0], mat[2][2] = xN[2][2];
	c[3] = -detMat33(mat);
	//=====d
	//d0
	mat[0][0] = 1, mat[0][1] = xN[1][0], mat[0][2] = xN[1][1];
	mat[1][0] = 1, mat[1][1] = xN[2][0], mat[1][2] = xN[2][1];
	mat[2][0] = 1, mat[2][1] = xN[3][0], mat[2][2] = xN[3][1];
	d[0] = -detMat33(mat);
	//d1
	mat[0][0] = 1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][1];
	mat[1][0] = 1, mat[1][1] = xN[2][0], mat[1][2] = xN[2][1];
	mat[2][0] = 1, mat[2][1] = xN[3][0], mat[2][2] = xN[3][1];
	d[1] = detMat33(mat);
	//d2
	mat[0][0] = 1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][1];
	mat[1][0] = 1, mat[1][1] = xN[1][0], mat[1][2] = xN[1][1];
	mat[2][0] = 1, mat[2][1] = xN[3][0], mat[2][2] = xN[3][1];
	d[2] = -detMat33(mat);
	//d3
	mat[0][0] = 1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][1];
	mat[1][0] = 1, mat[1][1] = xN[1][0], mat[1][2] = xN[1][1];
	mat[2][0] = 1, mat[2][1] = xN[2][0], mat[2][2] = xN[2][1];
	d[3] = detMat33(mat);
	//====total V==
	double V6 = 6.0 * (this->eleVolume(xN));

	//===================
	B->zero();
	for (int i = 0; i < 4; i++)
	{
		B->setCoeff(0, 3 * i, b[i]/V6);
		B->setCoeff(1, 3 * i + 1, c[i] / V6);
		B->setCoeff(2, 3 * i + 2, d[i] / V6);
		B->setCoeff(3, 3 * i, c[i] / V6);
		B->setCoeff(3, 3 * i + 1, b[i] / V6);
		B->setCoeff(4, 3 * i + 1, d[i] / V6);
		B->setCoeff(4, 3 * i + 2, c[i] / V6);
		B->setCoeff(5, 3 * i, d[i] / V6);
		B->setCoeff(5, 3 * i + 2, b[i] / V6);
	}
}

void pdfemEleTetra4N::eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3])
{
	// Nsigma[0][4]--sig_x[4],Nsigma[1]--sig_y....
	Matrix* B = new Matrix(6, 12);
	Vector* epsilon, * sigma;
	epsilon = new Vector(6);
	sigma = new Vector(6);
	BmatFEM(B, xN);
	matoperat.matMultiply(B, Ue, epsilon);
	matoperat.matMultiply(D, epsilon, sigma);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			//j: sx sy sz ...
			//i:node;
			Nsigma[j]->setCoeff(i, sigma->d_getCoeff(j));
		}
	}
	
	delete epsilon, delete sigma, delete B;
	epsilon = NULL; sigma = NULL, B = NULL;

}

void pdfemEleTetra4N::print_vtk(ofstream& fout, int* eleNodeID)
{
	fout << 4 << ' ';
	for (int i = 0; i < 4; i++)
	{
		fout << eleNodeID[i] - 1 << ' ';
	}
	fout << endl;
}

double pdfemEleTetra4N::eleVolume(double xN[][3])
{
	double mat[3][3];
	mat[0][0] = xN[1][0] - xN[0][0], mat[0][1] = xN[1][1] - xN[0][1], mat[0][2] = xN[1][2] - xN[0][2];
	mat[1][0] = xN[2][0] - xN[0][0], mat[1][1] = xN[2][1] - xN[0][1], mat[1][2] = xN[2][2] - xN[0][2];
	mat[2][0] = xN[3][0] - xN[0][0], mat[2][1] = xN[3][1] - xN[0][1], mat[2][2] = xN[3][2] - xN[0][2];
	return (abs(detMat33(mat))/6.0);
}
