#include "pdfemEleTri3N.h"
extern calMatrixOperations matoperat;
pdfemEleTri3N::pdfemEleTri3N(int id, int* nId, int algoType, dataLev2* p_datLev2):
	pdfemEles(id, 3, nId, algoType, p_datLev2)
{
	ci_eleType = 5;
}

pdfemEleTri3N::~pdfemEleTri3N()
{
}

void pdfemEleTri3N::eleCenter(double xc[], double xN[][3])
{
	for (int i = 0; i < 3; i++)
	{
		xc[i] = (xN[0][i] + xN[1][i] + xN[2][i]) / 3.0;
	}
}

void pdfemEleTri3N::eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3])
{
	Matrix* B = new Matrix(3, 6);
	Matrix* Bt, * BtD;
	Bt = new Matrix(6, 3);
	BtD = new Matrix(6, 3);
	double V = eleVolume(xN);
	BmatFEM(B, xN);
	matoperat.matTranspose(B, Bt);
	matoperat.matMultiply(Bt, D, BtD);
	matoperat.matMultiply(BtD, B, Ke);
	matoperat.matMultiply(Ke, V, Ke);
	delete B, delete Bt, delete BtD;
	B = NULL, Bt = NULL, BtD = NULL;
}

void pdfemEleTri3N::eleMassMat(Matrix* Me, double rho, double xN[][3])
{
	double ve, fac;
	ve = eleVolume(xN);
	fac = rho * ve / 12.0;

	Me->zero();
	for (int i = 0; i < 6; i++)
	{
		Me->setCoeff(i, i, fac * 2);
	}
	int col;
	for (int row = 0; row < 6; row++)
	{
		col = row + 2;
		while (col < 6)
		{
			Me->setCoeff(row, col, fac);
			Me->setCoeff(col, row, fac);
			col = col + 2;
		}
	}
}

void pdfemEleTri3N::eleEquivNodalForce(Vector* Fe, double t, double xN[][3])
{
}

double pdfemEleTri3N::detJacobi(double xN[][3], double p, double q, double r)
{
	return 0.0;
}

void pdfemEleTri3N::eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3])
{
	Matrix* B = new Matrix(3, 6);
	Vector* epsilon, * sigma;
	epsilon = new Vector(3);
	sigma = new Vector(3);
	BmatFEM(B, xN);
	matoperat.matMultiply(B, Ue, epsilon);
	matoperat.matMultiply(D, epsilon, sigma);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			//j: sx sy sxy ...
			//i:node;
			Nsigma[j]->setCoeff(i, sigma->d_getCoeff(j));
		}
	}

	delete epsilon, delete sigma, delete B;
	epsilon = NULL; sigma = NULL, B = NULL;
}

void pdfemEleTri3N::print_vtk(ofstream& fout, int* eleNodeID)
{
	if (ci_eleType==5)
	{
		fout << 3 << ' ';
		for (int i = 0; i < 3; i++)
		{
			fout << eleNodeID[i] - 1 << ' ';
		}
		fout << endl;
	}
	else
	{
		cout << "ERROR: element type is wrong\n";
		exit(0);
	}
}

double pdfemEleTri3N::eleVolume(double xN[][3])
{
	double mat[3][3];
	mat[0][0] =1, mat[0][1] = xN[0][0], mat[0][2] = xN[0][1];
	mat[1][0] = 1, mat[1][1] = xN[1][0], mat[1][2] = xN[1][1];
	mat[2][0] = 1, mat[2][1] = xN[2][0], mat[2][2] = xN[2][1];
	return (0.5 * abs(detMat33(mat)));
}

void pdfemEleTri3N::BmatFEM(Matrix* B, double xN[][3])
{
	double a[3], b[3], c[3];
	/*a[0] = xN[1][0] * xN[2][1] - xN[1][1] * xN[2][0];
	a[1] = xN[2][0] * xN[0][1] - xN[2][1] * xN[0][0];
	a[2] = xN[0][0] * xN[1][1] - xN[0][1] * xN[1][0];*/
	//==
	b[0] = xN[1][1] - xN[2][1];
	b[1] = xN[2][1] - xN[0][1];
	b[2] = xN[0][1] - xN[1][1];
	//==
	c[0] = -xN[1][0] + xN[2][0];
	c[1] = -xN[2][0] + xN[0][0];
	c[2] = -xN[0][0] + xN[1][0];

	//====total V==
	double V2 = 2.0 * (this->eleVolume(xN));
	//===================
	B->zero();
	for (int i = 0; i < 3; i++)
	{
		B->setCoeff(0, 2 * i, b[i] / V2);
		B->setCoeff(1, 2 * i + 1, c[i] / V2);
		B->setCoeff(2, 2 * i, c[i] / V2);
		B->setCoeff(2, 2 * i + 1, b[i] / V2);
	}
}
