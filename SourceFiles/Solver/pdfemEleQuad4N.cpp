#include "pdfemEleQuad4N.h"
extern calMatrixOperations matoperat;
extern pdGaussPt o_globGP;
extern pdGaussPt o_sigGP;
pdfemEleQuad4N::pdfemEleQuad4N(int id, int* nId, int algoType, dataLev2* p_datLev2) :
	pdfemEles(id, 4, nId, algoType,p_datLev2)
{
}

pdfemEleQuad4N::~pdfemEleQuad4N()
{

}

void pdfemEleQuad4N::eleCenter(double xc[], double xN[][3])
{
	xc[0] = 0; xc[1] = 0; xc[2] = 0;
	double N[4];
	shapeFunction(N, 0, 0, 0);
	for (int n = 0; n < ci_numNodes; n++)
	{
		for (int i = 0; i < 3; i++)
		{
			xc[i] = xc[i] + xN[n][i] * N[n];
		}
	}
}

void pdfemEleQuad4N::shapeFunction(double N[], double p, double q, double r)
{
	N[0] = 0.25 * (1 - p) * (1 - q);
	N[1] = 0.25 * (1 + p) * (1 - q);
	N[2] = 0.25 * (1 + p) * (1 + q);
	N[3] = 0.25 * (1 - p) * (1 + q);
}

double pdfemEleQuad4N::detJacobi(double xN[][3], double p, double q, double r)
{
	double J[2][2] = { 0.0 };
	double dNdp[4], dNdq[4];
	dNdp[0] = 0.25 * (-1 + q);
	dNdp[1] = 0.25 * (1 - q);
	dNdp[2] = 0.25 * (1 + q);
	dNdp[3] = 0.25 * (-1 - q);
	dNdq[0] = 0.25 * (-1 + p);
	dNdq[1] = 0.25 * (-1 - p);
	dNdq[2] = 0.25 * (1 + p);
	dNdq[3] = 0.25 * (1 - p);
	for (int i = 0; i < 4; i++)
	{
		J[0][0] = J[0][0] + dNdp[i] * xN[i][0];
		J[0][1] = J[0][1] + dNdp[i] * xN[i][1];
		J[1][0] = J[1][0] + dNdq[i] * xN[i][0];
		J[1][1] = J[1][1] + dNdq[i] * xN[i][1];
	}
	double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	return abs(detJ);
}

void pdfemEleQuad4N::eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3])
{
	if (flag==1)
	{
		//extrapolation method
		int nG;
		Vector* epsilon, * sigma, * Vgp_sig[3];//for gp's stress;
		Matrix* B;
		//s_x,sy,sz,s_xy,s_yz,s_zx;
		B = new Matrix(3, 8);
		epsilon = new Vector(3);
		sigma = new Vector(3);
		for (int i = 0; i < 3; i++)
		{
			Vgp_sig[i] = new Vector(4);
		}
		nG = o_sigGP.i_getNumPts();//==always =2;
		double p, q, r = 0;
		//get stresses of each gauss point and fit rhs;
		int count = 0;
		for (int nq = 0; nq < nG; nq++)
		{
			q = o_sigGP.d_getGaussPt(nq);
			for (int np = 0; np < nG; np++)
			{
				p = o_sigGP.d_getGaussPt(np);
				BmatFEM(B, xN, p, q, r);
				matoperat.matMultiply(B, Ue, epsilon);
				matoperat.matMultiply(D, epsilon, sigma);
				for (int i = 0; i < 3; i++)
				{
					//i: sx sy sz ...
					//count: gs point;
					Vgp_sig[i]->setCoeff(count, sigma->d_getCoeff(i));
				}
				count++;
			}
		}
		//=========
		for (int i = 0; i < 3; i++)
		{
			matoperat.matMultiply(L, Vgp_sig[i], Nsigma[i]);
		}
		//===
		delete  epsilon, sigma, B;
		B = NULL; sigma = NULL;  epsilon = NULL;
		for (int i = 0; i < 3; i++)
		{
			delete Vgp_sig[i];
			Vgp_sig[i] = NULL;
		}
	}
	else if (flag==2)
	{
		int nG;
		Vector* epsilon, * sigma;
		Matrix* B;
		Vector** rhsFIT, ** coeFIT;
		B = new Matrix(3, 8);
		epsilon = new Vector(3);
		sigma = new Vector(3);
		rhsFIT = new Vector * [3];
		coeFIT = new Vector * [3];
		for (int i = 0; i < 3; i++)
		{
			rhsFIT[i] = new Vector(3);
			coeFIT[i] = new Vector(3);
		}
		nG = o_globGP.i_getNumPts();
		double p, q, tempsigma;
		// initial rhsFIT;
		for (int i = 0; i < 3; i++)
		{
			rhsFIT[i]->zero();
		}
		//get stresses of each gauss point and fit rhs;
		for (int np = 0; np < nG; np++)
		{
			p = o_globGP.d_getGaussPt(np);
			for (int nq = 0; nq < nG; nq++)
			{
				q = o_globGP.d_getGaussPt(nq);
				BmatFEM(B, xN, p, q, 0);
				matoperat.matMultiply(B, Ue, epsilon);
				matoperat.matMultiply(D, epsilon, sigma);
				//get rhs for fitting;
				for (int ii = 0; ii < 3; ii++)
				{
					rhsFIT[ii]->addCoeff(0, sigma->d_getCoeff(ii));
					rhsFIT[ii]->addCoeff(1, (sigma->d_getCoeff(ii)) * p);
					rhsFIT[ii]->addCoeff(2, (sigma->d_getCoeff(ii)) * q);
				}
			}
		}
		//cal fit coefficients and nodal stress;
		double p_I[4] = { -1.0,1.,1.,-1. };
		double q_I[8] = { -1.0,-1.0,1.0,1.0 };
		for (int ii = 0; ii < 3; ii++)
		{
			matoperat.PLUSolve(L, rhsFIT[ii], coeFIT[ii]);
			for (int ni = 0; ni < 4; ni++)
			{
				tempsigma = coeFIT[ii]->d_getCoeff(0) + coeFIT[ii]->d_getCoeff(1) * p_I[ni] +
					+coeFIT[ii]->d_getCoeff(2) * q_I[ni];
				//Nsigma[ni][ii] = tempsigma;
				Nsigma[ii]->setCoeff(ni, tempsigma);
			}
		}
		
		delete  epsilon, sigma, B;
		B = NULL; sigma = NULL;  epsilon = NULL;
		for (int i = 0; i < 3; i++)
		{
			delete rhsFIT[i];
			delete coeFIT[i];
			rhsFIT[i] = NULL;
			coeFIT[i] = NULL;
		}
		delete[] rhsFIT, coeFIT;
		rhsFIT = NULL; coeFIT = NULL;
	}
	else
	{
		printf("ERROR: flag must be 1 or 2\n");
		exit(0);
	}
}

void pdfemEleQuad4N::print_vtk(ofstream& fout, int* eleNodeID)
{
	switch (ci_eleType)
	{
	default:
		cout << "ERROR: element type is wrong\n";
		exit(0);
	case 5:
		break;// 3; // triangle
	case 9:
		fout << 4 << ' ';
		for (int i = 0; i < 4; i++)
		{
			fout << eleNodeID[i] - 1 << ' ';
		}
		fout << endl;
		break;// 4; //quad
	}

}

double pdfemEleQuad4N::eleVolume(double xN[][3])
{
	double VolEle = 0;
	double wp, wq, detJ, p, q, r = 0;
	int nG = o_globGP.i_getNumPts();
	for (int mp = 0; mp < nG; mp++)
	{
		wp = o_globGP.d_getWeight(mp);
		p = o_globGP.d_getGaussPt(mp);
		for (int mq = 0; mq < nG; mq++)
		{
			wq = o_globGP.d_getWeight(mq);
			q = o_globGP.d_getGaussPt(mq);
			detJ = detJacobi(xN, p, q, r);
			VolEle = VolEle + wp * wq * detJ;
		}
	}
	return VolEle;
}

void pdfemEleQuad4N::JacobiMat(Matrix* J, double xN[][3], double p, double q, double r)
{
	J->zero();
	double dNdp[4], dNdq[4];

	dNdp[0] = 0.25 * (-1 + q);
	dNdp[1] = 0.25 * (1 - q);
	dNdp[2] = 0.25 * (1 + q);
	dNdp[3] = 0.25 * (-1 - q);

	dNdq[0] = 0.25 * (-1 + p);
	dNdq[1] = 0.25 * (-1 - p);
	dNdq[2] = 0.25 * (1 + p);
	dNdq[3] = 0.25 * (1 - p);

	for (int i = 0; i < 4; i++)
	{
		J->addCoeff(0, 0, dNdp[i] * xN[i][0]);
		J->addCoeff(0, 1, dNdp[i] * xN[i][1]);
		J->addCoeff(1, 0, dNdq[i] * xN[i][0]);
		J->addCoeff(1, 1, dNdq[i] * xN[i][1]);
	}
}

void pdfemEleQuad4N::BmatFEM(Matrix* B, double xN[][3], double p, double q, double r)
{
	//B matrix calculated from FEM with size 3*8;
	B->zero();
	double B1[4], B2[4], invJ[2][2];
	Matrix* J = new Matrix(2, 2);
	JacobiMat(J, xN, p, q, r);

	double detJ = J->d_getCoeff(0,0) * J->d_getCoeff(1,1) -
		J->d_getCoeff(0,1) * J->d_getCoeff(1,0);

	invJ[0][0] = J->d_getCoeff(1, 1) / detJ;
	invJ[0][1] = -J->d_getCoeff(0, 1) / detJ;
	invJ[1][0] = -J->d_getCoeff(1, 0) / detJ;
	invJ[1][1] = J->d_getCoeff(0, 0) / detJ;

	double dNdp[4], dNdq[4];

	dNdp[0] = 0.25 * (-1 + q);
	dNdp[1] = 0.25 * (1 - q);
	dNdp[2] = 0.25 * (1 + q);
	dNdp[3] = 0.25 * (-1 - q);

	dNdq[0] = 0.25 * (-1 + p);
	dNdq[1] = 0.25 * (-1 - p);
	dNdq[2] = 0.25 * (1 + p);
	dNdq[3] = 0.25 * (1 - p);

	for (int i = 0; i < 4; i++)
	{
		B1[i] = invJ[0][0] * dNdp[i] + invJ[0][1] * dNdq[i];
		B2[i] = invJ[1][0] * dNdp[i] + invJ[1][1] * dNdq[i];
	}

	for (int i = 0; i < 4; i++)
	{
		B->setCoeff(0, 2 * i, B1[i]);
		B->setCoeff(1, 2 * i + 1, B2[i]);
		B->setCoeff(2, 2 * i, B2[i]);
		B->setCoeff(2, 2 * i + 1, B1[i]);
	}
	delete J; J = NULL;
}

void pdfemEleQuad4N::eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3])
{
	Ke->zero();
	Matrix* B, * Bt, * BtD, * BtDB;
	B = new Matrix(3, 8);
	Bt = new Matrix(8, 3);
	BtD = new Matrix(8, 3);
	BtDB = new Matrix(8, 8);
	//calMatrixOperations matoperat;
	int nG = o_globGP.i_getNumPts();
	double Wm, Wn, p, q, r = 0, detJ;
	for (int m = 0; m < nG; m++)
	{
		Wm = o_globGP.d_getWeight(m);
		p = o_globGP.d_getGaussPt(m);
		for (int n = 0; n < nG; n++)
		{
			Wn = o_globGP.d_getWeight(n);
			q = o_globGP.d_getGaussPt(n);
			detJ = detJacobi(xN, p, q, r);
			BmatFEM(B, xN, p, q, r);
			matoperat.matTranspose(B, Bt);
			matoperat.matMultiply(Bt, D, BtD);
			matoperat.matMultiply(BtD, B, BtDB);
			matoperat.matMultiply(BtDB, Wm * Wn * detJ, BtDB);
			matoperat.matAdd(Ke, BtDB, Ke);
		}
	}


	delete B;
	delete Bt;

	delete BtD;
	delete BtDB;
	B = NULL;
	Bt = NULL;

	BtD = NULL;
}

void pdfemEleQuad4N::eleMassMat(Matrix* Me,double rho, double xN[][3])
{
	Me->zero();
	//calMatrixOperations matoperat;
	Matrix* NtN;
	NtN = new Matrix(8, 8);
	double detJ, p, q,r=0, wp, wq;
	
	int nG = o_globGP.i_getNumPts();
	for (int np = 0; np < nG; np++)
	{
		p = o_globGP.d_getGaussPt(np);
		wp = o_globGP.d_getWeight(np);
		for (int nq = 0; nq < nG; nq++)
		{

			q = o_globGP.d_getGaussPt(nq);
			wq = o_globGP.d_getWeight(nq);
			detJ = detJacobi(xN, p, q, r);
			shapFunMat_NtN(NtN, rho, p, q, r);
			matoperat.matMultiply(NtN, detJ * wp * wq, NtN);
			matoperat.matAdd(Me, NtN, Me);
		}
	}
	delete NtN;
	NtN = NULL;
}

void pdfemEleQuad4N::eleEquivNodalForce(Vector* Fe, double t, double xN[][3])
{
	// this function is to get euivalent nodal force of element quad 4N --3D;
	Fe->zero();
	Vector* mNtN = new Vector(12);
	int nG;
	double wp, wq, p, q;
	nG = o_globGP.i_getNumPts();
	for (int mp = 0; mp < nG; mp++)
	{
		wp = o_globGP.d_getWeight(mp);
		p = o_globGP.d_getGaussPt(mp);
		for (int mq = 0; mq < nG; mq++)
		{
			wq = o_globGP.d_getWeight(mq);
			q = o_globGP.d_getGaussPt(mq);
			vec_mNtN_NBC(mNtN, xN, p, q);
			matoperat.matAdd(Fe, mNtN, Fe);
		}
	}
	matoperat.matMultiply(Fe, t, Fe);
	delete mNtN;
	mNtN = NULL;
}

void pdfemEleQuad4N::vec_mNtN_NBC(Vector* mNtN, double xN[][3], double p, double q)
{
	// this function is a auxiliary function to get equivalent nodal force of NBCs
	// this function returns mathcal{N}^T*PxPp x PxPq;
	double N[4];
	//====mNt===;
	shapeFunction(N, p, q, 0);
	Matrix* mNt = new Matrix(12, 3);
	mNt->zero();
	for (int k = 0; k < 4; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			mNt->setCoeff(i + 3 * k, i, N[k]);
		}
	}
	//===N1 N2 N3===
	double dNdp[4], dNdq[4];
	dNdp[0] = 0.25 * (-1 + q);
	dNdp[1] = 0.25 * (1 - q);
	dNdp[2] = 0.25 * (1 + q);
	dNdp[3] = 0.25 * (-1 - q);
	dNdq[0] = 0.25 * (-1 + p);
	dNdq[1] = 0.25 * (-1 - p);
	dNdq[2] = 0.25 * (1 + p);
	dNdq[3] = 0.25 * (1 - p);
	double dXdp[3] = { 0 }, dXdq[3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dXdp[i] = dXdp[i] + dNdp[j] * xN[j][i];
			dXdq[i] = dXdq[i] + dNdq[j] * xN[j][i];
		}
	}
	Vector* normN = new Vector(3);
	normN->setCoeff(0, dXdp[1] * dXdq[2] - dXdp[2] * dXdq[1]);
	normN->setCoeff(1, dXdq[0] * dXdp[2] - dXdp[0] * dXdq[2]);
	normN->setCoeff(2, dXdp[0] * dXdq[1] - dXdq[0] * dXdp[1]);
	matoperat.matMultiply(mNt, normN, mNtN);
	delete mNt, normN;
	mNt = NULL, normN = NULL;
}

void pdfemEleQuad4N::shapFunMat_NtN(Matrix* NtN, double rho, double p, double q, double r)
{
	//Nt*rho*N;
	double N[4];
	shapeFunction(N, p, q, r);
	Matrix* NN, * NNt;
	NN = new Matrix(2, 8);
	NNt = new Matrix(8, 2);
	NN->zero();
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			NN->setCoeff(i, 2 * j + i, N[j]);
		}
	}
	//calMatrixOperations matoperat;
	matoperat.matTranspose(NN, NNt);
	matoperat.matMultiply(NNt, rho, NNt);
	matoperat.matMultiply(NNt, NN, NtN);
	delete NN, NNt;
	NN = NULL; NNt = NULL;
}
