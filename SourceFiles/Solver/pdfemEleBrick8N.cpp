#include "pdfemEleBrick8N.h"
#include<string>
extern calMatrixOperations matoperat;
extern pdGaussPt o_globGP;
extern pdGaussPt o_sigGP;
pdfemEleBrick8N::pdfemEleBrick8N(int id, int *nId, int algoType, dataLev2* p_datLev2) :
	pdfemEles(id, 8, nId, algoType,p_datLev2)
{
	/*if (nId[2]== nId[3]&& nId[4] == nId[5]
		&& nId[4] == nId[6] && nId[4] == nId[7])
	{
		ci_eleType = 10;
	}
	else if (nId[4] == nId[5] && nId[6] == nId[7] )
	{
		ci_eleType = 13;
	}
	else if (nId[4] == nId[5]&& nId[4] == nId[6] && nId[4] == nId[7])
	{
		ci_eleType = 14;
	}
	else*/
	{
		ci_eleType = 12;
	}
}

pdfemEleBrick8N::~pdfemEleBrick8N()
{
}

void pdfemEleBrick8N::eleCenter(double xc[], double xN[][3])
{
	xc[0] = 0; xc[1] = 0; xc[2] = 0;
	double N[8];
	shapeFunction(N, 0, 0, 0);
	for (int n = 0; n < ci_numNodes; n++)
	{
		for (int i = 0; i < 3; i++)
		{
			xc[i] = xc[i] + xN[n][i] * N[n];
		}
	}
}

void pdfemEleBrick8N::shapeFunction(double N[], double p, double q, double r)
{
	double p_I[8] = { -1.0,1.,1.,-1.,-1.,1.,1.,-1. };
	double q_I[8] = { -1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 };
	double r_I[8] = { -1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0 };
	for (int i = 0; i < 8; i++)
	{
		N[i] = 1.0 / 8 * (1 + p * p_I[i]) * (1 + q * q_I[i]) * (1 + r * r_I[i]);
	}
}

void pdfemEleBrick8N::shapFunMat_NtN(Matrix* NtN, double p, double q, double r)
{
	//Nt*rho*N;
	double N[8];
	shapeFunction(N, p, q, r);
	Matrix* NN, * NNt;
	NN = new Matrix(3, 24);
	NNt = new Matrix(24, 3);
	NN->zero();
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			NN->setCoeff(i, 3 * j + i, N[j]);
		}
	}
	
	matoperat.matTranspose(NN, NNt);
	//matoperat.matMultiply(NNt, rho, NNt);
	matoperat.matMultiply(NNt, NN, NtN);
	delete NN, NNt;
	NN = NULL; NNt = NULL;
}

double pdfemEleBrick8N::detJacobi(double xN[][3], double p, double q, double r)
{
	// 3D brick 8 nodes elements;
	Matrix* J = new Matrix(3, 3);
	JacobiMat(J, xN, p, q, r);
	//ofstream fout("test" + std::to_string(1) + ".out");
	//J->print(fout);
	//fout.close();
	//==det(J)===
	double det_J = 0.0;
	int rr, cc;
	double ans1, ans2;
	for (int i = 0; i < 3; i++) {
		rr = 0, cc = i;
		ans1 = 1., ans2 = 1.;
		for (int j = 0; j < 3; j++) {
			ans1 *= J->d_getCoeff(rr + j, (cc + j) % 3);
			ans2 *= J->d_getCoeff(rr + j, (cc - j + 3) % 3);
		}
		det_J = det_J + ans1 - ans2;
	}
	return (det_J); 
}

void pdfemEleBrick8N::eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3])
{
	// Nsigma[0][8]--sig_x[8],Nsigma[1]--sig_y....
	if (flag==1)
	{
		//===================extrapolation method
		Vector* epsilon, * sigma, * Vgp_sig[6];
		Matrix* B;
		//s_x,sy,sz,s_xy,s_yz,s_zx;
		B = new Matrix(6, 24);
		epsilon = new Vector(6);
		sigma = new Vector(6);
		for (int i = 0; i < 6; i++)
		{
			Vgp_sig[i] = new Vector(8);
		}
		//get stresses of each gauss point and fit rhs;
		double a = 1.0 / sqrt(3);
		double p[8] = { -a,a,a,-a,-a,a,a,-a };
		double q[8] = { -a,-a,a,a,-a,-a,a,a };
		double r[8] = { -a,-a, -a,-a, a,a,a,a };
		for (int i = 0; i < 8; i++)
		{
			BmatFEM(B, xN, p[i], q[i], r[i]);
			matoperat.matMultiply(B, Ue, epsilon);
			matoperat.matMultiply(D, epsilon, sigma);
			for (int j = 0; j < 6; j++)
			{
				//j: sx sy sz ...
				//i: gs point;
				Vgp_sig[j]->setCoeff(i, sigma->d_getCoeff(j));
			}
		}
		//=========
		for (int i = 0; i < 6; i++)
		{
			matoperat.matMultiply(L, Vgp_sig[i], Nsigma[i]);
		}

		//===
		delete  epsilon, sigma, B;
		B = NULL; sigma = NULL;  epsilon = NULL;
		for (int i = 0; i < 6; i++)
		{
			delete Vgp_sig[i];
			Vgp_sig[i] = NULL;
		}
	}
	else if (flag==2)
	{
		//LSM method
		int nG;
		Vector* epsilon, * sigma;
		Matrix* B;
		Vector** rhsFIT, ** coeFIT;
		//s_x,sy,sz,s_xy,s_yz,s_zx;
		B = new Matrix(6, 24);
		epsilon = new Vector(6);
		sigma = new Vector(6);
		rhsFIT = new Vector * [6];
		coeFIT = new Vector * [6];
		for (int i = 0; i < 6; i++)
		{
			//a b c d==
			rhsFIT[i] = new Vector(4);
			coeFIT[i] = new Vector(4);
		}
		nG = o_globGP.i_getNumPts();
		double p, q, r, tempsigma;
		// initial rhsFIT;
		for (int i = 0; i < 6; i++)
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
				for (int nr = 0; nr < nG; nr++)
				{
					r = o_globGP.d_getGaussPt(nr);
					BmatFEM(B, xN, p, q, r);
					matoperat.matMultiply(B, Ue, epsilon);
					matoperat.matMultiply(D, epsilon, sigma);
					//get rhs for fitting;
					for (int ii = 0; ii < 6; ii++)
					{
						rhsFIT[ii]->addCoeff(0, sigma->d_getCoeff(ii));
						rhsFIT[ii]->addCoeff(1, (sigma->d_getCoeff(ii)) * p);
						rhsFIT[ii]->addCoeff(2, (sigma->d_getCoeff(ii)) * q);
						rhsFIT[ii]->addCoeff(3, (sigma->d_getCoeff(ii)) * r);
					}
				}


			}
		}
		//cal fit coefficients and nodal stress;
		double p_I[8] = { -1.0,1.,1.,-1.,-1.,1.,1.,-1. };
		double q_I[8] = { -1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 };
		double r_I[8] = { -1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0 };
		for (int ii = 0; ii < 6; ii++)
		{
			matoperat.PLUSolve(L, rhsFIT[ii], coeFIT[ii]);
			for (int ni = 0; ni < 8; ni++)
			{
				tempsigma = coeFIT[ii]->d_getCoeff(0) + coeFIT[ii]->d_getCoeff(1) * p_I[ni] +
					+coeFIT[ii]->d_getCoeff(2) * q_I[ni] + coeFIT[ii]->d_getCoeff(3) * r_I[ni];
				//Nsigma[ni][ii] = tempsigma;
				Nsigma[ii]->setCoeff(ni, tempsigma);
			}
		}
		delete  epsilon, sigma, B;
		B = NULL; sigma = NULL;  epsilon = NULL;
		for (int i = 0; i < 6; i++)
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

void pdfemEleBrick8N::print_vtk(ofstream& fout, int* eleNodeID)
{
	switch (ci_eleType)
	{
	default:
		cout << "ERROR: element type is wrong\n";
		exit(0);
	case 10:
		// 0 1 2 =3; 4=5=6=7;
		fout << 4 << ' ';
		for (int i = 0; i < 3; i++)
		{
			fout << eleNodeID[i] - 1 << ' ';
		}
		fout << eleNodeID[4] - 1 << endl;
		break;// 4; //4n_tetra
	case 12:
		fout << 8 << ' ';
		for (int i = 0; i < 8; i++)
		{
			fout << eleNodeID[i] - 1 << ' ';
		}
		fout << endl;
		break;// 8;//hexahedron
	case 13:
		// 0 1 2 3; 4=5 6=7;
		fout << 6 << ' ';
		for (int i = 0; i < 4; i++)
		{
			fout << eleNodeID[i] - 1 << ' ';
		}
		fout << eleNodeID[4] - 1<<' ' << eleNodeID[6] - 1<< endl;
		break;// 6; //wedge;
	case 14:
		// 0 1 2 3; 4=5=6=7;
		fout << 5 << ' ';
		for (int i = 0; i < 4; i++)
		{
			fout << eleNodeID[i] - 1 << ' ';
		}
		fout << eleNodeID[4] - 1 << endl;
		break;// 5;//pyramid;
	}
}

double pdfemEleBrick8N::eleVolume(double xN[][3])
{
	double VolEle = 0;
	double wp, wq, wr, detJ, p, q, r;
	int nG = o_globGP.i_getNumPts();
	for (int mp = 0; mp < nG; mp++)
	{
		wp = o_globGP.d_getWeight(mp);
		p = o_globGP.d_getGaussPt(mp);
		for (int mq = 0; mq < nG; mq++)
		{
			wq = o_globGP.d_getWeight(mq);
			q = o_globGP.d_getGaussPt(mq);
			for (int mr = 0; mr < nG; mr++)
			{
				wr = o_globGP.d_getWeight(mr);
				r = o_globGP.d_getGaussPt(mr);
				detJ = detJacobi(xN, p, q, r);
				VolEle = VolEle + wp * wq * wr * detJ;
			}
		}
	}
	return VolEle;
}

void pdfemEleBrick8N::JacobiMat(Matrix* J, double xN[][3], double p, double q, double r)
{
	// 3D brick 8 nodes elements;
	double dNdp[8], dNdq[8], dNdr[8];
	double p_I[8] = { -1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0 };
	double q_I[8] = { -1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 };
	double r_I[8] = { -1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0 };
	for (int i = 0; i < 8; i++)
	{
		dNdp[i] = 1.0 / 8 * p_I[i] * (1 + q * q_I[i]) * (1 + r * r_I[i]);
		dNdq[i] = 1.0 / 8 * q_I[i] * (1 + p * p_I[i]) * (1 + r * r_I[i]);
		dNdr[i] = 1.0 / 8 * r_I[i] * (1 + p * p_I[i]) * (1 + q * q_I[i]);
	}
	// jacobi matrix===
	J->zero();
	for (int cc = 0; cc < 3; cc++)
	{
		for (int i = 0; i < 8; i++)
		{
			J->addCoeff(0, cc, dNdp[i] * xN[i][cc]);
			J->addCoeff(1, cc, dNdq[i] * xN[i][cc]);
			J->addCoeff(2, cc, dNdr[i] * xN[i][cc]);
		}
	}
}

void pdfemEleBrick8N::BmatFEM(Matrix* B, double xN[][3], double p, double q, double r)
{
	//printf("======================================================\n");
	//B mat size 6 *24;
	B->zero();
	Matrix* J = new Matrix(3, 3);
	Vector* dNdX[8], *dNdpqr[8];
	JacobiMat(J, xN, p, q, r);
	for (int i = 0; i < 8; i++)
	{
		dNdX[i] = new Vector(3);
		dNdpqr[i] = new Vector(3);
	}
	//for (int i = 0; i < 3; i++)
	//{
	//	for (int j = 0; j < 3; j++)
	//	{
	//		cout << J->d_getCoeff(i, j) << ' ';
	//	}
	//	cout << endl;
	//}
	//===dNdpqr====
	double p_I[8] = { -1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0 };
	double q_I[8] = { -1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 };
	double r_I[8] = { -1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0 };
	for (int i = 0; i < 8; i++)
	{
		dNdpqr[i]->setCoeff(0, 1.0 / 8 * p_I[i] * (1 + q * q_I[i]) * (1 + r * r_I[i]));
		dNdpqr[i]->setCoeff(1, 1.0 / 8 * q_I[i] * (1 + p * p_I[i]) * (1 + r * r_I[i]));
		dNdpqr[i]->setCoeff(2, 1.0 / 8 * r_I[i] * (1 + p * p_I[i]) * (1 + q * q_I[i]));
		matoperat.PLUSolve(J, dNdpqr[i], dNdX[i]);
	}
	/*for (int i = 0; i < 8; i++)
	{
		printf("dndx dndy dxnz = %e %e %e\n", dNdpqr[i]->d_getCoeff(0),
			dNdpqr[i]->d_getCoeff(1), dNdpqr[i]->d_getCoeff(2));
	}*/
	for (int i = 0; i < 8; i++)
	{
		B->setCoeff(0, 3 * i, dNdX[i]->d_getCoeff(0));
		B->setCoeff(1, 1 + 3 * i, dNdX[i]->d_getCoeff(1));
		B->setCoeff(2, 2 + 3 * i, dNdX[i]->d_getCoeff(2));
		B->setCoeff(3, 3 * i, dNdX[i]->d_getCoeff(1));
		B->setCoeff(3, 1 + 3 * i, dNdX[i]->d_getCoeff(0));
		B->setCoeff(4, 1 + 3 * i, dNdX[i]->d_getCoeff(2));
		B->setCoeff(4, 2 + 3 * i, dNdX[i]->d_getCoeff(1));
		B->setCoeff(5, 3 * i, dNdX[i]->d_getCoeff(2));
		B->setCoeff(5, 2 + 3 * i, dNdX[i]->d_getCoeff(0));
	}
	delete J; J = NULL;
	for (int i = 0; i < 8; i++)
	{
		delete dNdX[i], dNdpqr[i];
		dNdX[i]=NULL, dNdpqr[i]=NULL;
	}
}

void pdfemEleBrick8N::get_dNdX(Vector* dNdX[], double xN[][3], double p, double q, double r)
{
	//==dNdX [8][3]
	Matrix* J = new Matrix(3, 3);
	Vector * dNdpqr[8];
	JacobiMat(J, xN, p, q, r);
	for (int i = 0; i < 8; i++)
	{
		dNdpqr[i] = new Vector(3);
	}
	//===dNdpqr====
	double p_I[8] = { -1.0,1.,1.,-1.,-1.,1.,1.,-1. };
	double q_I[8] = { -1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 };
	double r_I[8] = { -1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0 };
	for (int i = 0; i < 8; i++)
	{
		dNdpqr[i]->setCoeff(0, 1.0 / 8 * p_I[i] * (1 + q * q_I[i]) * (1 + r * r_I[i]));
		dNdpqr[i]->setCoeff(1, 1.0 / 8 * q_I[i] * (1 + p * p_I[i]) * (1 + r * r_I[i]));
		dNdpqr[i]->setCoeff(2, 1.0 / 8 * r_I[i] * (1 + p * p_I[i]) * (1 + q * q_I[i]));
		matoperat.PLUSolve(J, dNdpqr[i], dNdX[i]);
	}
	delete J; J = NULL;
	for (int i = 0; i < 8; i++)
	{
		delete dNdpqr[i];
		dNdpqr[i] = NULL;
	}
}

void pdfemEleBrick8N::eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3])
{
	////==========optimized version====some where wrong, need debug;;
	////==initialize;
	//Ke->zero();
	//Vector* dNdX[8];
	//for (int i = 0; i < 8; i++)
	//{
	//	dNdX[i] = new Vector(3);
	//}
	//double beta = D->d_getCoeff(0, 0);
	//double mu = D->d_getCoeff(3, 3);
	//double lamda = beta - 2 * mu;
	//int cyclicAxes[5] = { 0,1,2,0,1 };
	////Gauss integration
	//int nG = o_globGP.i_getNumPts();
	//double Wp, Wq, Wr, p, q, r, detJ;
	//double temp, temp1, temp2, temp3;
	//int index_r, index_c;
	//for (int mp = 0; mp < nG; mp++)
	//{
	//	Wp = o_globGP.d_getWeight(mp);
	//	p = o_globGP.d_getGaussPt(mp);
	//	for (int mq = 0; mq < nG; mq++)
	//	{
	//		Wq = o_globGP.d_getWeight(mq);
	//		q = o_globGP.d_getGaussPt(mq);
	//		for (int mr = 0; mr < nG; mr++)
	//		{
	//			Wr = o_globGP.d_getWeight(mr);
	//			r = o_globGP.d_getGaussPt(mr);
	//			detJ = detJacobi(xN, p, q, r);
	//			get_dNdX(dNdX, xN, p, q, r);
	//			for (int m = 0; m < 8; m++)
	//			{
	//				//node m;
	//				for (int n = 0; n < 8; n++)
	//				{
	//					//node n;
	//					for (int i = 0; i < 3; i++)
	//					{
	//						//axis i
	//						for (int j = 0; j < 3; j++)
	//						{
	//							// axis j;
	//							index_r = 3 * m + i;
	//							index_c = 3 * n + j;
	//							if (i==j)
	//							{
	//								temp1 = dNdX[m]->d_getCoeff(i) * dNdX[n]->d_getCoeff(i);
	//								temp2 = dNdX[m]->d_getCoeff(cyclicAxes[i + 1]);
	//								temp3 = dNdX[m]->d_getCoeff(cyclicAxes[i + 2]);
	//								temp = beta * temp1 + mu * (temp2 * temp2 + temp3 * temp3);
	//								temp = temp * Wp * Wq * Wr * detJ;
	//								Ke->addCoeff(index_r, index_c, temp);
	//							}
	//							else
	//							{
	//								temp1 =lamda* dNdX[m]->d_getCoeff(i) * dNdX[n]->d_getCoeff(j);
	//								temp2 = mu * dNdX[m]->d_getCoeff(j) * dNdX[n]->d_getCoeff(i);
	//								temp = (temp1 + temp2) * Wp * Wq * Wr * detJ;
	//								Ke->addCoeff(index_r, index_c, temp);
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//for (int i = 0; i < 8; i++)
	//{
	//	delete dNdX[i];
	//	dNdX[i] = NULL;
	//}


	//=====================matrix multiplication method;
	Ke->zero();
	Matrix* B, * Bt, * BtD, * BtDB;
	B = new Matrix(6, 24);
	Bt = new Matrix(24, 6);
	BtD = new Matrix(24, 6);
	BtDB = new Matrix(24, 24);
	int nG = o_globGP.i_getNumPts();
	double Wp, Wq,Wr, p, q, r, detJ;
	for (int mr = 0; mr < nG; mr++)
	{
		Wr = o_globGP.d_getWeight(mr);
		r = o_globGP.d_getGaussPt(mr);
		for (int mq = 0; mq < nG; mq++)
		{
			Wq = o_globGP.d_getWeight(mq);
			q = o_globGP.d_getGaussPt(mq);
			for (int mp = 0; mp < nG; mp++)
			{
				Wp = o_globGP.d_getWeight(mp);
				p = o_globGP.d_getGaussPt(mp);
				detJ = detJacobi(xN, p, q, r);
				
				BmatFEM(B, xN, p, q, r);
				matoperat.matTranspose(B, Bt);
				matoperat.matMultiply(Bt, D, BtD);
				matoperat.matMultiply(BtD, B, BtDB);
				matoperat.matMultiply(BtDB, Wp * Wq*Wr * detJ, BtDB);
				matoperat.matAdd(Ke, BtDB, Ke);
			}
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

void pdfemEleBrick8N::eleMassMat(Matrix* Me, double rho, double xN[][3])
{
	Me->zero();
	Matrix* NtN;
	NtN = new Matrix(24, 24);
	double detJ, p, q, r, wp, wq, wr;
	int nG = o_globGP.i_getNumPts();
	for (int np = 0; np < nG; np++)
	{
		p = o_globGP.d_getGaussPt(np);
		wp = o_globGP.d_getWeight(np);
		for (int nq = 0; nq < nG; nq++)
		{
			q = o_globGP.d_getGaussPt(nq);
			wq = o_globGP.d_getWeight(nq);
			for (int nr = 0; nr < nG; nr++)
			{
				r = o_globGP.d_getGaussPt(nr);
				wr = o_globGP.d_getWeight(nr);
				detJ = detJacobi(xN, p, q, r);
				shapFunMat_NtN(NtN, p, q, r);
				matoperat.matMultiply(NtN, rho * detJ * wp * wq * wr, NtN);
				matoperat.matAdd(Me, NtN, Me);
			}
		}
	}
	delete NtN;
	NtN = NULL;
}

void pdfemEleBrick8N::eleEquivNodalForce(Vector* Fe, double t, double xN[][3])
{
}

