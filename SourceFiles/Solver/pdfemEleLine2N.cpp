#include "pdfemEleLine2N.h"

extern calMatrixOperations matoperat;
extern pdGaussPt o_globGP;


pdfemEleLine2N::pdfemEleLine2N(int id, int* nId, int algoType, dataLev2* p_datLev2) :
	pdfemEles(id, 2, nId, algoType,p_datLev2)
{

}

pdfemEleLine2N::~pdfemEleLine2N()
{
}

void pdfemEleLine2N::eleCenter(double xc[], double xN[][3])
{
}

void pdfemEleLine2N::shapeFunction(double N[], double p, double q, double r)
{
}

void pdfemEleLine2N::eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3])
{
}

void pdfemEleLine2N::eleMassMat(Matrix* Me, double rho, double xN[][3])
{
}

void pdfemEleLine2N::eleEquivNodalForce(Vector* Fe, double t, double xN[][3])
{
	Fe->zero();
	//gauss integration
	int nG;
	double p, wp, N[2], dNdp[2], temp1, temp2, temp;
	nG = o_globGP.i_getNumPts();
	for (int m = 0; m < nG; m++)
	{
		p = o_globGP.d_getGaussPt(m);
		wp = o_globGP.d_getWeight(m);
		N[0] = 0.5 * (1 - p);
		N[1] = 0.5 * (1 + p);
		dNdp[0] = -0.5;
		dNdp[1] = 0.5;
		for (int ii = 0; ii < 2; ii++)
		{
			//ii is node index
			//get dNdp[i]y[i] and dNdp[i]x[i]
			temp1 = 0; temp2 = 0;
			for (int jj = 0; jj < 2; jj++)
			{
				temp1 = temp1 + dNdp[jj] * xN[jj][1];
				temp2 = temp2 + dNdp[jj] * xN[jj][0];
			}
			/*fN[ii][0] = fN[ii][0] + Wm * N[ii] * (tn * temp1 + ts * temp2);
			fN[ii][1] = fN[ii][1] + Wm * N[ii] * (-tn * temp2 + ts * temp1);*/
			temp = wp * N[ii] * (t * temp1);
			Fe->addCoeff(ii * 2, temp);
			temp = wp * N[ii] * (-t * temp2);
			Fe->addCoeff(ii * 2 + 1, temp);

		}
	}
}

double pdfemEleLine2N::detJacobi(double xN[][3], double p, double q, double r)
{
	return 0.0;
}

void pdfemEleLine2N::eleFitStresses(int flag,Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3])
{
}


void pdfemEleLine2N::print_vtk(ofstream& fout, int* eleNodeID)
{
	fout << 2<<' ' << eleNodeID[0] - 1 << ' ' << eleNodeID[1] - 1 << endl;
}

double pdfemEleLine2N::eleVolume(double xN[][3])
{
	return 0.0;
}


void pdfemEleLine2N::vec_mNtN_NBC(Vector* mNtN, double xN[][3], double p, double q)
{
}

