#pragma once
#include"pdfemEles.h"
using namespace std;
class pdfemEleBrick8N:
	public pdfemEles
{
public:
	pdfemEleBrick8N(int id, int *nId, int algoType);
	~pdfemEleBrick8N();
	void shapeFunction(double N[], double p, double q, double r);
	void eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3]);
	void eleMassMat(Matrix* Me, double rho, double xN[][3]);
	void eleEquivNodalForce(Vector* Fe, double t, double xN[][3]);
	double detJacobi(double xN[][3], double p, double q, double r);
	void eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3]);
	void print_vtk(ofstream& fout, int* eleNodeID);
	double eleVolume(double xN[][3]);
private:
	void shapFunMat_NtN(Matrix* NtN,  double p, double q, double r);
	void JacobiMat(Matrix* J, double xN[][3], double p, double q, double r);
	void BmatFEM(Matrix* B, double xN[][3], double p, double q, double r);
	void get_dNdX(Vector* pNpX[], double xN[][3], double p, double q, double r);

	

};

