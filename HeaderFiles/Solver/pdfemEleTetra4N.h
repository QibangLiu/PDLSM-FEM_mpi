#pragma once
#include "pdfemEles.h"
class pdfemEleTetra4N :
	public pdfemEles
{
public:
	pdfemEleTetra4N(int id, int* nId, int algoType, dataLev2* p_datLev2);
	~pdfemEleTetra4N();
	void eleCenter(double xc[], double xN[][3]);
	void eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3]);
	void eleMassMat(Matrix* Me, double rho, double xN[][3]);
	void eleEquivNodalForce(Vector* Fe, double t, double xN[][3]);
	void eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3]);
	void print_vtk(ofstream& fout, int* eleNodeID);
	double eleVolume(double xN[][3]);
private:
	void BmatFEM(Matrix* B, double xN[][3]);
	//void get_dNdX(Vector* pNpX[], double xN[][3], double p, double q, double r);
	//void shapeFunction(double N[], double p, double q, double r);
	//double detJacobi(double xN[][3], double p, double q, double r);
	//void shapFunMat_NtN(Matrix* NtN, double p, double q, double r); //for Me
	//void JacobiMat(Matrix* J, double xN[][3], double p, double q, double r);
};

