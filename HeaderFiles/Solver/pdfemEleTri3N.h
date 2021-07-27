#pragma once
#include "pdfemEles.h"
class pdfemEleTri3N :
	public pdfemEles
{
public:
	pdfemEleTri3N(int id, int* nId, int algoType, dataLev2* p_datLev2);
	~pdfemEleTri3N();
	void eleCenter(double xc[], double xN[][3]);
	void eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3]);
	void eleMassMat(Matrix* Me, double rho, double xN[][3]);
	void eleEquivNodalForce(Vector* Fe, double t, double xN[][3]);
	double detJacobi(double xN[][3], double p, double q, double r);
	void eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3]);
	void print_vtk(ofstream& fout, int* eleNodeID);
	double eleVolume(double xN[][3]);
private:
	void BmatFEM(Matrix* B, double xN[][3]);
};

