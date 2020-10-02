#pragma once
#include"pdfemEles.h"
class pdfemEleLine2N :
	public pdfemEles
{
public:
	pdfemEleLine2N(int id, int* nId, int algoType);
	~pdfemEleLine2N();
	void shapeFunction(double N[], double p, double q, double r);
	void eleStiffMatFEM(Matrix* Ke, Matrix* D, double xN[][3]);
	void eleMassMat(Matrix* Me, double rho, double xN[][3]);
	void eleEquivNodalForce(Vector* Fe, double t, double xN[][3]);
	double detJacobi(double xN[][3], double p, double q, double r);
	void eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3]);
	void print_vtk(ofstream& fout, int* eleNodeID);
	double eleVolume(double xN[][3]);
private:
	void vec_mNtN_NBC(Vector* mNtN, double xN[][3], double p, double q);
};

