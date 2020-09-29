#pragma once
#include<fstream>
using namespace std;
class pdMaterial
{
public:
	pdMaterial(double e,double nu, double rho,double KIc, double sigult);
	~pdMaterial();
	double getE()const;
	double getnu()const;
	double getrho()const;
	double getKIc()const;
	void print(ofstream &fout);
private:
	pdMaterial();
	double cd_e;
	double cd_nu;
	double cd_rho;
	double cd_KIc;
	double cd_sigmaUlt;
};

