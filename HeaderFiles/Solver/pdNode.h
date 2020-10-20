// material points data
#pragma once
#include<fstream>
//#include"pdSetsolve.h"
#include"pdDof.h"
#include"dataLev2.h"
#include<iomanip>

using namespace std;
class pdNode
{
public:
	pdNode(int id, double x[], dataLev2* p_datLev2);
	~pdNode();
	int getId()const;
	void getcoor(double x[])const;
	double getvolume()const;
	void setVolume(double dv);
	void addvolume(double dv);
	pdDof *op_getDof(int i)const;
	void print(ofstream &fout);
	void printStress(ofstream &fout);
	void printStressTensor_vtk(ofstream& fout);
	void printFinalU(ofstream&PDout);
	void printDamage(ofstream &fout);
	void setLocalDamage(double phi);
	void setNodeType(int val);
	
	int getNodeType();

	void setFamID(int famID);
	int getFamID()const;
	//stress;
	void setStress(int i, double val);
	void getStress(double sigma[]);
	void addStress(int i, double val);
	void calAverageStress(int count);
private:
	pdNode();
	int ci_NodeId;//Node id;
	double cd_dv;//volume ;
	//double cd_x[3];//coordinate
	pdDof *cop_dof[3];
	int ci_famID;
	//double *cdp_sigma;
	double cd_localDamage;
	int ci_nodeType; // 0--pure fem node; 1-- interface node; 2-- pure PD node;
	dataLev2* cop_datLev2;
};
