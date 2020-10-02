// material points data
#pragma once
#include<fstream>
//#include"pdSetsolve.h"
#include"pdDof.h"
#include "dataLev2.h"
#include<iomanip>
using namespace std;
class pdNode
{
public:
	pdNode(int id, dataLev2* op_datLev2);
	~pdNode();
	int getId()const;
	void getcoor(double x[])const;
	double getvolume()const;
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
	dataLev2* cop_datLev2;
	int ci_NodeId;//Node id;
	double cd_dv;//volume ;
	//double *cd_x;//coordinate
	pdDof *cop_dof[3];
	int ci_famID;
	double cd_sigma[6];
	double cd_localDamage;
	int ci_nodeType; // 0--pure fem node; 1-- interface node; 2-- pure PD node;

};
