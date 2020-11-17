#pragma once
#include<string>
using namespace std;
class pdfemEssentialBCs
{
public:
	pdfemEssentialBCs(int numNode, string dof, double val);
	~pdfemEssentialBCs();
	int* cip_NID;
	int getNumNODE();
	double getValue();
	int getDOF();
	bool cb_varing;// this ebc is varing or not, during dynamic solving
	double cd_velocity;// velocity during dynamic solving;
	double cd_value;
private:
	pdfemEssentialBCs();
	int ci_numNODE;
	int ci_dof;
	
	
};

