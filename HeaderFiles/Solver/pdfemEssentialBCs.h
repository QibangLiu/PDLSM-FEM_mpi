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
private:
	pdfemEssentialBCs();
	int ci_numNODE;
	int ci_dof;
	double cd_value;
	
};

