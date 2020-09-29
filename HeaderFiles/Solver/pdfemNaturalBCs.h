#pragma once
#include"pdfemEleQuad4N.h"
#include"pdfemEleLine2N.h"
#include"Matrix.h"
#include"calMatrixOperations.h"
#include"Vector.h"
using namespace std;
class pdfemNaturalBCs
{
public:
	pdfemNaturalBCs(int numEle, double val);
	~pdfemNaturalBCs();
	int getNumEle();
	double getValue();
	pdfemEles* op_getNBCsEle(int ele);
	void initialEles(int ele, int conNID[],int eleType);

	//int** cip2_nbcNID;
	//int* cip_nbcNormDire;
private:
	pdfemNaturalBCs();
	pdfemEles** cop2_NBCsEle;
	int ci_numEle;
	double cd_valule;
};

