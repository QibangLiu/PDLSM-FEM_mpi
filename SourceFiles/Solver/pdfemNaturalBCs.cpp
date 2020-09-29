#include "pdfemNaturalBCs.h"

pdfemNaturalBCs::pdfemNaturalBCs(int numEle, double val)
{
	//========eleType:
	//2--1D-2node;
	//3--1D -3nodes;
	//4--2D--4Nodes;
	//8--2D--8Nodes;
	//==================
	ci_numEle = numEle;
	cd_valule = val;
	cop2_NBCsEle = new pdfemEles * [numEle];
	//cip2_nbcNID = new int* [numEle];
	//ci_numNodeOfEle = numNode;
	//cip_nbcNormDire = new int[numEle];

}

pdfemNaturalBCs::~pdfemNaturalBCs()
{
	for (int i = 0; i < ci_numEle; i++)
	{
		delete cop2_NBCsEle[i];
		cop2_NBCsEle[i] = nullptr;
	}
	delete cop2_NBCsEle;
	cop2_NBCsEle = nullptr;
}

int pdfemNaturalBCs::getNumEle()
{
	return ci_numEle;
}

double pdfemNaturalBCs::getValue()
{
	return cd_valule;
}

pdfemEles* pdfemNaturalBCs::op_getNBCsEle(int ele)
{
	return cop2_NBCsEle[ele];
}

void pdfemNaturalBCs::initialEles(int ele, int conNID[],int eleType)
{
	switch (eleType)
	{
	case 2:
	{
		cop2_NBCsEle[ele] = new pdfemEleLine2N(ele + 1, conNID, -1);
		break;
	}
	case 3:
	{
		break;
	}
	case 4:
	{
		cop2_NBCsEle[ele] = new pdfemEleQuad4N(ele + 1, conNID, -1);
		break;
	}
	default:
		break;
	}

}

