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
	cd_value = val;
	cop2_NBCsEle = new pdfemEles * [numEle];
	cb_varing = false;
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
	delete[] cop2_NBCsEle;
	cop2_NBCsEle = nullptr;
}

int pdfemNaturalBCs::getNumEle()
{
	return ci_numEle;
}

double pdfemNaturalBCs::getValue()
{
	return cd_value;
}

pdfemEles* pdfemNaturalBCs::op_getNBCsEle(int ele)
{
	return cop2_NBCsEle[ele];
}

void pdfemNaturalBCs::initialEles(int ele, int conNID[],int eleType,dataLev2 *p_datLev2)
{
	switch (eleType)
	{
	case 2:
	{
		cop2_NBCsEle[ele] = new pdfemEleLine2N(ele + 1, conNID, -1, p_datLev2);
		break;
	}
	case 3:
	{
		break;
	}
	case 4:
	{
		cop2_NBCsEle[ele] = new pdfemEleQuad4N(ele + 1, conNID, -1, p_datLev2);
		break;
	}
	default:
		break;
	}

}

