#include "pdfemEssentialBCs.h"

pdfemEssentialBCs::pdfemEssentialBCs(int numNode, string dof, double val)
{
	ci_numNODE = numNode;
	cd_value = val;
	cip_NID = new int[ci_numNODE];
	if (dof=="UX")
	{
		ci_dof = 0;
	}
	else if (dof=="UY")
	{
		ci_dof = 1;
	}
	else if (dof=="UZ")
	{
		ci_dof = 2;
	}
	else
	{
		printf("DOFs must be UX, UY, or UZ\n");
		exit(0);
	}
}

pdfemEssentialBCs::~pdfemEssentialBCs()
{
	delete[] cip_NID;
	cip_NID = NULL;
}

int pdfemEssentialBCs::getNumNODE()
{
	return ci_numNODE;
}

double pdfemEssentialBCs::getValue()
{
	return cd_value;
}

int pdfemEssentialBCs::getDOF()
{
	return ci_dof;
}

