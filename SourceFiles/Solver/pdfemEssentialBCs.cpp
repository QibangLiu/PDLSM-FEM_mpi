/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdfemEssentialBCs.h"

pdfemEssentialBCs::pdfemEssentialBCs(int numNode, string dof, double val)
{
	cb_varing = false;
	cd_velocity = 0;
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

