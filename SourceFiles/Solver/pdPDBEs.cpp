/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */


#include "pdPDBEs.h"

pdPDBEs::pdPDBEs(int nodeID[], int numNode)
{
	ci_numNode = numNode;
	ci_nodeID = new int[ci_numNode];
	for (int i = 0; i < ci_numNode; i++)
	{
		ci_nodeID[i] = nodeID[i];
	}
}

pdPDBEs::~pdPDBEs()
{
	delete[] ci_nodeID;
	ci_nodeID = nullptr;
}

void pdPDBEs::getNodeID(int nID[])
{
	for (int i = 0; i < ci_numNode; i++)
	{
		nID[i] = ci_nodeID[i];
	}

}

int pdPDBEs::i_getNumnode()
{
	return ci_numNode;
}
