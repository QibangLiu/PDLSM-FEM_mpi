/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */


#include "pdPointBC.h"

pdPointBC::pdPointBC(int NID, int fdof, double value)
{
	ci_NodeID = NID;
	ci_fDof = fdof;
	cd_value = value;
}

int pdPointBC::i_getfDOF() const
{
	return ci_fDof;
}


double pdPointBC::d_getValue() const
{
	return cd_value;
}

int pdPointBC::i_getNID() const
{
	return ci_NodeID;
}

void pdPointBC::setForceV(double f)
{
	cd_value = f;
}
