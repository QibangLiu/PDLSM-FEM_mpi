/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdVaryEssentialBCs.h"



pdVaryEssentialBCs::pdVaryEssentialBCs(int Nid, int direc)
{
	ci_Nid = Nid;
	ci_direction = direc;
}


pdVaryEssentialBCs::~pdVaryEssentialBCs()
{
}

int pdVaryEssentialBCs::getNid() const
{
	return ci_Nid;
}

int pdVaryEssentialBCs::getDirec() const
{
	return ci_direction;
}
