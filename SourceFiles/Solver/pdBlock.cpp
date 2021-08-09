/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdBlock.h"


pdBlock::pdBlock(int id)
{
	ci_ID = id;
}

pdBlock::~pdBlock()
{
}

void pdBlock::putNodeInBlock(int id)
{
	cv_Nodeid.push_back(id);
}

void pdBlock::putEleInBlock(int eleID)
{
	cv_eleID.push_back(eleID);
}


int pdBlock::getNodeoB(int i)
{
	return cv_Nodeid[i];
}



int pdBlock::getNumNodeoB() const
{
	return cv_Nodeid.size();
}

int pdBlock::getEleoB(int i)
{
	return cv_eleID[i];
}

int pdBlock::getNumEleoB() const
{
	return cv_eleID.size();
}



void pdBlock::print(ofstream & fout)
{
	fout << "Block:\t" << ci_ID << endl;
	fout<<"numNode:\t"<< cv_eleID.size() << " *** ";
	/*for (int i = 0; i < ci_numNode; i++)
	{
		fout << setw(6) << cv_Nodeid[i];
	}*/
	for (int i = 0; i < cv_eleID.size(); i++)
	{
		fout << ' ' << cv_eleID[i];
	}
	fout << endl;

}
