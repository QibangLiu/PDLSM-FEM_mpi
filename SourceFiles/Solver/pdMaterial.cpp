/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdMaterial.h"

pdMaterial::pdMaterial(double e, double nu, double rho,double KIc,double sigult)
{
	cd_e = e;
	cd_nu = nu;
	cd_sigmaUlt = sigult;
	cd_rho = rho;
	cd_KIc = KIc;
}

pdMaterial::~pdMaterial()
{
}

double pdMaterial::getE() const
{
	return cd_e;
}

double pdMaterial::getnu() const
{
	return cd_nu;
}

double pdMaterial::getrho() const
{
	return cd_rho;
}

double pdMaterial::getKIc() const
{
	return cd_KIc;
}

double pdMaterial::getSigult() const
{
	return cd_sigmaUlt;
}



void pdMaterial::print(ofstream & fout)
{
	fout << "          Young's modules=" << cd_e << endl;
	fout << "             poison ratio=" << cd_nu << endl;
	fout << "                Density=" << cd_rho << endl;
}
