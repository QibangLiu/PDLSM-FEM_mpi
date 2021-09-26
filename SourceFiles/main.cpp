/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */



#include"datModel.h"
#include<iostream>
#include<pdfem_mpi.h>
//#include<omp.h>
using namespace std;
//==some global variables
calMatrixOperations matoperat;
pdGaussPt o_globGP(2);
pdGaussPt o_sigGP(2);// this GP is for stress extraplation;
pdGaussPt o_sifGP(4);// this GP is for SIFs calculation;
int main(int argc, char* argv[])
{
	//omp_set_num_threads(4);
	MPI_Init(&argc, &argv);// message passing interface
	pdfem_mpi po_pdfem(argc, argv);
	MPI_Finalize();
	return 0;
}