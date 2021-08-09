/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdfem_mpi.h"
#include"datModel.h"
#include"fioFile.h"
#include "pdsolve.h"
#include<mpi.h>
pdfem_mpi::pdfem_mpi(int argc, char* argv[])
{

	MPI_Comm_rank(MPI_COMM_WORLD, &ci_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ci_numProce);

	//reading data to data model==========================
	fioFiles o_files(ci_rank);
	datModel o_modeldata;
	//===========input file=======
	ifstream fin;		//fin for input file
	string ifName; 
	if (argc>1)
	{
		ifName = argv[1];
	}
	else
	{
		ifName = "PDFEMcmd.in";
	}
	fin.open(ifName);
	if (!fin.is_open())
	{
		if (ci_rank==0)
		{
			printf("ERROR: Command file \"%s\" is not exist.\n", ifName.c_str());
		}
		exit(0);
	}
	o_files.CMDfile(o_modeldata, fin);
	fin.close();
	//o_modeldata.writeData();
	//===========solving===========================
	double t1, t2;
	if (ci_rank==0)
	{
		const char* proName = o_modeldata.cs_title.c_str();
		printf("PDLSM-FEM solving project of %s.........\n", proName);
		
	}
	pdsolve o_sol(o_modeldata, ci_rank, ci_numProce);
	if (ci_rank==0)
	{
		t1 = MPI_Wtime();
	}
	o_sol.pdfemSolver(o_modeldata, o_files,argv);
	if (ci_rank==0)
	{
		double t2 = MPI_Wtime();
		const char *proName = o_modeldata.cs_title.c_str();
		printf("Elapsed time of PDLSM-FEM mpi solving of %s is %f s\n", proName, t2 - t1);
		cout << "PDLSM-FEM program finished." << endl;
	}
}
