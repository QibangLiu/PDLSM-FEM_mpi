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
	if (ci_rank==0)
	{
		if (!fin.is_open())
		{
			printf("ERROR: Command file \"%s\" is not exist.\n", ifName);
			exit(0);
		}
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
		printf("Elapsed time of PDLSM-FEM mpi solving of %s is %f\n", proName, t2 - t1);
		cout << "PDLSM-FEM program finished." << endl;
	}
}
