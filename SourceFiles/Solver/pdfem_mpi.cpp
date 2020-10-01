#include "pdfem_mpi.h"
#include"datModel.h"
#include"fioFile.h"
#include "pdsolve.h"
#include<mpi.h>
pdfem_mpi::pdfem_mpi(int argc, char* argv[])
{

	MPI_Comm_rank(MPI_COMM_WORLD, &ci_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ci_numProce);
	//=====OpenMP thread seting
	int NUM_THREAD;
	if (argc == 3)
	{
		 NUM_THREAD = atoi(argv[2]);
		omp_set_num_threads(NUM_THREAD);
	}
	else
	{
		NUM_THREAD = omp_get_num_procs();
	}
	if (ci_rank==0)
	{

		printf("Number of OpenMP thread: %d\n", NUM_THREAD);
	}
	
	//===========input file=======
	ifstream fin;		//fin for input file
	char* ifName = argv[1];
	//string ifName = "Cubic_ebcFEM.lis";
	fin.open(ifName);
	if (ci_rank==0)
	{
		if (!fin.is_open())
		{
			printf("File %s is not exist\n", ifName);
			exit(0);
		}
	}
	

	//====write results flag;
	int wflag = 3;
	//reading data to data model==========================
	datModel o_modeldata;
	printf("Reading data on %d-th core of %d....\n", ci_rank, ci_numProce);
	o_modeldata.readdata(fin);
	pdsolve o_sol(o_modeldata,ci_rank,ci_numProce);
	//o_modeldata.writeData();
	//===========solving===========================
	double t1, t2;
	if (ci_rank==0)
	{
		cout << "PDLSM-FEM solving........." << endl;
		t1= MPI_Wtime();
	}	
	o_sol.pdfemStaticSolver_CSRformat(o_modeldata);
	if (ci_rank==0)
	{
		double t2 = MPI_Wtime();
		printf("Elapsed time of PDLSM-FEM mpi solving of %s is %f\n", t2 - t1, argv[1]);
		//===========post processing=====================
		printf("writing results.......\n");
		fioFiles o_files;
		o_files.writeResults(o_modeldata, wflag);

		cout << "PDLSM-FEM program finished." << endl;
	}
	
}
