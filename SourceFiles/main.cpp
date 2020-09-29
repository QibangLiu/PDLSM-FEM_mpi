#include"datModel.h"
#include<iostream>
#include<pdfem_mpi.h>
#include<omp.h>
using namespace std;
//==some global variables
calMatrixOperations matoperat;
pdGaussPt o_globGP(2);
pdGaussPt o_sigGP(2);// this GP is for stress extraplation;
int main(int argc, char* argv[])
{
	omp_set_num_threads(4);
	MPI_Init(&argc, &argv);
	pdfem_mpi* po_pdfem = new pdfem_mpi(argc, argv);
	MPI_Finalize();
	return 0;
}