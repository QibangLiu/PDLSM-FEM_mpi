#pragma once
#include<mpi.h>
class pdfem_mpi
{
public:
	pdfem_mpi(int argc,char *argv[]);

public:
	int ci_rank, ci_numProce;
};

