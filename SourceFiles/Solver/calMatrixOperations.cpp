#include "calMatrixOperations.h"

void calMatrixOperations::matGaussJordanInverse( Matrix  *op_a,  Matrix * op_aInv)
{
	op_aInv->zero();
	int nr = op_a->i_getNumRows();
	int nc = op_a->i_getNumCols();
	for (int i = 0; i <nr; i++)
	{
		op_aInv->setCoeff(i, i, 1);
	}
	for (int i = 0; i < nr; i++)
	{
		if (abs(op_a->d_getCoeff(i, i)) < 1e-15)
		{
			cout << "can't use gaussin meathod" << endl;
			return;
		}
	}
	double temCoe, dFac, temCoeInv;// dFacInv;
	for (int i = 0; i <(nr - 1); i++)// eliminate
	{
		for (int j = (i + 1); j < nr; j++)
		{
			dFac = (op_a->d_getCoeff(j, i)) / (op_a->d_getCoeff(i, i));
			for (int k = 0; k < nr; k++)
			{
				temCoe = op_a->d_getCoeff(j, k) - (op_a->d_getCoeff(i, k)) *dFac;
				op_a->setCoeff(j, k, temCoe);
				temCoeInv = op_aInv->d_getCoeff(j, k) - (op_aInv->d_getCoeff(i, k)) *dFac;
				op_aInv->setCoeff(j, k, temCoeInv);
			}
		}
	}

	//diagonal normalized
	for (int i = 0; i < nr; i++)
	{
		dFac = op_a->d_getCoeff(i, i);
		for (int j = 0; j<nr; j++)
		{
			temCoe = op_a->d_getCoeff(i, j) / dFac;
			op_a->setCoeff(i, j, temCoe);
			temCoeInv = op_aInv->d_getCoeff(i, j) / dFac;
			op_aInv->setCoeff(i, j, temCoeInv);
		}
	}


	for (int i = nr - 1; i >= 1; i--)
	{
		for (int j = (i - 1); j >= 0; j--)
		{
			dFac = op_a->d_getCoeff(j, i);
			for (int k = 0; k < nr; k++)
			{
				temCoe = op_a->d_getCoeff(j, k) - (op_a->d_getCoeff(i, k)) *dFac;
				op_a->setCoeff(j, k, temCoe);
				temCoeInv = op_aInv->d_getCoeff(j, k) - (op_aInv->d_getCoeff(i, k)) *dFac;
				op_aInv->setCoeff(j, k, temCoeInv);
			}
		}

	}
}

   void calMatrixOperations::matGaussJordan( Matrix * op_a,  Vector * op_b,  Vector * op_x)
{
	int nc = op_a->i_getNumCols();
	int nr = op_a->i_getNumRows();
	 Matrix op_aInv(nr, nc);
	 Vector op_bb(nc), op_xx(nr);
	op_bb = *op_b;
	matGaussJordanInverse(op_a, &op_aInv);
	matMultiply(&op_aInv, &op_bb, &op_xx);
	*op_x = op_xx;//op_x = *op_xx;
}

   void calMatrixOperations::LupDescomposition(Matrix * op_a, Matrix * L, Matrix * U, int *P)
   {
	   int tmp;
	   double tmp2, u, l,temp;

	   L->zero();
	   U->zero();
	   int DN = op_a->i_getNumRows();
	   int row = 0;
	   for (int i = 0; i < DN; i++)
	   {
		   P[i] = i;
	   }
	   for (int i = 0; i < DN-1; i++)
	   {
		   double pp = 0.0;
		   for (int j = i; j < DN; j++)
		   {
			   if (fabs(op_a->d_getCoeff(j,i))>pp)
			   {
				   pp = fabs(op_a->d_getCoeff(j, i));
				   row = j;
			   }
		   }

		   if (pp==0.0)
		   {
			   cout << "Sigularity." << endl;
			   return;
		   }

		   //change P[i] and P[row]
		   tmp = P[i];
		   P[i] = P[row];
		   P[row] = tmp;

		   tmp2 = 0.0;
		   for (int j = 0; j < DN; j++)
		   {
			   //change A[i][j] and A[row][j]
			   tmp2 = op_a->d_getCoeff(i, j);
			   op_a->setCoeff(i, j, op_a->d_getCoeff(row, j));
			   op_a->setCoeff(row, j, tmp2);
		   }

		   //below same as LU decomposition
		    u = op_a->d_getCoeff(i, i);
		    l = 0.0;
			for (int j = i+1; j < DN; j++)
			{
				l = op_a->d_getCoeff(j, i) / u;
				op_a->setCoeff(j, i, l);
				for (int k = i+1; k < DN; k++)
				{
					temp = op_a->d_getCoeff(j, k) - op_a->d_getCoeff(i, k)*l;
					op_a->setCoeff(j, k, temp);
				}
			}
	   }

	   //construct U and L
	   for (int i = 0; i < DN; i++)
	   {
		   for (int j = 0; j <=i; j++)
		   {
			   if (i!=j)
			   {
				   L->setCoeff(i, j, op_a->d_getCoeff(i, j));
			   }
			   else
			   {
				   L->setCoeff(i, j, 1.0);
			   }
		   }
		   for (int k = i; k < DN; k++)
		   {
			   U->setCoeff(i, k, op_a->d_getCoeff(i, k));
		   }
	   }
   }

   void calMatrixOperations::setSolveMatParaANDSolv(Matrix * o_A, Vector * o_b,Vector *o_x)
   {
	   //set PARDISO solver input data*****
	   //  MKL_INT mkli_n number of equation;
	   //MKL_INT *mkli_ia CSR3 format, ia[i] (i<n) points to the first column index of row i in the array ja;one based
	   //MKL_INT *mkli_ja  CSR3 format, array ja contains column indices of the sparse matrix A
		//doub_A for store non-zero coefficient of o_A;
	   //RHS for store components of o_b;
	   long long int mkli_n;//number of equation;
	   long long int*mkli_ia;//CSR3 format, ia[i] (i<n) points to the first column index of row i in the array ja;
	   long long int*mkli_ja;// CSR3 format, array ja contains column indices of the sparse matrix A
	   double *doub_A;//doub_A for store non-zero coefficient of o_A;
	   double *RHS; //RHS for store components of o_b;
	   double *v_Results;

	   mkli_n = o_b->i_getNumRows();
	   int numROW = mkli_n;
	   v_Results = new double[numROW];
	   RHS = new double[numROW];
	   mkli_ia = new long long int[mkli_n+1];

	   int numNonZeroCoef = 0;
	   
	   vector<int>vec_ja; //store value of ja
	   vector<double>vec_a;//store non-zero coef for o_A;
	   double temp;
	   mkli_ia[0] = 0; //zero-based

	   for (int i = 0; i < mkli_n; i++)
	   {
		   //set RHS cd_FP;
		   RHS[i] = o_b->d_getCoeff(i);

		   for (int j = 0; j < mkli_n; j++)
		   {
			   temp = o_A->d_getCoeff(i, j);
			   if (temp != 0)
			   {
				   //set ja,cr;
				   vec_ja.push_back(j); //zero-based;
				   vec_a.push_back(temp);
				   numNonZeroCoef = numNonZeroCoef + 1;
			   }
		   }
		   //set ia;
		   mkli_ia[i + 1] = numNonZeroCoef; //zero-based;
	   }

	   //set mkli_ja, cd_cr;
	   mkli_ja = new long long int[numNonZeroCoef];
	   doub_A = new double[numNonZeroCoef];

	   for (int i = 0; i < numNonZeroCoef; i++)
	   {
		   mkli_ja[i] = vec_ja[i];
		   doub_A[i] = vec_a[i];
	   }

	   //PARDISO solver
	   PARDISO_64Solver(mkli_n, mkli_ia, mkli_ja, doub_A, RHS, v_Results);

	   
	   // store results to o_x;
	   for (int i = 0; i < numROW; i++)
	   {
		   o_x->setCoeff(i, v_Results[i]);
	   }
	   
	   delete[] doub_A;
	   delete[] RHS;
	   delete[] v_Results;
	   doub_A = NULL;
	   RHS = NULL;
	   v_Results = NULL;
	   
	   delete[] mkli_ja;
	   mkli_ja = NULL;
	   //delete [] mkli_a?
	   delete[] mkli_ia;
	   mkli_ia = NULL;
	
   }

   void calMatrixOperations::PARDISO_64Solver(long long int& mkli_n, long long int* mkli_ia, long long int* mkli_ja, double * doub_A, double * RHS, double * v_Results)
   {
	   long long int mtype = 11;       /* Real unsymmetric matrix */
		// Descriptor of main sparse matrix properties
	   //struct matrix_descr descrA;
	   // Structure with sparse matrix stored in CSR format
	   //sparse_matrix_t       csrA;
	   //sparse_operation_t    transA;

	   long long int nrhs = 1;     /* Number of right hand sides. */
	   /* Internal solver memory pointer pt, */
	   /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	   /* or void *pt[64] should be OK on both architectures */
	   void *pt[64];
	   /* Pardiso control parameters. */
	   long long int iparm[64];
	   long long int maxfct, mnum, phase, error, msglvl;
	   /* Auxiliary variables. */

	   double ddum;          /* Double dummy */
	   long long int idum;         /* Integer dummy. */
   /* -------------------------------------------------------------------- */
   /* .. Setup Pardiso control parameters. */
   /* -------------------------------------------------------------------- */
	   for (int i = 0; i < 64; i++)
	   {
		   iparm[i] = 0;
	   }
	   //iparm[59] = 2;         //out of core mode;
	   iparm[0] = 1;         /* No solver default */
	   iparm[1] = 2;         /* Fill-in reordering from METIS */
	   iparm[3] = 0;         /* No iterative-direct algorithm */
	   iparm[4] = 0;         /* No user fill-in reducing permutation */
	   iparm[5] = 0;         /* Write solution into x */
	   iparm[6] = 0;         /* Not in use */
	   iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	   iparm[8] = 0;         /* Not in use */
	   iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	   iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	   iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	   iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	   iparm[13] = 0;        /* Output: Number of perturbed pivots */
	   iparm[14] = 0;        /* Not in use */
	   iparm[15] = 0;        /* Not in use */
	   iparm[16] = 0;        /* Not in use */
	   iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	   iparm[18] = -1;       /* Output: Mflops for LU factorization */
	   iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	   iparm[34] = 1;		//me: set zero-based, default is one-based;
	   maxfct = 1;           /* Maximum number of numerical factorizations. */
	   mnum = 1;         /* Which factorization to use. */
	   msglvl = 0;           /* Print statistical information  */
	   error = 0;            /* Initialize error flag */
   /* -------------------------------------------------------------------- */
   /* .. Initialize the internal solver memory pointer. This is only */
   /* necessary for the FIRST call of the PARDISO solver. */
   /* -------------------------------------------------------------------- */
	   for (int i = 0; i < 64; i++)
	   {
		   pt[i] = 0;
	   }
	   /* -------------------------------------------------------------------- */
	   /* .. Reordering and Symbolic Factorization. This step also allocates */
	   /* all memory that is necessary for the factorization. */
	   /* -------------------------------------------------------------------- */
	   phase = 11;
	   PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, doub_A, mkli_ia, mkli_ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	   if (error != 0)
	   {
		   printf("\nERROR during symbolic factorization: %lld", error);
		   //printf("\niparm[14]: %d, imparm[15]: %d , imparm[16]: %d\n", iparm[14], iparm[15], iparm[16]);
		   exit(1);
	   }
	  // printf("\nReordering completed ... ");
	  // printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	  // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	   /* -------------------------------------------------------------------- */
	   /* .. Numerical factorization. */
	   /* -------------------------------------------------------------------- */
	   phase = 22;
	   PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, doub_A, mkli_ia, mkli_ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	   if (error != 0)
	   {
		   printf("\nERROR during numerical factorization: %lld", error);
		   //printf("\niparm[14]: %d, imparm[15]: %d , imparm[16]: %d\n", iparm[14], iparm[15], iparm[16]);
		   exit(2);
	   }
	   //printf("\nFactorization completed ... ");
	   /* -------------------------------------------------------------------- */
	   /* .. Back substitution and iterative refinement. */
	   /* -------------------------------------------------------------------- */
	   phase = 33;
	   //descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	   //descrA.mode = SPARSE_FILL_MODE_UPPER;
	   //descrA.diag = SPARSE_DIAG_NON_UNIT;
	   //mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, mkli_n, mkli_n, mkli_ia, mkli_ia + 1, mkli_ja, doub_A);
	   iparm[11] = 0;        //Non Conjugate transposed/transpose solve 
	   //transA = SPARSE_OPERATION_NON_TRANSPOSE;
	   //printf("\n\nSolving system with iparm[11] = %d ...\n", (int)iparm[11]);
	   PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, doub_A, mkli_ia, mkli_ja, &idum, &nrhs, iparm, &msglvl, RHS, v_Results, &error);
	   if (error != 0)
	   {
		   printf("\nERROR during solution: %lld", error);
		  //printf("\niparm[14]: %d, imparm[15]: %d , imparm[16]: %d\n", iparm[14], iparm[15], iparm[16]);
		   exit(3);
	   }



	  // mkl_sparse_destroy(csrA);

	   /* -------------------------------------------------------------------- */
	   /* .. Termination and release of memory. */
	   /* -------------------------------------------------------------------- */
	   phase = -1;           /* Release internal memory. */
	   PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, &ddum, mkli_ia, mkli_ja, &idum, &nrhs,
		   iparm, &msglvl, &ddum, &ddum, &error);
   }

   void calMatrixOperations::cluster_PARDISO_64Solver(long long int& mkli_n, long long int* mkli_ia, 
	   long long int* mkli_ja, double* doub_A, double* RHS, double* v_Results, const int *comm)
   {
	   long long int mtype = 11;       /* Real unsymmetric matrix */
	   // Descriptor of main sparse matrix properties
	  // Structure with sparse matrix stored in CSR format

	   long long int nrhs = 1;     /* Number of right hand sides. */
	   /* Internal solver memory pointer pt, */
	   /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	   /* or void *pt[64] should be OK on both architectures */
	   void* pt[64];
	   /* Pardiso control parameters. */
	   long long int iparm[64];
	   long long int maxfct, mnum, phase, error, msglvl;
	   /* Auxiliary variables. */

	   double ddum;          /* Double dummy */
	   long long int idum;         /* Integer dummy. */
   /* -------------------------------------------------------------------- */
   /* .. Setup Pardiso control parameters. */
   /* -------------------------------------------------------------------- */
	   for (int i = 0; i < 64; i++)
	   {
		   iparm[i] = 0;
	   }
	   //iparm[59] = 2;         //out of core mode;
	   iparm[0] = 1;         /* No solver default */
	   iparm[1] = 2;         /* Fill-in reordering from METIS */
	   iparm[3] = 0;         /* No iterative-direct algorithm */
	   iparm[4] = 0;         /* No user fill-in reducing permutation */
	   iparm[5] = 0;         /* Write solution into x */
	   iparm[6] = 0;         /* Not in use */
	   iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	   iparm[8] = 0;         /* Not in use */
	   iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	   iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	   iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	   iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	   iparm[13] = 0;        /* Output: Number of perturbed pivots */
	   iparm[14] = 0;        /* Not in use */
	   iparm[15] = 0;        /* Not in use */
	   iparm[16] = 0;        /* Not in use */
	   iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	   iparm[18] = -1;       /* Output: Mflops for LU factorization */
	   iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	   iparm[34] = 1;		//me: set zero-based, default is one-based;
	   maxfct = 1;           /* Maximum number of numerical factorizations. */
	   mnum = 1;         /* Which factorization to use. */
	   msglvl = 0;           /* Print statistical information  */
	   error = 0;            /* Initialize error flag */
   /* -------------------------------------------------------------------- */
   /* .. Initialize the internal solver memory pointer. This is only */
   /* necessary for the FIRST call of the PARDISO solver. */
   /* -------------------------------------------------------------------- */
	   for (int i = 0; i < 64; i++)
	   {
		   pt[i] = 0;
	   }

	   /* -------------------------------------------------------------------- */
	   /* .. Reordering and Symbolic Factorization. This step also allocates */
	   /* all memory that is necessary for the factorization. */
	   /* -------------------------------------------------------------------- */
	   phase = 11;
	   cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, doub_A, mkli_ia, mkli_ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, comm, &error);
	   if (error != 0)
	   {
		   printf("\nERROR during symbolic factorization: %lld", error);
		   //printf("\niparm[14]: %d, imparm[15]: %d , imparm[16]: %d\n", iparm[14], iparm[15], iparm[16]);
		   exit(1);
	   }
	   // printf("\nReordering completed ... ");
	   // printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	   // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
		///* -------------------------------------------------------------------- */
		///* .. Numerical factorization. */
		///* -------------------------------------------------------------------- */
	   phase = 22;
	   cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, doub_A, mkli_ia, mkli_ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum,comm, &error);
	   if (error != 0)
	   {
		   printf("\nERROR during numerical factorization: %lld", error);
		   //printf("\niparm[14]: %d, imparm[15]: %d , imparm[16]: %d\n", iparm[14], iparm[15], iparm[16]);
		   exit(2);
	   }
	 //  //printf("\nFactorization completed ... ");
	 //  /* -------------------------------------------------------------------- */
	 //  /* .. Back substitution and iterative refinement. */
	 //  /* -------------------------------------------------------------------- */
	   phase = 33;
	   //descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
	   //descrA.mode = SPARSE_FILL_MODE_UPPER;
	   //descrA.diag = SPARSE_DIAG_NON_UNIT;
	   //mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, mkli_n, mkli_n, mkli_ia, mkli_ia + 1, mkli_ja, doub_A);
	   iparm[11] = 0;        //Non Conjugate transposed/transpose solve 
	   //transA = SPARSE_OPERATION_NON_TRANSPOSE;
	   //printf("\n\nSolving system with iparm[11] = %d ...\n", (int)iparm[11]);
	   cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, doub_A, mkli_ia, mkli_ja, &idum, &nrhs, iparm, &msglvl, RHS, v_Results,comm, &error);
	   if (error != 0)
	   {
		   printf("\nERROR during solution: %lld", error);
		   //printf("\niparm[14]: %d, imparm[15]: %d , imparm[16]: %d\n", iparm[14], iparm[15], iparm[16]);
		   exit(3);
	   }



	   // mkl_sparse_destroy(csrA);

		/* -------------------------------------------------------------------- */
		/* .. Termination and release of memory. */
		/* -------------------------------------------------------------------- */
	   phase = -1;           /* Release internal memory. */
	   cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase,
		   &mkli_n, &ddum, mkli_ia, mkli_ja, &idum, &nrhs,
		   iparm, &msglvl, &ddum, &ddum, comm, &error);
   }

   void calMatrixOperations::LupSolve(Matrix * op_a, Vector * op_b,  Vector * op_x)
   {
	   double temp,temp1;
	   int DN = op_a->i_getNumRows();
	   int *P;
	   P = new int[DN];

	   Matrix *L, *U;
	   L = new Matrix(DN, DN);
	   U = new Matrix(DN, DN);

	   double *y;
	   
	   y = new double[DN];

	   LupDescomposition(op_a, L, U, P);
	   //change
	   for (int i = 0; i < DN; i++)
	   {
		   y[i] = op_b->d_getCoeff(P[i]);
		   for (int j = 0; j < i; j++)
		   {
			   y[i] = y[i] - L->d_getCoeff(i, j)*y[j];
		   }
	   }

	   for (int i = DN-1; i >=0; i--)
	   {
		   op_x->setCoeff(i, y[i]);
		   for (int j = DN-1; j >i; j--)
		   {
			   temp = op_x->d_getCoeff(i) - U->d_getCoeff(i, j)*op_x->d_getCoeff(j);
			   op_x->setCoeff(i, temp);
		   }
		   temp1 = op_x->d_getCoeff(i) / U->d_getCoeff(i, i);
		   op_x->setCoeff(i, temp1);
	   }

	   delete L;
	   delete U;
	   delete[] y;
	   delete[] P;
	   L = NULL;
	   U = NULL;
	   P = NULL;
	   y = NULL;

   }

   void calMatrixOperations::LupSolveInverse(Matrix * op_a, Matrix * op_a_inv)
   {
	   int DN = op_a->i_getNumRows();
	   Matrix *a_mirror;
	   a_mirror = new Matrix(DN, DN);
	 
	   Vector *p_x,*b;
	   p_x = new Vector(DN);
	   b = new Vector(DN);

	   for (int i = 0; i < DN; i++)
	   {
		   for (int j = 0; j < DN; j++)
		   {
			   b->setCoeff(j, 0);
		   }
		   b->setCoeff(i, 1.0);

		   for (int j = 0; j < DN; j++)
		   {
			   for (int k = 0; k < DN; k++)
			   {
				   a_mirror->setCoeff(j, k, op_a->d_getCoeff(j, k));
			   }
		   }


		   LupSolve(a_mirror, b, p_x);
		   for (int j = 0; j < DN; j++)
		   {
			   op_a_inv->setCoeff(j, i, p_x->d_getCoeff(j));
		   }

	   }
	   delete a_mirror, p_x, b;
	   a_mirror = NULL; p_x = NULL; b = NULL;
   }

   void calMatrixOperations::PLUSolve(Matrix * op_a, Vector * op_b, Vector * op_x)
   {
	   /*   The routine solves for X the system of linear equations A*X = B,
			where A is an n-by-n matrix, the columns of matrix B are individual
			right-hand sides, and the columns of X are the corresponding
			solutions.

			 The LU decomposition with partial pivoting and row interchanges is
			 used to factor A as A = P*L*U, where P is a permutation matrix, L
			is unit lower triangular, and U is upper triangular. The factored
			form of A is then used to solve the system of equations A*X = B.
			LAPACKE_dgesv (row-major, high-level) Example Program Results.*/
	   MKL_INT n = op_a->i_getNumCols(), nrhs = 1, lda = n, ldb = 1, info;
	   MKL_INT *ipiv;
	   ipiv = new MKL_INT[n];
	   double *a, *b;
	   a = new double[n*n];
	   b = new double[n];

	   for (int i = 0; i < n; i++)
	   {
		   b[i] = op_b->d_getCoeff(i);
		   for (int j = 0; j < n; j++)
		   {
			   a[i*n + j] = op_a->d_getCoeff(i, j);
		   }
	   }

	   /* Solve the equations A*X = B */
	   info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv,
		   b, ldb);
	   /* Check for the exact singularity */
	   if (info > 0) {
		   printf("The diagonal element of the triangular factor of A,\n");
		   printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
		   printf("the solution could not be computed.\n");
		   exit(1);
	   }


	   //save result;
	   for (int i = 0; i < n; i++)
	   {
		   op_x->setCoeff(i, b[i]);
	   }
	   delete[] a;
	   delete[] b;
	   delete[] ipiv;
	   a = NULL;
	   b = NULL;
	   ipiv = NULL;
   }

   void calMatrixOperations::dSymeEigenV(char jobz,Matrix * op_a, Vector * op_eigValu, Matrix * op_eigVector)
   {
	   // jobz must be 'V' or 'N', V-- get eigV, N-- don't get eigV;
	   MKL_INT n, lda, info;
	   n= op_a->i_getNumRows();
	   lda = n;

	   double *w, *a;
	   w = new double[n];
	   a = new double[n*n];

	   for (int i = 0; i < n; i++)
	   {
		   for (int j = i; j < n; j++)
		   {
			   a[i*n + j] = op_a->d_getCoeff(i, j);
		   }
	   }
	 
	   /* Solve eigenproblem */
	   info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, 'U', n, a, lda, w);
	   /* Check for convergence */
	   if (info > 0) {
		   printf("The algorithm failed to compute eigenvalues.\n");
		   exit(1);
	   }
	   //save results
	   for (int i = 0; i < n; i++)
	   {
		   op_eigValu->setCoeff(i, w[i]);
		   if (jobz=='V')
		   {
			   for (int j = 0; j < n; j++)
			   {
				   op_eigVector->setCoeff(i, j, a[i * lda + j]);
			   }
		   }
		   
	   }

	   /* Free workspace */
	   delete[] w, a;
	   w = NULL; a = NULL;
   }

   

   void calMatrixOperations::PARDISOsolveSparse(Matrix * o_A, Vector * o_b, Vector * o_x)
   {
	  /* printf("\n\n\n");
	   printf("\n****************************************\n");
	   printf("****************************************\n");
	   printf("\n    Solving Sparse Matrix......\n\n");*/
	  
	  
	   setSolveMatParaANDSolv(o_A, o_b, o_x);
	   
	  /* printf("\nSolving Finished. \n");
	   printf("\n****************************************\n");
	   printf("****************************************\n");
	   printf("\n\n\n");*/
   }

   void calMatrixOperations::matTranspose( Matrix * o_a,  Matrix* o_at)
{
	int nr = o_a->i_getNumRows();
	int nc = o_a->i_getNumCols();
	double temp;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			temp = o_a->d_getCoeff(i, j);
			o_at->setCoeff(j, i, temp);
		}
	}
}

   void calMatrixOperations::matMultiply( Matrix * o_a,  Matrix * o_b,  Matrix * o_c)
{
	// note initialize
	int a_nr, a_nc, b_nr, b_nc;
	a_nr = o_a->i_getNumRows();
	a_nc = o_a->i_getNumCols();
	b_nr = o_b->i_getNumRows();
	b_nc = o_b->i_getNumCols();


	for (int i = 0; i < a_nr; i++)
	{

		for (int j = 0; j < b_nc; j++)
		{
			double temp = 0.0;
			for (int k = 0; k < a_nc; k++)
			{
				temp =( o_a->d_getCoeff(i, k))*(o_b->d_getCoeff(k, j)) + temp;
			}
			o_c->setCoeff(i, j, temp);
		}
	}
}

   void calMatrixOperations::matMultiply( Matrix * o_a,  Vector * o_x,  Vector * o_y)
{
	//note initialize
	int a_nr, a_nc;
	a_nr = o_a->i_getNumRows();
	a_nc = o_a->i_getNumCols();

	for (int i = 0; i < a_nr; i++)
	{
		double temp = 0.0;
		for (int j = 0; j < a_nc; j++)
		{
			temp = o_a->d_getCoeff(i, j)*o_x->d_getCoeff(j) + temp;

			o_y->setCoeff(i, temp);
		}
	}
}

   void calMatrixOperations::matAdd( Matrix * o_a,  Matrix * o_b,  Matrix * o_c)
{
	int nr = o_a->i_getNumRows();
	int nc = o_a->i_getNumCols();
	double temp;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			temp = o_a->d_getCoeff(i, j) + o_b->d_getCoeff(i, j);
			o_c->setCoeff(i, j, temp);
		}
	}
}

   void calMatrixOperations::matMinus(Matrix * o_a, Matrix * o_b, Matrix * o_c)
   {
	   int nr = o_a->i_getNumRows();
	   int nc = o_a->i_getNumCols();
	   double temp;
	   for (int i = 0; i < nr; i++)
	   {
		   for (int j = 0; j < nc; j++)
		   {
			   temp = o_a->d_getCoeff(i, j) - o_b->d_getCoeff(i, j);
			   o_c->setCoeff(i, j, temp);
		   }
	   }
   }

   void calMatrixOperations::matMinus(Vector * o_a, Vector * o_b, Vector * o_c)
   {
	   int nr = o_a->i_getNumRows();
	  
	   double temp;
	   for (int i = 0; i < nr; i++)
	   {
		   temp = o_a->d_getCoeff(i) - o_b->d_getCoeff(i);
		   o_c->setCoeff(i,temp);
			 
	   }
   }

   void calMatrixOperations::matAdd(Vector * o_x, Vector * o_y, Vector * o_z)
   {
	   int nr = o_x->i_getNumRows();
	   double temp;
	   for (int i = 0; i < nr; i++)
	   {
		   temp = o_x->d_getCoeff(i) + o_y->d_getCoeff(i);
		   o_z->setCoeff(i, temp);
	   }
   }

   void calMatrixOperations::matAdd(Vector* o_x, double factor, Vector* o_y, Vector* o_z)
   {
	   int nr = o_x->i_getNumRows();
	   double temp;
	   for (int i = 0; i < nr; i++)
	   {
		   temp = o_x->d_getCoeff(i) + factor * (o_y->d_getCoeff(i));
		   o_z->setCoeff(i, temp);
	   }
   }

   void calMatrixOperations::matAdd(Vector * o_x, int index, double value)
   {
	   double temp;
	   temp = o_x->d_getCoeff(index);
	   temp = temp + value;
	   o_x->setCoeff(index, temp);
   }

   void calMatrixOperations::matMultiply( Matrix *o_a, double value,  Matrix *o_b)
{
	int nr = o_a->i_getNumRows();
	int nc = o_a->i_getNumCols();
	double temp;
	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < nc; j++)
		{
			temp = o_a->d_getCoeff(i, j)*value;
			o_b->setCoeff(i, j, temp);
		}
	}
}

   void calMatrixOperations::matMultiply(Vector * o_x, double value, Vector * o_y)
   {
	   int nr = o_x->i_getNumRows();
	   for (int i = 0; i < nr; i++)
	   {
		 o_y->setCoeff(i,o_x->d_getCoeff(i)*value);
	   }
   }

   void calMatrixOperations::matMultiply(Vector * o_x, Vector * o_y, double &c)
   {
	   int	nr = o_x->i_getNumRows();
	   c = 0;
	   for (int i = 0; i < nr; i++)
	   {
		   c = c + (o_x->d_getCoeff(i))*(o_y->d_getCoeff(i));
	   }
   }

   void calMatrixOperations::vecDivede(Vector * o_x, Vector * o_y, Vector * o_z)
   {
	   for (int i = 0; i <o_x->i_getNumRows(); i++)
	   {
		   o_z->setCoeff(i, o_y->d_getCoeff(i) / o_x->d_getCoeff(i));
	   }
   }

   void calMatrixOperations::vecMultiply(Vector * o_x, Vector * o_y, Vector * o_z)
   {
	   for (int i = 0; i <o_x->i_getNumRows() ; i++)
	   {
		   o_z->setCoeff(i, o_x->d_getCoeff(i)*(o_y->d_getCoeff(i)));
	   }
   }
