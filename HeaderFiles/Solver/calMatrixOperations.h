/* calMatrixOperations.h
this class is for matrix operations such as matrix multiplication, inversion

*/

#pragma once

//#include "stdafx.h"
#include "Matrix.h"
#include "Vector.h"
#include<vector>
#include<iostream>
#include<cmath>
#include "mkl_pardiso.h"
#include"mkl_cluster_sparse_solver.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>
#include "mkl_lapacke.h"
using namespace std;

class calMatrixOperations
{
public:
	// global function. it solves [a]{x}={b} using Gauss-Jordan elimination. 
	// When it is done, {x} is the solution (=[a]^(-1) {b})
	// and [a] becomes the inverse of the original [a]
	void matGaussJordanInverse( Matrix *op_a,  Matrix* op_aInv);
	void matGaussJordan( Matrix* op_a,  Vector* op_b,  Vector* op_x);
	

	//Using PLU descomposition to solve Ax=b or inver of matrix A
	void LupSolve(Matrix *op_a, Vector *op_b, Vector* op_x);
	void LupSolveInverse(Matrix *op_a, Matrix *op_a_inv);
	
	//============function LAPACKE_dgesv of MKL===========================
	//===Using PLU descomposition to solve Ax=b 
	void PLUSolve(Matrix *op_a, Vector *op_b, Vector* op_x);
	//===get square real symmetrical matrix eigenvalue and vector
	void dSymeEigenV(char jobz, Matrix *op_a, Vector *op_eigValu, Matrix*op_eigVector);
	//===PARDISO solver==========
	void PARDISOsolveSparse(Matrix * o_A, Vector * o_b, Vector* o_x);
	void PARDISO_64Solver(long long int& mkli_n, long long int* mkli_ia, long long int* mkli_ja, 
		double* doub_A, double* RHS, double* v_Results);
	void cluster_PARDISO_64Solver(long long int& mkli_n, long long int* mkli_ia, long long int* mkli_ja,
		double* doub_A, double* RHS, double* v_Results, const int *comm);
	//===========================================================================
	//transpose matrix a and store the result in matrix at
	void matTranspose( Matrix *o_a,  Matrix*o_at);
	
	//multiply matrix a with the matrix b, store the result in matrix c
	void matMultiply( Matrix* o_a,  Matrix* o_b,  Matrix* o_c);
	//multiply matrix a with vector x, store the result in vector y
	void matMultiply( Matrix* o_a,  Vector* o_x,  Vector* o_y);
	void matMultiply(Matrix *o_a, double value, Matrix*o_b);
	void matMultiply(Vector *o_x, double value, Vector *o_y);
	void matMultiply(Vector *o_x, Vector *o_y, double &c);
	
	void matAdd( Matrix *o_a,  Matrix *o_b,  Matrix *o_c);
	void matAdd(Vector *o_x, Vector*o_y, Vector *o_z);
	void matAdd(Vector* o_x,double factor, Vector* o_y, Vector* o_z);
	void matAdd(Vector *o_x, int index, double value);
	void matMinus(Matrix *o_a, Matrix *o_b, Matrix *o_c);
	void matMinus(Vector *o_a, Vector *o_b, Vector *o_c);

	void vecDivede(Vector *o_x, Vector*o_y, Vector *o_z);//o_y diveded by o_x( corresponding element)
	void vecMultiply(Vector *o_x, Vector*o_y, Vector *o_z);
private:
	//Using LUP descomposition to solve Ax=b or inver of matrix A
	void LupDescomposition(Matrix *op_a, Matrix *L, Matrix*U, int *P);
	
	//PARDISO solver
	void setSolveMatParaANDSolv(Matrix * o_A, Vector * o_b,Vector *o_x);
};

