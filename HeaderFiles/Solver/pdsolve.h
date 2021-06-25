#pragma once
#include"fioFile.h"
#include "datModel.h"
#include"Matrix.h"
#include"Vector.h"
#include"calMatrixOperations.h"
#include<iomanip>
#include<cmath>
#include <stdlib.h>
#include <stdio.h>
#include<mpi.h>
#include<omp.h>
#include<limits>
#include<queue>
//#define pi  3.141592653589793
using namespace std;
struct critStru 
{
	// criterion structure;
	// for top K criterion.
	double sd_value;
	int si_famk, si_m;
	critStru(double val, int famk, int m) { sd_value = val; si_famk = famk; si_m = m; }
	bool operator<(const critStru& s_a) const  //operator <
	{
		return sd_value > s_a.sd_value; //min top heap;
	}
};
class pdsolve
{
public:
	pdsolve(datModel&o_dat, int rank, int numProce);
	~pdsolve();
	
	//====data model setting=============================================
	void setDatModel(datModel& o_dat); // set the data model, must be called before solving;
	void setFEID_PDEID(datModel& o_dat);//set the FE id and PDE id for each cores;
	void findDomainDimen(datModel &o_dat);// get the domain dimension size;
	void calVolumeOfNode(datModel& o_dat);// calculate volume of pd node;
	void setPDNODEandnumFami(datModel& o_dat);// set the node type and memory allocate for family;
	void Setdof_Index(datModel& o_dat);//set each dof's equation position;
	void setDeltaMaxMin(datModel& o_dat);// find out the max and min Delta;
	void setBlockAndFami(datModel& o_dat);// initialize block;
	void initialBondState(datModel& o_dat);
	void setNoFailRegion(datModel& o_dat);
	//============PD algorithem===============================================;
	//====some auxiliary functions
	long long int findCSRIndexOfMat(int rowIndex, int colIndex);
	double calArea(double vec1[], double vec2[]);
	void shapeFunctionQuad4N(double N[], double p, double q);
	void matMathcalNt(Matrix* mNt, double p, double q);
	void matN_trans(Matrix* Nmat, double xN[][3], double p, double q);
	void matMathcalNtN(Matrix* mNtN[],Matrix *mNt[],double xN[][3]);//PDBEs of triangle elements
	void matMathcalNt_trian(Matrix* mNt[]);//PDBEs of triangle elements
	//functions for 2D===
	void shapTens2D(Matrix* A, pdFamily* p_fami, datModel& o_dat);
	void vec_gd2D(double g[], double d[], Matrix* A, pdFamily* p_fami, double xi[], datModel& o_dat);
	void matG2D(Matrix* G, Matrix* A, pdFamily* p_fami, int m, datModel& o_dat);
	void matH2D(Matrix* H, pdFamily* p_fami, datModel& o_dat);
	void matC2D(Matrix* C, pdFamily* p_fami, datModel& o_dat);
	void DispGrad(double DG[], int NID, datModel& o_dat);//for calculate displacement gradient;
	//functions for 3D====
	void shapTens3D(Matrix* A, pdFamily* p_fami, datModel& o_dat);
	void shapTens3D3rd(Matrix* A, pdFamily* p_fami, datModel& o_dat);
	void vec_gd3D(double g[], double d[], Matrix* A, pdFamily* p_fami, double xi[], datModel& o_dat);
	void vec_gd3D3rd(double g[], double d[], Matrix* A, pdFamily* p_fami, double xi[], datModel& o_dat);
	void matG3D(Matrix* G, Matrix* A, pdFamily* p_fami, int m, datModel& o_dat);
	void matH3D(Matrix* H, pdFamily* p_fami, datModel& o_dat);
	void matC3D(Matrix* C, pdFamily* p_fami, datModel& o_dat);

	//functions for 2D && 3D====
	double inflFunc(double xi[], pdFamily* p_fami, datModel& o_dat);
	void assembleInterWorkPD(datModel& o_dat);
	void assemblePDBEwork(datModel& o_dat);
	void assembleMassMatPD(datModel& o_dat); //PD node mass;
	void assembleLumpedMass(datModel& o_dat,int numEqua);
	//====CSR format===
	void assembleInterWorkPD_CSRformat(datModel& o_dat,double *U_N);
	void assemblePDBEwork_CSRformat(datModel& o_dat, double* U_N);
	void assemblePDBEworkQuad_CSRformat(datModel& o_dat, double* U_N);
	void assemblePDBEworkTetrahe_CSRformat(datModel& o_dat, double* U_N);
	void assembleMassMatPD_CSRformat(datModel& o_dat); //PD node mass;
	//=======================================================================
	//========================FEM algorithem=================================
	//=======================================================================
	//=========functions for 2D====
	//=========functions for 3D====
	//====functions for 3D && 2D======
	void assembleSEDbyFEM(datModel&o_dat);  // assemble strain energy density by fem ;
	void assembleElemassMatFEM(datModel&o_dat, ofstream&test); // assemble element mass ;
	void calExternalForce(datModel& o_dat); // calcule equivalent extern nodal force;
	void assembleSEDbyFEM_CSRformat(datModel& o_dat, double* U_N);//CSR ---assemble strain energy density by fem ;
	void calExternalForce_CSRformat(datModel& o_dat);//CSR ---equivalent extern nodal force
	//==========================================================================
	//===========Solvers==========Solvers========Solvers====================
	//==========================================================================
	void pdfemSolver(datModel& o_dat, fioFiles& o_files, char* argv[]);
	void calinternalForce_CSRformat(datModel& o_dat, int numEq, double* U_N);
	void calReacForc(datModel& o_dat, ofstream& fout, int stepN);
	//===================================
	//====static solver==================
	//===================================
	//==full matrix format;
	void pdfemAssembleEquaSys(datModel& o_dat);
	void pdfemStaticSolver(datModel& o_dat);
	//===CSR format
	void setCSRIndexes_gloStiffMat(datModel& o_dat);
	void pdfemAssembleEquasSys_CSRformat(datModel& o_dat, int numEq);
	void pdfemStaticSolver_CSRformat(datModel& o_dat, fioFiles& o_files, char* argv[]);
	//===================================
	//===dynamical solver================
	//===================================
	//==CSR format
	void calAcceleration(int numEq, double* dp_A);
	void timeIntegration(datModel& o_dat,Vector * Vu_n,Vector* Vu_nm1,Vector* Vu_np1,int numEq);
	void assembleElemassMatFEM_CSRformat(datModel& o_dat);//CSR ---element mass ;
	void storeDisplacementResult(datModel& o_dat, Vector* U);
	void setCSRIndexes_gloMassMat(datModel& o_dat);
	void pdfemDynamicSolver_CSRformat(datModel& o_dat, fioFiles &o_files, char* argv[]);
	//===Newmark's method
	void pdfemDynamicNewmarkSolver_CSRformat(datModel& o_dat, fioFiles& o_files, char* argv[]);
	void pdfemAssembleKN_CSRformat(datModel& o_dat, int numEq, int n);
	void pdfemAssembleKNFR_CSRformat(datModel& o_dat, int numEq, Vector* a_n, Vector* V_n,double *U_N, int n);
	void updateDispVelo(Vector* Vu_n, Vector* Vv_n, Vector* Va_n, Vector* Va_np1, int numEq, datModel& o_dat);
	void printArr(double* V, int numEq, ofstream& fout);
	//===================================
	//===quasi-static solver================
	//===================================
	//==CSR format
	void pdfemQuasiStaticAssembleEquasSys_CSRformat(datModel& o_dat, int numEq, bool AddLoad);
	void pdfemQuasiStaticSolver_CSRformat(datModel& o_dat, fioFiles& o_files);

	//==failure criterion===========================
	double failureProcess(datModel& o_dat, int Tk, bool& addLoad,ofstream&cracPath);
	double failureCriterion_stretch(datModel& o_dat, int Tk, bool& addLoad); // add 2d Sc later;
	double failureCriterion_stress(datModel& o_dat, int Tk, bool& addLoad);
	double failureCriterion_TopK_Stress(datModel& o_dat,int Tk,bool &addLoad);// TopK is the largest k elements in a set.
	bool Finally_TopK_andUpdate_bondStaus(datModel& o_dat, int Tk, priority_queue<critStru>& TopK, int flag);
	//===maximum principal stress;
	double failureCriterion_maxPriSig(datModel& o_dat, int Tk, bool& addLoad);
	//==maximum circumferential tensile stress:Keq, KI, KII=========
	bool b_cracPropag_qusiaStatic_BYKeq(datModel& o_dat, ofstream& cracPath, double& KEQ);
	void SIFsAndPropaDire(double& KI, double& KII, double& theta_c, double R, double m_R, datModel& o_dat, int ck, ofstream& fout);
	void Jintegrand(double& I_m1, double& I_m2, double& J_s1, int NID1, int NID2, double norV[], int ck, datModel& o_dat, ofstream& fout);
	//handle damage;
	void updateBondstate(double xN[][3], datModel& o_dat);
	void updateBondState_CirclePlane(double xc[], double norm[], double maxDelta, datModel& o_dat);
	bool intersection(double L1X1[], double L1X2[], double L2X1[], double L2X2[]);//judge two line segments intersect or not;
	bool segmentPlaneIntersection(double xp1[], double xp2[], double xN[][3]);
	bool segCirclePlanIntersect(double xp1[], double xp2[], double xc[], double norm[], double Del);
	void calLocalDamage(datModel& o_dat);
	//==vary ebc=================;
	void setVaryEssentialBC(datModel&o_dat);
	void setVaryNaturalBC(datModel& o_dat);
	void set_step_DispBC_qusia_static(datModel& o_dat);
	void resetDispBC(datModel& o_dat, double Multip);



	
	void PDsolve(datModel&o_dat, ofstream&test);
	
	//quasi-static solve;
	void PDsolveQusiaStatic(datModel& o_dat, ofstream& test);
	
	
	void DispGrad(double DG[], int NID, datModel& o_dat, ofstream& test);//for calculate displacement gradient;

	void cracPropag_qusiaStatic(datModel& o_dat, ofstream& test, int solType, double& maxeigv);
	void SIFsAndPropaDire_GaussIntegr(double& KI, double& KII, double& theta_c, double R, datModel& o_dat, double* cracTip, ofstream& test);
	void Jintegr_GaussIntegr(double& I_m1, double& I_m2, int eleID, double* dp_q, double* cracTip, datModel& o_dat, ofstream& test);
	void DispGradAndStressAndPqPx_FEM(double DG[], double sigma[], double PqPx[2], double* dp_q, int ele, double p, double q, datModel& o_dat, ofstream& test);
	void judgeBondState(datModel &o_dat, ofstream&test, int solType,double &maxeigv);
	



	//=====================post-processing==================================;
	//===stresses===;
	void calPDNodeStresses(datModel& o_dat, int* count);
	void calFEMNodeStresses_EXP(datModel& o_dat, int* count);
	void calFEMNodeStresses_LSM(datModel& o_dat, int* count);
	void calGlobalNodeStresses(datModel &o_dat);
	void postProcessing(datModel&o_dat, ofstream&test);

	
private:
	pdsolve();
	//===MPI rank, number of processor;
	int ci_rank, ci_numProce;
	double cd_blockFac;//block size factor;
	//============================================
	////==============FLAGs ========================
	////===solver;
	//int ci_solvFlag; // 0---dynamic solver; 1--static solver; 2 ---quasi-static solver;
	//// PD node on the interface, interact with node in fem domain or not
	//int ci_PDBN_ITA_flag; //0-----NO, 1-----YES;
	//bool cb_InteralForce;
	////=============END flags=========================
	//============================================
	double cd_beta, cd_gamma;
	//==material constants;
	Matrix* cop_D;
	double cd_lambda;
	double cd_mu;
	//For full matrix solving =====
	Matrix *cop_M;//mass matrix;
	Matrix *cop_Ku; //corresponding active dof,unkown displacement;;
	Vector *cop_F;
	//==========Storing matrix in CSR format=========
	double* cdp_Ku,*cdp_KuGlo;// corresponding active dof,unkown displacement; in CRS format;
	double* cdp_M, * cdp_MGlo;// mass and acceleration in CRS format;
	double* cdp_F,*cdp_FGlo, * cdp_Ug;//RHS force, and global displament; mass matrix
	long long int*cip_ja;//for K and M mat; ja is colume index in CSR format; index of ja may be long long int;
	long long int*cip_ia; //for K and M mat; ia is row index in CSR format; one-based format;
	//========
	vector<double*> cv_cracTips;//  old tip- current tip -new tip;

};