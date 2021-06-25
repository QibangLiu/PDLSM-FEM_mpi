#pragma once

#include"pdVaryEssentialBCs.h"
#include "pdMaterial.h"
#include "pdNode.h"
#include"pdBlock.h"
#include"pdGeomP.h"
#include"pdFamily.h"
#include "pdfemEleQuad4N.h"
#include "pdfemEleBrick8N.h"
#include"pdfemEleTetra4N.h"
#include<iostream>
#include<fstream>
#include<string>
#include<iostream>
#include<vector>
#include"pdPDBEs.h"
#include"pdfemNaturalBCs.h"
#include "pdPointBC.h"
#include "pdGaussPt.h"
#include "pdReactionForceNode.h"
#include"pdfemEssentialBCs.h"
#include"dataLev2.h"
#include"pdReacForceEle.h"
#include<vector>
#include<algorithm>

//#define pi  3.141592653589793
//const double pi = acos(-1.0);

using namespace std;


class datModel
{
public:
	//useful
	datModel();
	~datModel();
	void readdata(ifstream &fin);
	void writeData();
	pdFamily *op_getFami(int famk);
	pdNode *op_getNode(int N);
	pdfemEles * op_getEles(int ele);
	pdPDBEs * op_getPDBE(int k);
	pdfemNaturalBCs *op_getNaturalBC(int nBC);
	pdfemEssentialBCs* op_getEssenBC(int ebc);
	pdPointBC* op_getPointBC(int nBC);
	pdVaryEssentialBCs *op_getVaryEssenBC(int k);
	pdReactionForceNode* op_getReacForcNode(int k);
	pdBlock* op_getBlock(int index);
	pdMaterial* op_getmaterial()const;
	pdGeomP* op_getGeomP()const;
	pdGaussPt* op_getGaussPt(int k);

	int getnumTstep()const;
	double getTstep()const;
	int getSaveFreq()const;
	int getNumCrack()const;

	int getProType()const;
	int getTotnumVaryEssenBC()const;
	int getTotnumReacForcNode()const;
	int getTotnumEle()const;
	int getTotnumEssentialBCs()const;
	int getTotnumPointBCs()const;
	int getTotnumNode()const;
	int getTotnumNaturalBCs()const;
	int getTotnumPDBEs()const;
	

	void writeLocalDamage(ofstream&fout);
	//===blocks===;
	void getNumOfBLock(int numBlocks[]);
	void allocaMemoryBLock();
	void deleteBLOCK();
	void setNumBlocs(int numBlocks[]);
	//==families====
	int getTotnumFami()const;
	void SetNumFamilies(int numFami);
	void allocaMemoryFami();

	void setEssentialBC(int id, double val);

	dataLev2* cop_datLev2;
	string cs_title;
	string cs_fileName;
	int ci_Numdimen;
	int ci_eleType;//element Type;
	vector<int>civ_feIDX;// finite elements' Index ;//for SED assembling and FE stress; 
						//or pure FE node volume, or family setting, or FE mass
	vector<int>civ_pdeIDX;//PD elements' Index; for PD node volume only; cleared after volume calculated;
	vector<int>civ_pdNodeIDX;// PD node Index; for family setting and Max min delta;initial in set PD node function;
	//====for reaction force=========================
	vector<int>civ_reacForceOfessBCId;
	vector<int>civ_reaForceNID;
	pdReacForceEle** cop2_reacForceEle;
	int ci_numReacForceEle;
	void setReacForcNode();
	//====================================
	//==crack===========================
	int ci_numCrack;
	double (*cdp_crack)[3][3];
	//==solving setting=================;
	double cd_dt;
	int ci_numTstep;
	int ci_savefrequence;
	double cd_NLF;//non-local factor;
	double cd_gamma, cd_beta;
	//============================================
	//==============FLAGs ========================
	//===solver;
	int ci_solvFlag; // 0---dynamic solver; 1--static solver; 2 ---quasi-static solver;
	// PD node on the interface, interact with node in fem domain or not
	int ci_PDBN_ITA_flag; //0-----NO, 1-----YES;
	/*====failure criterion flags:
	0--critical stretch; 
	1--equivalent stress; 
	2--maximum circumferential tensile stress;
	3-- maximum principal stress*/
	int ci_failFlag;
	bool cb_lumpedMass;
	int ci_TESflag;//2: 2nd TES, 3:3rd TES
	bool cb_Newmark;// nNewmark's method.
	int ci_proType;
	int ci_topk;
	bool cb_vtkBinary;
	//=============END flags=========================

	//=========NO fail===========
	int ci_numNOFAILnode;
	int *cip_NOFailNode;
private:
	datModel(const datModel&);// never using copy constructor;
	string cs_label;
	
	
	int ci_numVaryEssenBC;
	int ci_numNode;//number of material points
	int ci_numEle;
	int ci_numPDBEs;
	int ci_numNaturalBCs;//traction
	int ci_numEssentialBCs;// Displacement
	int ci_numReaForceNode;
	int ci_numPointBCs;// node force;
	int ci_numBlocks[3];//total number of  block in X direction
	int ci_numFami;




	pdMaterial *cop_material;//store material data
	pdGeomP *cop_geomp;//store geometry data
	pdNode **cop2_Node;//store MP data
	pdfemEles** cop2_Eles;//store element data
	pdfemNaturalBCs **cop2_NaturalBC;
	pdBlock **cop2_Block;//store block data
	pdFamily **cop2_FamiOfNode;//store family data
	pdGaussPt **cop_Gauss;
	pdfemEssentialBCs** cop2_EssenBCs;
	pdVaryEssentialBCs **cop2_VaryessentialBC;
	pdPointBC** cop2_pointBC;
	pdReactionForceNode** cop2_ReaForceNode;
	pdPDBEs **cop2_PDBE;
	

	
	// set famlily;
	

};

