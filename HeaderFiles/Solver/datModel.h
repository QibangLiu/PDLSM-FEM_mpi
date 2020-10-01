#pragma once

#include"pdVaryEssentialBCs.h"
#include "pdMaterial.h"
#include "pdNode.h"
#include"pdBlock.h"
#include"pdGeomP.h"
#include"pdFamily.h"
#include "pdfemEleQuad4N.h"
#include "pdfemEleBrick8N.h"
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
#include<vector>

#define pi  3.141592653589793
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
	double* op_getcrack(int i);

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


	
	string cs_title;
	int ci_Numdimen;
	vector<int>civ_feID;// finite elements' ID of each core;
	vector<int>civ_pdeID;//PD elements' ID of each core
	vector<int>civ_pdNodeID;// PD node ID of each core
private:
	datModel(const datModel&);// never using copy constructor;
	string cs_label;
	int ci_proType;
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
	int ci_numCrack;
	double **cdp2_crack;

	double cd_dt;
	int ci_numTstep;
	int ci_savefrequence;


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
	//set reaction force node;
	void setReacForcNode();

};

