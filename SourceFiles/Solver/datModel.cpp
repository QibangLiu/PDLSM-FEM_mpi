/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */


#include "datModel.h"
#include<unordered_set>
using namespace std;
datModel::datModel()
{
	cop_Gauss = new pdGaussPt*[2];
	cop_Gauss[0] = new pdGaussPt(2);
	cop_Gauss[1] = new pdGaussPt(5);
	//=====inital data=====
	cop2_Block = NULL;
	cop_datLev2 = new dataLev2();
	ci_numCrack = 0;
	ci_numFami = 0;
	cd_NLF = 0.33;
	cd_dt = 1.0;
	ci_numTstep = 1;
	ci_savefrequence = 1;
	cd_gamma = 0.5;
	cd_beta = 0.25;
	cd_mr = 6.0;
	cd_dcf = 1.0;
	//=====initial flages;
	ci_topk = 0;
	ci_solvFlag = -1;
	cb_FENSF = true;
	ci_failFlag = 0;
	cb_lumpedMass = false;
	ci_TESflag = 2;
	cb_vtkBinary = true;
	//==
	cb_Newmark = false;
}

datModel::~datModel()
{
}

void datModel::readdata(ifstream & fin)
{
	double factor;
	string s_dimen;
	fin >> cs_title;
	getline(fin, cs_label);
	getline(fin, cs_label);
	fin >> s_dimen>> ci_proType;
	//fin >> cd_dtf >> ci_numTstep >> ci_savefrequence >>factor;
	factor = 3;
	cop_geomp = new pdGeomP( factor);
	double e, nu, rho, KIc, sigult;
	fin >> e >> nu >> rho >> KIc >> sigult; //read material parameters;
	cop_material = new pdMaterial(e, nu, rho, KIc, sigult);
	if (s_dimen=="3D")
	{
		ci_Numdimen = 3;
	}
	else if (s_dimen=="2D")
	{
		ci_Numdimen = 2;
	}
	else
	{
		printf("Error: The problem dimension is not 3D or 2D\n");
		exit(0);
	}
	
	//read the material point data 
	int Nid;
	double x[3];
	fin >> ci_numNode >> ci_numEle;
	cop2_Node = new pdNode*[ci_numNode];
	cop_datLev2->cdp_sigma = new double[6 * ci_numNode];
	cop_datLev2->cdp_X = new double[3 * ci_numNode];
	for (int i = 0; i < ci_numNode; i++)
	{

		fin >> Nid >> x[0] >> x[1] >> x[2];
		//cout << id << ' ' << xvalue << endl;
		cop2_Node[i] = new pdNode(Nid,x, cop_datLev2);
	}

	//====read the element data ==================
	int Eid, eleNid[8], algoType, eleNidTetra[4];
	cop2_Eles = new pdfemEles *[ci_numEle];
	for (int i = 0; i < ci_numEle; i++)
	{
		fin >> Eid >> algoType;
		if (ci_Numdimen==3)
		{
			for (int j = 0; j < 8; j++)
			{
				fin >> eleNid[j];
			}
			if (eleNid[2]!=eleNid[3]&&eleNid[4]!=eleNid[5])
			{
				//8N brick element
				cop2_Eles[i] = new pdfemEleBrick8N(Eid, eleNid, algoType, cop_datLev2);
			}
			else
			{
				//4N tetrahedron 
				for (int kk = 0; kk < 3; kk++)
				{
					eleNidTetra[kk] = eleNid[kk];
				}
				eleNidTetra[3] = eleNid[4];
				cop2_Eles[i] = new pdfemEleTetra4N(Eid, eleNidTetra, algoType, cop_datLev2);
			}
		}
		else
		{
			//2D ====
			for (int j = 0; j < 4; j++)
			{
				fin >> eleNid[j];
			}
			if (eleNid[3]!=eleNid[2])
			{
				cop2_Eles[i] = new pdfemEleQuad4N(Eid, eleNid, algoType, cop_datLev2);
			}
			else
			{
				cop2_Eles[i] = new pdfemEleTri3N(Eid, eleNid, algoType, cop_datLev2);
			}
			
		}
		
	}
	//=============read PD boundary================
	string s_blank;
	getline(fin, s_blank);
	getline(fin, s_blank);
	int bNid[4], tn=0;
	fin >> ci_numPDBEs;
	cop2_PDBE = new pdPDBEs*[ci_numPDBEs];
	if (ci_Numdimen==3)
	{
		tn = 4;
	}
	else if (ci_Numdimen == 2)
	{
		tn = 2;
	}
	for (int i = 0; i < ci_numPDBEs; i++)
	{
		for (int j = 0; j < tn; j++)
		{
			fin >> bNid[j];
		}
		cop2_PDBE[i] = new pdPDBEs(bNid, tn);
	}
	//=========read the boundry condition=========
	//=====read essensital BCs;
	getline(fin, s_blank);
	getline(fin, s_blank);
	fin >> ci_numEssentialBCs;
	cop2_EssenBCs = new pdfemEssentialBCs * [ci_numEssentialBCs];
	string s_dDof;
	int numNODE, i_dof;
	double val;
	for (int i = 0; i < ci_numEssentialBCs; i++)
	{
		fin >> numNODE >> s_dDof >> val;
		cop2_EssenBCs[i] = new pdfemEssentialBCs(numNODE, s_dDof, val);
	}
	for (int i = 0; i < ci_numEssentialBCs; i++)
	{
		numNODE = cop2_EssenBCs[i]->getNumNODE();
		val = cop2_EssenBCs[i]->getValue();
		i_dof = cop2_EssenBCs[i]->getDOF();
		for (int j = 0; j < numNODE; j++)
		{
			fin >> Nid;
			cop2_EssenBCs[i]->cip_NID[j] = Nid;
			if (cop2_Node[Nid - 1]->op_getDof(i_dof)->b_isActive())
			{
				cop2_Node[Nid - 1]->op_getDof(i_dof)->setValue(val);
				cop2_Node[Nid - 1]->op_getDof(i_dof)->setNotActive();
			}
			else
			{
				printf("ERROR: The DOF of node %d is repeatly defined!\n", Nid);
				exit(0);
			}
		}
	}
	//=========== read Point BC===================
	ci_numPointBCs = 0;//add future;
	cop2_pointBC = new pdPointBC * [ci_numPointBCs];
	double vForce;
	int pBCid, fDof = -1;
	string s_fDof;

	for (int i = 0; i < ci_numPointBCs; i++)
	{
		fin >> pBCid >> Nid >> s_fDof >> vForce;
		if (s_fDof == "FX")
		{
			fDof = 0;
		}
		else if (s_fDof == "FY")
		{
			fDof = 1;
		}
		else if (s_fDof == "FZ")
		{
			fDof = 2;
		}
		else
		{
			cout << "ERROR: input point force DOF is not FX or FY." << endl;
		}
		cop2_pointBC[i] = new pdPointBC(Nid, fDof, vForce);
	}
	//=========read natural BCs=========
	getline(fin, s_blank);
	getline(fin, s_blank);
	fin >> ci_numNaturalBCs;
	cop2_NaturalBC = new pdfemNaturalBCs * [ci_numNaturalBCs];
	int numEle, nbcConNID[4];
	for (int i = 0; i < ci_numNaturalBCs; i++)
	{
		fin >> numEle >> val;
		cop2_NaturalBC[i] = new pdfemNaturalBCs(numEle, val);
	}
	for (int i = 0; i < ci_numNaturalBCs; i++)
	{
		numEle = cop2_NaturalBC[i]->getNumEle();
		for (int j = 0; j < numEle; j++)
		{
			for (int k = 0; k < tn; k++)
			{
				fin >> nbcConNID[k];
			}
			cop2_NaturalBC[i]->initialEles(j, nbcConNID, tn,cop_datLev2);
		}
	}
	//read no fail region
	getline(fin, s_blank);
	getline(fin, s_blank);
	fin >> ci_numNOFAILnode;
	cip_NOFailNode = new int[ci_numNOFAILnode];
	for (int i = 0; i < ci_numNOFAILnode; i++)
	{
		fin >> cip_NOFailNode[i];
	}
	//// ===========read crack data;
	getline(fin, s_blank);
	getline(fin, s_blank);
	fin >> ci_numCrack;
	cdp_crack = new double[ci_numCrack][3][3];
	if (ci_Numdimen==3)
	{
		for (int i = 0; i < ci_numCrack; i++)
		{
			for (int j = 0; j < 3; j++)//point 
			{
				for (int k = 0; k < 3; k++)// x, y,z
				{
					fin >> cdp_crack[i][j][k];
				}
				
			}
		}
	}
	else if (ci_Numdimen==2)
	{
		for (int i = 0; i < ci_numCrack; i++)
		{
			for (int j = 0; j < 2; j++)//point 
			{
				for (int k = 0; k < 3; k++)// x, y,z
				{
					fin >> cdp_crack[i][j][k];
				}
			}
		}
	}

	//fin >> ci_numCrack;
	//cdp2_crack = new double* [ci_numCrack];
	//for (int i = 0; i < ci_numCrack; i++)
	//{
	//	cdp2_crack[i] = new double[4];
	//	for (int j = 0; j < 4; j++)
	//	{
	//		fin >> cdp2_crack[i][j];
	//	}
	//}
	fin.close();
	
	//setReacForcNode();

	
	
}

void datModel::writeData()
{
	ofstream fout("fout.out");		//fout for out file 	
	fout << "*********************PD program*********" << endl;
	fout << endl;
	fout << cs_title << endl;
	fout << endl;

	fout << "********Material parameter*********" << endl;
	fout << endl;
	cop_material->print(fout);
	fout << endl;
	fout << "************Geometry Model data*****" << endl;
	cop_geomp->print(fout);
	fout << endl;


	fout << "*********Node data********" << endl;
	fout << setw(12) << "Node ID" << setw(16) << "x coordinate" 
		<< setw(16) << "y coordinate" << endl;
	
	for (int i = 0; i <ci_numNode; i++)
	{
		cop2_Node[i]->print(fout);
		
	}

	fout << "************************************" << endl;
	fout << "*******Data of Blocks**********" << endl;
	fout << "************************************" << endl;
	fout << setw(5) << "X ID" << setw(5) << "Y ID"
		<< setw(16) << "Number of Ele" <<setw(24) << " ***  ID of Ele" << endl;
	for (int i = 0; i < ci_numBlocks[0]* ci_numBlocks[1]* ci_numBlocks[2]; i++)
	{
		cop2_Block[i]->print(fout);
	}
	fout << "****************************************************" << endl;
	fout << "***************Data of Families******************" << endl;
	fout << "****************************************************" << endl;
	fout << "****************************************************" << endl;
	fout << "******************Node of Families***************" << endl;
	fout << "****************************************************" << endl;
	fout << setw(12) << "Family ID" << setw(13) << "Number of Node"
		<< setw(18) << " ** ID of NODE" <<  endl;
	
	for (int i = 0; i < ci_numFami; i++)
	{
		cop2_FamiOfNode[i]->printNODE(fout);
		
	}
	fout.close();
}



int datModel::getProType() const
{
	return ci_proType;
}

int datModel::getTotnumVaryEssenBC() const
{
	return ci_numVaryEssenBC;
}

int datModel::getTotnumReacForcNode() const
{
	return ci_numReaForceNode;
}

int datModel::getTotnumFami() const
{
	return ci_numFami;
}

void datModel::SetNumFamilies(int numFami)
{
	ci_numFami = numFami;
}

void datModel::allocaMemoryFami()
{
	/*allocate memory for famlilys*/
	cop2_FamiOfNode = new pdFamily * [ci_numFami];
	//double fac = cop_geomp->getFactor();
	//int nd;
	//double delta_k;
	for (int famk = 0; famk < ci_numFami; famk++)
	{
		cop2_FamiOfNode[famk] = new pdFamily();
		cop2_FamiOfNode[famk]->setID(famk + 1);
		/*nd = civ_pdNodeIDX[famk];
		delta_k = fac * pow((cop2_Node[nd]->getvolume()), 1.0 / ci_Numdimen);
		cop2_FamiOfNode[famk]->sethorizon(delta_k);*/

	}
	if (cop2_FamiOfNode == NULL)
	{
		cout << "Error in memory allocation for families." << endl;
		exit(0);
	}
}

void datModel::setEssentialBC(int id,double val)
{
	
	int numNODE = cop2_EssenBCs[id]->getNumNODE();
	cop2_EssenBCs[id]->cd_value = val;
	int i_dof = cop2_EssenBCs[id]->getDOF();
	int Nid;
	for (int j = 0; j < numNODE; j++)
	{
		Nid = cop2_EssenBCs[id]->cip_NID[j];
		cop2_Node[Nid - 1]->op_getDof(i_dof)->setValue(val);
	}
	
}

int datModel::getTotnumEle() const
{
	return ci_numEle;
}

int datModel::getTotnumEssentialBCs() const
{
	return ci_numEssentialBCs;
}

int datModel::getTotnumPointBCs() const
{
	return ci_numPointBCs;
}

int datModel::getTotnumNode() const
{
	return ci_numNode;
}

int datModel::getTotnumNaturalBCs() const
{
	return ci_numNaturalBCs;
}

int datModel::getTotnumPDBEs() const
{
	return ci_numPDBEs;
}

void datModel::getNumOfBLock(int numBlocks[])
{
	for (int i = 0; i < 3; i++)
	{
		numBlocks[i] = ci_numBlocks[i];
	}
}

void datModel::allocaMemoryBLock()
{
	/*allocate memory for block*/
	int totNumBlocks = ci_numBlocks[0] * ci_numBlocks[1] * ci_numBlocks[2];
	cop2_Block = new pdBlock *[totNumBlocks];
	// i=iz*Nx*Ny+iy*Nx+ix;
	//iz=i/(Nx*Ny); iy=(i%(Nx*Ny))/Nx; ix=i%Nx;
	for (int i = 0; i < totNumBlocks; i++)
	{
		cop2_Block[i] = new pdBlock(i + 1);
	}
	if (cop2_Block == NULL)
	{
		cout << "Error in memory allocation for Blocks." << endl;
		exit(0);
	}
}

void datModel::deleteBLOCK()
{
	if (cop2_Block)
	{
		for (int i = 0; i < ci_numBlocks[0] * ci_numBlocks[1] * ci_numBlocks[2]; i++)
		{
			delete cop2_Block[i];
			cop2_Block[i] = NULL;
		}
		delete[]cop2_Block;
		cop2_Block = NULL;
	}
}



pdFamily * datModel::op_getFami(int famk)
{
	return cop2_FamiOfNode[famk];
}

pdNode * datModel::op_getNode(int N)
{
	return cop2_Node[N];
}

pdfemEles * datModel::op_getEles(int ele)
{
	return cop2_Eles[ele];
}

pdPDBEs * datModel::op_getPDBE(int k)
{
	return cop2_PDBE[k];
}







pdMaterial * datModel::op_getmaterial() const
{
	return cop_material;
}

pdGeomP * datModel::op_getGeomP() const
{
	return cop_geomp;
}

pdGaussPt * datModel::op_getGaussPt(int k)
{
	return cop_Gauss[k];
}




pdfemNaturalBCs* datModel::op_getNaturalBC(int nBC)
{
	return cop2_NaturalBC[nBC];
}

pdfemEssentialBCs* datModel::op_getEssenBC(int ebc)
{
	return cop2_EssenBCs[ebc];
}

pdPointBC* datModel::op_getPointBC(int nBC)
{
	return cop2_pointBC[nBC];
}

pdVaryEssentialBCs * datModel::op_getVaryEssenBC(int k)
{
	return cop2_VaryessentialBC[k];
}

pdReactionForceNode* datModel::op_getReacForcNode(int k)
{
	return cop2_ReaForceNode[k];
}

pdBlock* datModel::op_getBlock(int index)
{

	return cop2_Block[index];
}

int datModel::getnumTstep() const
{
	return ci_numTstep;
}

double datModel::getTstep() const
{
	/*double E, delt_min, rho;
	delt_min = cop_geomp->getminDelta();
	E = cop_material->getE();
	rho = cop_material->getrho();
	double dt = cd_dtf * delt_min / sqrt(E / rho);*/
	return cd_dt;
}

int datModel::getSaveFreq() const
{
	return ci_savefrequence;
}

int datModel::getNumCrack() const
{
	return ci_numCrack;
}



void datModel::writeLocalDamage(ofstream & fout)
{
	//cal phi
	int numNodeOfFam, NID_k;
	double tempDam, locaDama;
	for (int k = 0; k < ci_numFami; k++)
	{

		tempDam = 0;
		numNodeOfFam = cop2_FamiOfNode[k]->getNumNode();
		NID_k = cop2_FamiOfNode[k]->getNodeID(0);
		for (int m = 0; m < numNodeOfFam; m++)
		{
			tempDam = tempDam + cop2_FamiOfNode[k]->getbondstate(m);
		}
		locaDama = 1.0 - tempDam / numNodeOfFam;
		
		cop2_Node[NID_k - 1]->setLocalDamage(locaDama);
	}

	

	fout << "*****The local damage of  point*****" << endl;
	fout << std::left << setw(8) << "eID" << setw(15)
		<< "xc" << setw(15) << "yc" <<
		setw(18) << "phi" << endl;
	fout << setiosflags(ios::scientific) << setprecision(5);
	for (int k = 0; k < ci_numNode; k++)
	{
		cop2_Node[k]->printDamage(fout);
	}
}

void datModel::setNumBlocs(int numBlocks[])
{
	for (int i = 0; i < 3; i++)
	{
		ci_numBlocks[i] = numBlocks[i];
	}
}

void datModel::setReacForcNode()
{
	/*====find the elements for these node that need to get node force;
	This function may be parallelize in future;
	binary_search(begin, end, val) find exist or not;
	upper_bound(begin, end , value); upper_bound to find the first one larger than value, return iterator position*/
	unordered_set<int> reaForNID;
	vector<int>::iterator iter;
	int i_temp, numN;
	for (iter = civ_reacForceOfessBCId.begin(); iter!= civ_reacForceOfessBCId.end(); iter++)
	{
		if (*iter > -1 && *iter < ci_numEssentialBCs)
		{
			numN = cop2_EssenBCs[*iter]->getNumNODE();
			for (int i = 0; i < numN; i++)
			{
				i_temp = cop2_EssenBCs[*iter]->cip_NID[i];
				reaForNID.insert(i_temp);
			}
		}
	}
	civ_reaForceNID.assign(reaForNID.begin(), reaForNID.end());
	reaForNID.clear();
	if (!civ_reaForceNID.empty())
	{
		//sort
		sort(civ_reaForceNID.begin(), civ_reaForceNID.end());
		//=========find element;
		int* conNId, reaFNid, nNE;
		vector<int>reaForceEle;
		for (int ele = 0; ele < ci_numEle; ele++)
		{
			nNE = cop2_Eles[ele]->ci_numNodes;
			conNId = new int[nNE];
			cop2_Eles[ele]->getConNid(conNId);
			for (int i = 0; i < nNE; i++)
			{
				if (binary_search(civ_reaForceNID.begin(), civ_reaForceNID.end(), conNId[i]))
				{
					reaForceEle.push_back(ele);
					break;
				}
			}
			delete[] conNId;
			conNId = NULL;
		}
		ci_numReacForceEle = reaForceEle.size();
		cop2_reacForceEle = new pdReacForceEle * [ci_numReacForceEle];
		for (int k = 0; k < ci_numReacForceEle; k++)
		{
			cop2_reacForceEle[k] = new pdReacForceEle(reaForceEle.at(k));
			nNE = cop2_Eles[reaForceEle.at(k)]->ci_numNodes;
			conNId = new int[nNE];
			cop2_Eles[reaForceEle.at(k)]->getConNid(conNId);
			for (int i = 0; i < nNE; i++)
			{
				if (binary_search(civ_reaForceNID.begin(), civ_reaForceNID.end(), conNId[i]))
				{
					cop2_reacForceEle[k]->civ_reaForceEleNidx.push_back(i);
					i_temp = upper_bound(civ_reaForceNID.begin(), civ_reaForceNID.end(), conNId[i]) - civ_reaForceNID.begin() - 1;
					cop2_reacForceEle[k]->civ_reaForceCPNDnodeIDX.push_back(i_temp);

				}
			}
			delete[] conNId;
			conNId = NULL;
		}
	}
	
}



