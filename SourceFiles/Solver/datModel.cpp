#include "datModel.h"
using namespace std;
datModel::datModel()
{
	cop_Gauss = new pdGaussPt*[2];
	cop_Gauss[0] = new pdGaussPt(2);
	cop_Gauss[1] = new pdGaussPt(5);
	//=====inital data=====
	cop2_Block = NULL;
	cop_datLev2 = new dataLev2();
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
	fin >> cd_dt >> ci_numTstep >> ci_savefrequence >>factor;
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
	fin >> ci_numNode >> ci_numEle >> ci_eleType;
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
	int Eid, eleNid[8], algoType;
	cop2_Eles = new pdfemEles *[ci_numEle];
	for (int i = 0; i < ci_numEle; i++)
	{
		fin >> Eid >> algoType;
		if (ci_eleType==12)//8N brick element
		{
			for (int j = 0; j < 8; j++)
			{
				fin >> eleNid[j];
			}
			cop2_Eles[i] = new pdfemEleBrick8N(Eid, eleNid, algoType, cop_datLev2);
		}
		else if (ci_eleType==10)
		{
			for (int j = 0; j < 4; j++)
			{
				fin >> eleNid[j];
			}
			cop2_Eles[i] = new pdfemEleTetra4N(Eid, eleNid, algoType, cop_datLev2);
		}
		else
		{
			printf("ERROR: element type is wrong.\n");
			exit(0);
		}
		
	}
	//=============read PD boundary================
	string s_blank;
	getline(fin, s_blank);
	getline(fin, s_blank);
	int bNid[4];
	fin >> ci_numPDBEs;
	cop2_PDBE = new pdPDBEs*[ci_numPDBEs];
	for (int i = 0; i < ci_numPDBEs; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			fin >> bNid[j];
		}
		//fin >> normDire;
		cop2_PDBE[i] = new pdPDBEs(bNid, 4);
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
			cop2_Node[Nid - 1]->op_getDof(i_dof)->setValue(val);
			cop2_Node[Nid - 1]->op_getDof(i_dof)->setNotActive();
		}
	}
	//=========== read Point BC===================
	ci_numPointBCs = 0;
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
			for (int k = 0; k < 4; k++)
			{
				fin >> nbcConNID[k];
			}
			cop2_NaturalBC[i]->initialEles(j, nbcConNID, 4,cop_datLev2);
		}
	}
	//// ===========read crack data;
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
	double fac = cop_geomp->getFactor();
	int nd;
	double delta_k;
	for (int famk = 0; famk < ci_numFami; famk++)
	{
		cop2_FamiOfNode[famk] = new pdFamily();
		cop2_FamiOfNode[famk]->setID(famk + 1);
		nd = civ_pdNodeIDX[famk];
		delta_k = fac * pow((cop2_Node[nd]->getvolume()), 1.0 / ci_Numdimen);
		cop2_FamiOfNode[famk]->sethorizon(delta_k);

	}
	if (cop2_FamiOfNode == NULL)
	{
		cout << "Error in memory allocation for families." << endl;
		exit(0);
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

double* datModel::op_getcrack(int i)
{
	return cdp2_crack[i];
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

//void datModel::setReacForcNode()
//{
//	int dire;
//	cout << "please input reaction force direction (0/1)" << endl;
//	cin >> dire;
//	if (dire!=0&&dire!=1)
//	{
//		cout << "Error: please input 0 or 1" << endl;
//	}
//
//	int conNId[4], reaFNid;
//	int algoType;
//	for (int ele = 0; ele < ci_numEle; ele++)
//	{
//		algoType = cop2_Ele4N[ele]->getAlgoType();
//		if (algoType==2)
//		{
//			cop2_Ele4N[ele]->getConNid(conNId);
//
//
//			for (int k = 0; k < ci_numReaForceNode; k++)
//			{
//				reaFNid = cop2_ReaForceNode[k]->getNid();
//				for (int i = 0; i < 4; i++)
//				{
//					if (conNId[i]== reaFNid)
//					{
//						cop2_ReaForceNode[k]->putEleInside(ele + 1,2*i+dire);
//						break;
//					}
//				}
//				
//			}
//			
//		}
//	}
//}



