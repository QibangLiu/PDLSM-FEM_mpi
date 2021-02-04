#include "pdsolve.h"
#include<unordered_set>
#include<algorithm>
extern calMatrixOperations matoperat;
extern pdGaussPt o_globGP;
pdsolve::pdsolve(datModel & o_dat,int rank,int numProce)
{
	ci_rank = rank;
	ci_numProce = numProce;
	//===initial the pointers============
	cop_M = NULL;
	cop_Ku = NULL; cop_F = NULL;
	cip_ia = NULL, cip_ja = NULL, cdp_F = NULL, cdp_Ku = NULL, cdp_Ug = NULL, cdp_M = NULL;
	cdp_FGlo = NULL, cdp_KuGlo = NULL;
	cb_InteralForce = false;
	//===set data model ==================================
	
	setDatModel(o_dat);
	//============================D matrix=================
	int prob_type = o_dat.getProType();
	//D matrix is material stiffness with size 3*3
	double E, nu;
	E = o_dat.op_getmaterial()->getE();
	nu = o_dat.op_getmaterial()->getnu();
	cd_mu = 0.5*E / (1 + nu);
	cd_lambda = 0;
	if (o_dat.ci_Numdimen==2)
	{
		cop_D = new Matrix(3, 3);
		cop_D->zero();
		if (prob_type == 1)
		{
			//plane stress
			cd_lambda = E * nu / (1 - nu * nu); // this lambda is a equivalent value;
			cop_D->setCoeff(0, 0, E / (1 - nu * nu));
			cop_D->setCoeff(0, 1, E * nu / (1 - nu * nu));
			cop_D->setCoeff(1, 0, E * nu / (1 - nu * nu));
			cop_D->setCoeff(1, 1, E / (1 - nu * nu));
			cop_D->setCoeff(2, 2, 0.5 * E / (1 + nu));
		}
		else if (prob_type == 2)
		{
			//plane strain;
			cd_lambda = E * nu / ((1 - 2 * nu) * (1 + nu));
			double multp = E / ((1 - 2 * nu) * (1 + nu));
			cop_D->setCoeff(0, 0, multp * (1 - nu));
			cop_D->setCoeff(0, 1, multp * nu);
			cop_D->setCoeff(1, 0, multp * nu);
			cop_D->setCoeff(1, 1, multp * (1 - nu));
			cop_D->setCoeff(2, 2, multp * 0.5 * (1 - 2 * nu));
		}
		else
		{
			cout << "Error: problem type should be 1 or 2." << endl;
			exit(0);
		}
	}
	else if (o_dat.ci_Numdimen==3)
	{
		//3D;
		cd_lambda = E * nu / ((1 - 2 * nu) * (1 + nu));
		cop_D = new Matrix(6, 6);
		cop_D->zero();
		double multp = E / ((1.0 - 2.0 * nu) * (1.0 + nu));
		cop_D->setCoeff(0, 0, multp * (1 - nu));
		cop_D->setCoeff(0, 1, multp * nu);
		cop_D->setCoeff(0, 2, multp * nu);
		cop_D->setCoeff(1, 0, multp * nu);
		cop_D->setCoeff(1, 1, multp * (1 - nu));
		cop_D->setCoeff(1, 2, multp * nu);
		cop_D->setCoeff(2, 0, multp * nu);
		cop_D->setCoeff(2, 1, multp * nu);
		cop_D->setCoeff(2, 2, multp * (1 - nu));
		cop_D->setCoeff(3, 3, multp * 0.5 * (1 - 2 * nu));
		cop_D->setCoeff(4, 4, multp * 0.5 * (1 - 2 * nu));
		cop_D->setCoeff(5, 5, multp * 0.5 * (1 - 2 * nu));
	}
}

pdsolve::~pdsolve()
{
	delete cop_D;
	cop_D = NULL;
}

void pdsolve::setDatModel(datModel& o_dat)
{
	if (ci_rank==0)
	{
		printf("Setting data model......\n");
	}
	
	// below is order dependent;
	setFEID_PDEID(o_dat);
	Setdof_Index(o_dat);
	setPDNODEandnumFami(o_dat);
	if (o_dat.getTotnumFami()!=0)// if not pure FEM
	{
		findDomainDimen(o_dat);
		calVolumeOfNode(o_dat);
		setDeltaMaxMin(o_dat);
		setBlockAndFami(o_dat);
		initialBondState(o_dat);
		setNoFailRegion(o_dat);
	}
	
	if (ci_rank == 0)
	{
		printf("Finished setting data model.\n");
	}
}

void pdsolve::findDomainDimen(datModel& o_dat)
{
	int numNode, startP, endP;
	numNode = o_dat.getTotnumNode();
	startP = ci_rank * numNode / ci_numProce;
	endP = (ci_rank + 1) * numNode / ci_numProce;

	double lbc[3], rtc[3], d_x[3], lbcGlo[3], rtcGlo[3];
	o_dat.op_getNode(startP)->getcoor(lbc);
	o_dat.op_getNode(startP)->getcoor(rtc);
	o_dat.op_getNode(startP)->getcoor(lbcGlo);
	o_dat.op_getNode(startP)->getcoor(rtcGlo);
	for (int i = startP+ 1; i < endP; i++)
	{
		o_dat.op_getNode(i)->getcoor(d_x);
		for (int j = 0; j < 3; j++)
		{
			if (d_x[j] < lbc[j])
			{
				lbc[j] = d_x[j];
			}
			if (d_x[j] > rtc[j])
			{
				rtc[j] = d_x[j];
			}
		}
	}
	MPI_Reduce(&lbc[0], &lbcGlo[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&rtc[0], &rtcGlo[0], 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (ci_rank == 0)
	{
		for (int i = 0; i < 3; i++)
		{
			lbcGlo[i] = lbcGlo[i] - 1.0E-14 * abs(lbcGlo[i]);
			rtcGlo[i] = rtcGlo[i] + 1.0E-14 * abs(rtcGlo[i]);
		}
	}
	MPI_Bcast(&lbcGlo[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rtcGlo[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	o_dat.op_getGeomP()->setLBC(lbcGlo);
	o_dat.op_getGeomP()->setRTC(rtcGlo);
}

void pdsolve::setPDNODEandnumFami(datModel& o_dat)
{
	//====maybe MPI parallize future; or omp parallize;
	int algoType, * conNID, numNele;
	for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		numNele = o_dat.op_getEles(ele)->getNumNodes();
		if (algoType == 1)
		{
			conNID = new int[numNele];
			o_dat.op_getEles(ele)->getConNid(conNID);
			for (int i = 0; i < numNele; i++)
			{
				o_dat.op_getNode(conNID[i] - 1)->setNodeType(2);
			}
			delete[] conNID; conNID = NULL;
		}
	}
	for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		numNele = o_dat.op_getEles(ele)->getNumNodes();
		if (algoType == 2)
		{
			conNID = new int[numNele];
			o_dat.op_getEles(ele)->getConNid(conNID);
			for (int i = 0; i < numNele; i++)
			{
				if (o_dat.op_getNode(conNID[i] - 1)->getNodeType()==2)
				{
					o_dat.op_getNode(conNID[i] - 1)->setNodeType(1);
				}
				
			}
			delete[] conNID; conNID = NULL;
		}
	}
	/*====get pd node list and the number====*/
	/*Caution: also set the node famyID here*/
	int numFami = 0;
	for (int k = 0; k < o_dat.getTotnumNode(); k++)
	{
		if (o_dat.op_getNode(k)->getNodeType())
		{
			numFami = numFami + 1;
			o_dat.civ_pdNodeIDX.push_back(k);
			o_dat.op_getNode(k)->setFamID(numFami);
		}
	}
	o_dat.SetNumFamilies(numFami);
	o_dat.allocaMemoryFami();
}

void pdsolve::Setdof_Index(datModel& o_dat)
{
	// dependency, can not parallize;
	int countEq, numDimen;
	numDimen = o_dat.ci_Numdimen;
	countEq = 0;
	for (int k = 0; k < o_dat.getTotnumNode(); k++)
	{
		for (int i = 0; i < numDimen; i++)
		{
			if (o_dat.op_getNode(k)->op_getDof(i)->b_isActive())
			{
				o_dat.op_getNode(k)->op_getDof(i)->setEqInde(countEq);
				countEq = countEq + 1;
			}
			else
			{
				// if dof is not active,eqindex=-1
				o_dat.op_getNode(k)->op_getDof(i)->setEqInde(-1);
			}
		}
	}
}

void pdsolve::calVolumeOfNode(datModel& o_dat)
{
	
	/*MPI parallized*/
	//====initiallize;
	int numNode = o_dat.getTotnumNode();
	double* loc_V = new double[numNode];
	double* glo_V = new double[numNode];
	//#pragma omp parallel for
	for (int i = 0; i < numNode; i++)
	{
		loc_V[i] = 0;
		glo_V[i] = 0;
	}
	//=======
	int totNumPDE, startP, endP;
	totNumPDE = o_dat.civ_pdeIDX.size();
	startP = ci_rank * totNumPDE / ci_numProce;
	endP = (ci_rank + 1) * totNumPDE / ci_numProce;
	//======
	int* conNID, numNodeELE, ele;
	double(*xN)[3], VolEle=0;
	for (int pde = startP; pde < endP; pde++)
	{
		
		ele = o_dat.civ_pdeIDX[pde];
		numNodeELE = o_dat.op_getEles(ele)->ci_numNodes;
		conNID = new int[numNodeELE];
		xN = new double[numNodeELE][3];
		o_dat.op_getEles(ele)->getConNid(conNID);
		for (int i = 0; i < numNodeELE; i++)
		{
			o_dat.op_getNode(conNID[i] - 1)->getcoor(xN[i]);
		}
		//Gauss integration to get element volumn
		VolEle = o_dat.op_getEles(ele)->eleVolume(xN);
		for (int i = 0; i < numNodeELE; i++)
		{
			//o_dat.op_getNode(conNID[i] - 1)->addvolume(VolEle / numNodeELE);
			loc_V[conNID[i] - 1] = loc_V[conNID[i] - 1] + VolEle / numNodeELE;
		}
		delete[] conNID, delete[] xN;
		conNID = NULL, xN = NULL;
		
	}
	if (o_dat.ci_PDBN_ITA_flag==1)
	{
		int totNumFE;
		totNumFE = o_dat.civ_feIDX.size();
		startP = ci_rank * totNumFE / ci_numProce;
		endP = (ci_rank + 1) * totNumFE / ci_numProce;
		for (int fe = startP; fe < endP; fe++)
		{
			ele= o_dat.civ_feIDX[fe];
			numNodeELE = o_dat.op_getEles(ele)->ci_numNodes;
			conNID = new int[numNodeELE];
			xN = new double[numNodeELE][3];
			o_dat.op_getEles(ele)->getConNid(conNID);
			for (int i = 0; i < numNodeELE; i++)
			{
				o_dat.op_getNode(conNID[i] - 1)->getcoor(xN[i]);
			}
			//Gauss integration to get element volumn
			VolEle = o_dat.op_getEles(ele)->eleVolume(xN);
			for (int i = 0; i < numNodeELE; i++)
			{
				if (o_dat.op_getNode(conNID[i] - 1)->getNodeType() == 0)
				{
					//only for pure fem node;
					//o_dat.op_getNode(conNID[i] - 1)->addvolume(VolEle / numNodeELE);
					loc_V[conNID[i] - 1] = loc_V[conNID[i] - 1] + VolEle / numNodeELE;
				}
			}
			delete[] conNID, delete[] xN;
			conNID = NULL, xN = NULL;
			
		}
	}
	MPI_Reduce(loc_V, glo_V, numNode, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(glo_V, numNode, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	for (int i = 0; i < numNode; i++)
	{
		o_dat.op_getNode(i)->setVolume(glo_V[i]);
	}
	
	//==release memory;
	delete[] loc_V, delete[] glo_V;
	loc_V = NULL; glo_V = NULL;
	o_dat.civ_pdeIDX.clear();

	//==set horizon of PD NODE;
	double fac = o_dat.op_getGeomP()->getFactor();
	int nd, numDime = o_dat.ci_Numdimen;
	double delta_k;
	int numFami = o_dat.getTotnumFami();
	for (int famk = 0; famk < numFami; famk++)
	{
		nd = o_dat.civ_pdNodeIDX[famk];
		delta_k = fac * pow((o_dat.op_getNode(nd)->getvolume()), 1.0 / numDime);
		o_dat.op_getFami(famk)->sethorizon(delta_k);
	}
}

void pdsolve::setDeltaMaxMin(datModel& o_dat)
{
	double volume_max, volume_min, dv, minDelta, maxDelta, gloVmax, gloVmin;
	//====initialize==
	volume_max = -(numeric_limits<double>::max());
	volume_min = numeric_limits<double>::max();
	gloVmax = volume_max;
	gloVmin = volume_min;
	int totNumPDN, startP, endP, nd;
	totNumPDN = o_dat.civ_pdNodeIDX.size();
	startP = ci_rank * totNumPDN / ci_numProce;
	endP = (ci_rank + 1) * totNumPDN / ci_numProce;
	//======find out max min volumn;
	for (int k = startP; k < endP; k++)
	{
		//if (o_dat.op_getNode(i)->getNodeType())
		{
			nd = o_dat.civ_pdNodeIDX[k];
			dv = o_dat.op_getNode(nd)->getvolume();
			if (dv > volume_max)
			{
				volume_max = dv;
			}
			if (dv < volume_min)
			{
				volume_min = dv;
			}
		}
	}
	MPI_Reduce(&volume_max, &gloVmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&volume_min, &gloVmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	if (ci_rank == 0)
	{
		int numDime = o_dat.ci_Numdimen;
		minDelta = pow(gloVmin, 1.0 / numDime);
		maxDelta = pow(gloVmax, 1.0 / numDime);
	}
	MPI_Bcast(&minDelta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxDelta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	o_dat.op_getGeomP()->setMinDelta(minDelta);
	o_dat.op_getGeomP()->setMaxDelta(maxDelta);
}

void pdsolve::setBlockAndFami(datModel& o_dat)
{
	//=============================set blocks======================
	double fac, maxDelta, blockSize, lbc[3], rtc[3];
	fac = o_dat.op_getGeomP()->getFactor();
	maxDelta = o_dat.op_getGeomP()->getmaxDelta();
	blockSize = 2.5 * fac * maxDelta;
	o_dat.op_getGeomP()->setBlockSize(blockSize);
	o_dat.op_getGeomP()->getlbc(lbc);
	o_dat.op_getGeomP()->getrtc(rtc);
	//==initia numblocks;
	int numBlocks[3];
	for (int i = 0; i < 3; i++)
	{
		numBlocks[i] = (rtc[i] - lbc[i]) / blockSize + 1;
	}
	if (o_dat.ci_Numdimen==2)
	{
		numBlocks[2] = 1;
	}
	o_dat.setNumBlocs(numBlocks);
	o_dat.allocaMemoryBLock();
	//===allocate Nodes into block
	int blockIndex, i_bIndex[3];
	double xnode[3];
	if (o_dat.ci_PDBN_ITA_flag==0)
	{
		for (int k = 0; k < o_dat.getTotnumNode(); k++)
		{
			//only pd node;
			if (o_dat.op_getNode(k)->getNodeType())
			{
				o_dat.op_getNode(k)->getcoor(xnode);
				for (int i = 0; i < 3; i++)
				{
					i_bIndex[i] = (xnode[i] - lbc[i]) / blockSize;
				}
				blockIndex = i_bIndex[0] + i_bIndex[1] * numBlocks[0] +
					i_bIndex[2] * numBlocks[0] * numBlocks[1];
				o_dat.op_getBlock(blockIndex)->putNodeInBlock(k + 1);
			}
		}
	}
	else if (o_dat.ci_PDBN_ITA_flag==1)
	{
		for (int k = 0; k < o_dat.getTotnumNode(); k++)
		{
			// all node;
			o_dat.op_getNode(k)->getcoor(xnode);
			for (int i = 0; i < 3; i++)
			{
				i_bIndex[i] = (xnode[i] - lbc[i]) / blockSize;
			}
			blockIndex = i_bIndex[0] + i_bIndex[1] * numBlocks[0] +
				i_bIndex[2] * numBlocks[0] * numBlocks[1];

			o_dat.op_getBlock(blockIndex)->putNodeInBlock(k + 1);
		}
	}
	else
	{
		printf("ERROR: PD node iteraction flag must be 1 or 1.\n");
		exit(0);
	}
	
	//==allocate element into block; only for quasi-static solver;
	int AlgoType, * conNID, numNele;
	double eleCen[3], (*xN)[3], *N;
	if (o_dat.ci_solvFlag==2)
	{
		for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
		{
			AlgoType = o_dat.op_getEles(ele)->getAlgoType();
			if (AlgoType == 1)
			{
				//eleCen[0] = 0; eleCen[1] = 0; eleCen[2] = 0;
				numNele = o_dat.op_getEles(ele)->getNumNodes();
				conNID = new int[numNele];
				xN = new double[numNele][3];
				N = new double[numNele];
				o_dat.op_getEles(ele)->getConNid(conNID);
				//o_dat.op_getEles(ele)->shapeFunction(N, 0, 0, 0);
				for (int n = 0; n < numNele; n++)
				{
					o_dat.op_getNode(conNID[n] - 1)->getcoor(xN[n]);
					/*for (int i = 0; i < 3; i++)
					{
						eleCen[i] = eleCen[i] + xN[n][i] * N[n];
					}*/
				}
				o_dat.op_getEles(ele)->eleCenter(eleCen, xN);
				for (int i = 0; i < 3; i++)
				{
					i_bIndex[i] = (eleCen[i] - lbc[i]) / blockSize;
				}
				blockIndex = i_bIndex[0] + i_bIndex[1] * numBlocks[0] +
					i_bIndex[2] * numBlocks[0] * numBlocks[1];
				o_dat.op_getBlock(blockIndex)->putEleInBlock(ele + 1);
				delete[]conNID, delete[] xN, delete[] N;
				conNID = NULL, xN = NULL, N = NULL;
			}
		}
	}

	//===============set families===============================
	/*===parallize===*/
	int totNumPDN, startP, endP, nd;
	totNumPDN = o_dat.civ_pdNodeIDX.size();
	startP = ci_rank * totNumPDN / ci_numProce;
	endP = (ci_rank + 1) * totNumPDN / ci_numProce;
	//=====
	double delta_k, delta_j;
	double MPx[3], cx[3];//MPx, cx,are coordinate of  center point of node j and k	 respectively
	double dist;//distance between MPx, cx
	int numNodeoB, NodeIDoB, numDimen, maxDel;
	numDimen = o_dat.ci_Numdimen;
	maxDel = o_dat.op_getGeomP()->getmaxDelta();
	//int countFam = 0;
	for (int famk = startP; famk < endP; famk++)
	{
		nd = o_dat.civ_pdNodeIDX[famk];
		//delta_k = fac * pow((o_dat.op_getNode(nd)->getvolume()), 1.0 / numDimen);
		delta_k = o_dat.op_getFami(famk)->gethorizon();
		o_dat.op_getNode(nd)->getcoor(cx);
		o_dat.op_getFami(famk)->putNodeIntoFami(nd + 1);
		//set fam ID of node; moved to setPDNODEandnumFami()function;
		//o_dat.op_getNode(nd)->setFamID(famk + 1);
		for (int i = 0; i < 3; i++)
		{
			i_bIndex[i] = (cx[i] - lbc[i]) / blockSize;
		}
		for (int xIdex = i_bIndex[0] - 1; xIdex < i_bIndex[0] + 2; xIdex++)
		{
			for (int yIdex = i_bIndex[1] - 1; yIdex < i_bIndex[1] + 2; yIdex++)
			{
				for (int zIdex = i_bIndex[2] - 1; zIdex < i_bIndex[2] + 2; zIdex++)
				{
					if (xIdex >= 0 && xIdex < (numBlocks[0]) &&
						yIdex >= 0 && yIdex < (numBlocks[1]) &&
						zIdex >= 0 && zIdex < (numBlocks[2]))
					{
						blockIndex = xIdex + yIdex * numBlocks[0] +
							zIdex * numBlocks[0] * numBlocks[1];
						numNodeoB = o_dat.op_getBlock(blockIndex)->getNumNodeoB();
						for (int ii = 0; ii < numNodeoB; ii++)
						{
							NodeIDoB = o_dat.op_getBlock(blockIndex)->getNodeoB(ii);
							o_dat.op_getNode(NodeIDoB - 1)->getcoor(MPx);
							delta_j = fac * pow((o_dat.op_getNode(NodeIDoB - 1)->getvolume()), 1.0 / numDimen);
							dist = sqrt((MPx[0] - cx[0]) * (MPx[0] - cx[0]) + (MPx[1] - cx[1]) * (MPx[1] - cx[1])
								+ (MPx[2] - cx[2]) * (MPx[2] - cx[2]));
							if (dist > maxDel * 1.0E-16 &&( dist < (delta_k + maxDel * 1.0E-15)
								 || dist < (delta_j + maxDel * 1.0E-15)))
							{
								o_dat.op_getFami(famk)->putNodeIntoFami(NodeIDoB);
							}
						}
					}
				}
			}
		}
		
	}

	//==== sent and receive data;
	int LocSP, LocEP;
	int* numNDofFam = new int[totNumPDN];
	for (int rank = 0; rank < ci_numProce; rank++)
	{
		LocSP = rank * totNumPDN / ci_numProce;
		LocEP = (rank + 1) * totNumPDN / ci_numProce;
		for (int famk = LocSP; famk < LocEP; famk++)
		{
			//==== sent and receive data;
			numNDofFam[famk] = o_dat.op_getFami(famk)->civ_NID.size();
			MPI_Bcast(&(numNDofFam[famk]), 1, MPI_INT, rank, MPI_COMM_WORLD);
			o_dat.op_getFami(famk)->civ_NID.resize(numNDofFam[famk]);
			MPI_Bcast(&(o_dat.op_getFami(famk)->civ_NID[0]), numNDofFam[famk],
				MPI_INT, rank, MPI_COMM_WORLD);
		}
	}
	//release memory;
	delete[]numNDofFam; numNDofFam = NULL;
	o_dat.civ_pdNodeIDX.clear();
	//initial bond state;
	for (int k = 0; k < o_dat.getTotnumFami(); k++)
	{
		o_dat.op_getFami(k)->initialBondState();
	}
}

void pdsolve::setFEID_PDEID(datModel& o_dat)
{
	int totNumEle, algoType;
	totNumEle = o_dat.getTotnumEle();
	for (int ele = 0; ele < totNumEle; ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		if (algoType == 2)
		{
			o_dat.civ_feIDX.push_back(ele);
		}
		else if (algoType==1)
		{
			o_dat.civ_pdeIDX.push_back(ele);
		}
		else
		{
			printf("Error: Element algorithem type must be 1 or 2.\n");
			exit(0);
		}
	}
}

void pdsolve::setCSRIndexes_gloStiffMat(datModel& o_dat)
{
	//This function is set indexes for global stiffness matrix in CSR format;
	// output is ia and ja;

	//allocate ia zero-based format;
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode, numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = totNumNode * numDime - numPreDof;
	cip_ia = new long long int[numEq + 1];
	for (int i = 0; i < numEq+1; i++)
	{
		//zero-based, ia[0]=0 always.
		cip_ia[i] = 0;
	}
	//============ find interactions of each node=========;
	vector<int>* inteNode;
	inteNode = new vector<int>[totNumNode];
	//====PD family===;
	int NID_k;
	vector<int> famV;
	for (int famk = 0; famk < o_dat.getTotnumFami(); famk++)
	{
		NID_k = o_dat.op_getFami(famk)->getNodeID(0);
		famV = o_dat.op_getFami(famk)->GetvecNid();
		inteNode[NID_k - 1].assign(famV.begin(), famV.end());
	}
	famV.clear();
	//==PDBE interactions==
	int numNele, * conNID, famID, numNodeOfFam, NID_m;
	for (int pdbe = 0; pdbe < o_dat.getTotnumPDBEs(); pdbe++)
	{
		numNele = o_dat.op_getPDBE(pdbe)->i_getNumnode();
		conNID = new int[numNele];
		o_dat.op_getPDBE(pdbe)->getNodeID(conNID);
		for (int nd = 0; nd < numNele; nd++)
		{
			famID = o_dat.op_getNode(conNID[nd] - 1)->getFamID();
			numNodeOfFam = o_dat.op_getFami(famID - 1)->getNumNode();
			for (int m = 0; m < numNodeOfFam; m++)
			{
				NID_m = o_dat.op_getFami(famID - 1)->getNodeID(m);
				for (int i = 0; i < numNele; i++)
				{
					if (i!=nd)
					{
						inteNode[conNID[i] - 1].push_back(NID_m);
					}
				}
			}
		}
		delete[]conNID;
		conNID = NULL;
	}
	//===fem interactions;
	int algoType;
	for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		if (algoType==2)
		{
			numNele = o_dat.op_getEles(ele)->getNumNodes();
			conNID = new int[numNele];
			o_dat.op_getEles(ele)->getConNid(conNID);
			for (int nd = 0; nd < numNele; nd++)
			{
				for (int i = 0; i < numNele; i++)
				{
					inteNode[conNID[nd] - 1].push_back(conNID[i]);
				}
			}
			delete[]conNID;
			conNID = NULL;
		}
	}
	//=====delete the duplicated nodes and sort them in ascending order;
	unordered_set<int> inteSet;
	for (int nd = 0; nd < totNumNode; nd++)
	{
		for (int NODEID : inteNode[nd])
		{
			inteSet.insert(NODEID);
		}
		inteNode[nd].assign(inteSet.begin(), inteSet.end());
		//sort
		sort(inteNode[nd].begin(), inteNode[nd].end());
		inteSet.clear();
	}
	//============get ia, ja==========================================
	int eqInd_row, eqInd_col;
	vector<int>vec_ja;
	for (int nd = 0; nd < totNumNode; nd++)
	{
		for (int i = 0; i < numDime; i++)
		{
			eqInd_row = o_dat.op_getNode(nd)->op_getDof(i)->i_getEqInde();
			if (eqInd_row !=-1)
			{
				for (int itN = 0; itN < inteNode[nd].size(); itN++)
				{
					for (int j = 0; j < numDime; j++)
					{
						eqInd_col = o_dat.op_getNode(inteNode[nd][itN] - 1)->op_getDof(j)->i_getEqInde();
						if (eqInd_col!=-1)
						{
							vec_ja.push_back(eqInd_col);//zero-based;
							// cip_ia[0]=0, always for zero-based, so eqInd_row+1;
							//count the num of non-zero coeff of each row
							// will be revised later;
							cip_ia[eqInd_row+1] = cip_ia[eqInd_row+1] + 1;
						}
					}
				}
			}
		}
	}
	for (int nd = 0; nd < totNumNode; nd++)
	{
		inteNode[nd].clear();
	}
	delete[] inteNode;
	inteNode = NULL;
	//====get real ia=== caution: zero-based;
	for (int i = 0; i < numEq; i++)
	{
		cip_ia[i + 1] = cip_ia[i] + cip_ia[i+1];
	}
	//====store ja in 1D array===;
	long long int size_ja = vec_ja.size();
	cip_ja = new long long int[size_ja];
	for (int i = 0; i < size_ja; i++)
	{
		cip_ja[i] = vec_ja[i];
	}
	vec_ja.clear();
	//===initializa cdp_ku;
	cdp_Ku = new double[size_ja];
	cdp_KuGlo = new double[size_ja];
	for (long long int i = 0; i < size_ja; i++)
	{
		cdp_Ku[i] = 0;
		cdp_KuGlo[i] = 0;
	}
}

bool pdsolve::segmentPlaneIntersection(double xp1[], double xp2[], double xN[][3])
{
	double norm[3], v1[3], v2[3];
	for (int i = 0; i < 3; i++)
	{
		v1[i] = xN[1][i] - xN[0][i];
		v2[i] = xN[2][i] - xN[0][i];
	}
	norm[0] = v1[1] * v2[2] - v1[2] * v2[1];
	norm[1] = v1[2] * v2[0] - v1[0] * v2[2];
	norm[2] = v1[0] * v2[1] - v1[1] * v2[0];
	double temp1, temp2;
	for (int i = 0; i < 3; i++)
	{
		v1[i] = xp1[i] - xN[0][i];
		v2[i] = xp2[i] - xN[0][i];
	}
	temp1 = v1[0] * norm[0] + v1[1] * norm[1] + v1[2] * norm[2];
	temp2 = v2[0] * norm[0] + v2[1] * norm[1] + v2[2] * norm[2];
	if (temp1*temp2>0)
	{
		return false;
	}
	else if (temp1*temp2==0)
	{
		printf("Warning: node on crack plane\n");
		return true;
	}
	else
	{
		double tempV1[3], tempV2[3];
		for (int i = 0; i < 3; i++)
		{
			tempV1[i] = xN[0][i] - xp1[i];
			tempV2[i] = xp2[i] - xp1[i];
		}
		double d = (tempV1[0] * norm[0] + tempV1[1] * norm[1] + tempV1[2] * norm[2]) /
			(tempV2[0] * norm[0] + tempV2[1] * norm[1] + tempV2[2] * norm[2]);
		double Xisc[3];
		for (int i = 0; i < 3; i++)
		{
			Xisc[i] = xp1[i] + d * tempV2[i];
		}

		double A, A0, A1, A2;
		A = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
		//A0
		for (int i = 0; i < 3; i++)
		{
			v1[i] = xN[1][i] - Xisc[i];
			v2[i] = xN[2][i] - Xisc[i];
		}
		norm[0] = v1[1] * v2[2] - v1[2] * v2[1];
		norm[1] = v1[2] * v2[0] - v1[0] * v2[2];
		norm[2] = v1[0] * v2[1] - v1[1] * v2[0];
		A0 = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
		//A1
		for (int i = 0; i < 3; i++)
		{
			v1[i] = Xisc[i] - xN[0][i];
			v2[i] = xN[2][i] - xN[0][i];
		}
		norm[0] = v1[1] * v2[2] - v1[2] * v2[1];
		norm[1] = v1[2] * v2[0] - v1[0] * v2[2];
		norm[2] = v1[0] * v2[1] - v1[1] * v2[0];
		A1 = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
		//A2
		for (int i = 0; i < 3; i++)
		{
			v1[i] = xN[1][i] - xN[0][i];
			v2[i] = Xisc[i] - xN[0][i];
		}
		norm[0] = v1[1] * v2[2] - v1[2] * v2[1];
		norm[1] = v1[2] * v2[0] - v1[0] * v2[2];
		norm[2] = v1[0] * v2[1] - v1[1] * v2[0];
		A2 = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
		if (abs(A0+A1+A2-A)>A*1.0e-14)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	
}

long long int pdsolve::findCSRIndexOfMat(int rowIndex, int colIndex)
{
	//This function is using binary method to 
	// find the  CRS index of the non-zero coeff Mat[rowindex][colindex];
	// in the asscending ja;
	// colIndex is target; 
	// zero -based;
	long long int low = cip_ia[rowIndex];
	long long int high = cip_ia[rowIndex + 1];
	long long int mid;
	while (low<high)
	{
		mid = (low + high) / 2;
		if (colIndex==cip_ja[mid])
		{
			return mid;
		}
		else if (colIndex<cip_ja[mid])
		{
			high = mid;
		}
		else
		{
			low = mid + 1;
		}
	}
	//==if not find;
	printf("Error: can not build matrix in CSR format\n");
	exit(0);
}

double pdsolve::calArea(double vec1[], double vec2[])
{
	return (0.5 * abs(vec1[0] * vec2[1] - vec1[1] * vec2[0]));
}

void pdsolve::shapeFunctionQuad4N(double N[], double p, double q)
{
	N[0] = 0.25 * (1 - p) * (1 - q);
	N[1] = 0.25 * (1 + p) * (1 - q);
	N[2] = 0.25 * (1 + p) * (1 + q);
	N[3] = 0.25 * (1 - p) * (1 + q);
}

pdsolve::pdsolve()
{
}



void pdsolve::initialBondState(datModel& o_dat)
{
	if (o_dat.ci_numCrack==0)
	{
		return;
	}
	//==crack 
	int numDime = o_dat.ci_Numdimen;
	int numCrac;
	numCrac = o_dat.ci_numCrack;
	double(*crack)[3][3] = o_dat.cdp_crack;
	double Xmin[3], Xmax[3], xk[3], xm[3];
	//===block and family
	double blockSize, lbc[3], temp[3];
	blockSize = o_dat.op_getGeomP()->getBlockSize();
	o_dat.op_getGeomP()->getlbc(lbc);
	int i_bIndexMin[3], i_bIndexMax[3], numBlocks[3], blockIndex, numNodeoB, NodeIDoB;
	int famID, numNodeFami, NID_k, NID_m;
	bool crosCrack;
	o_dat.getNumOfBLock(numBlocks);
	if (numDime==3)
	{
		for (int cidx = 0; cidx < numCrac; cidx++)
		{
			for (int kk = 0; kk < 3; kk++)
			{
				temp[0] = crack[cidx][0][kk];
				temp[1] = crack[cidx][1][kk];
				temp[2] = crack[cidx][2][kk];
				Xmin[kk] = *min_element(temp, temp + 3);
				Xmax[kk] = *max_element(temp, temp + 3);
			}

			for (int i = 0; i < 3; i++)
			{
				i_bIndexMin[i] = (Xmin[i] - lbc[i]) / blockSize;
				i_bIndexMax[i] = (Xmax[i] - lbc[i]) / blockSize;
			}
			for (int xIdex = i_bIndexMin[0] - 1; xIdex < i_bIndexMax[0] + 1; xIdex++)
			{
				for (int yIdex = i_bIndexMin[1] - 1; yIdex < i_bIndexMax[1] + 1; yIdex++)
				{
					for (int zIdex = i_bIndexMin[2] - 1; zIdex < i_bIndexMax[2] + 1; zIdex++)
					{
						if (xIdex >= 0 && xIdex < (numBlocks[0]) &&
							yIdex >= 0 && yIdex < (numBlocks[1]) &&
							zIdex >= 0 && zIdex < (numBlocks[2]))
						{
							blockIndex = xIdex + yIdex * numBlocks[0] +
								zIdex * numBlocks[0] * numBlocks[1];
							numNodeoB = o_dat.op_getBlock(blockIndex)->getNumNodeoB();
							for (int ii = 0; ii < numNodeoB; ii++)
							{
								NodeIDoB = o_dat.op_getBlock(blockIndex)->getNodeoB(ii);
								famID = o_dat.op_getNode(NodeIDoB - 1)->getFamID();
								numNodeFami = o_dat.op_getFami(famID - 1)->getNumNode();
								NID_k = o_dat.op_getFami(famID - 1)->getNodeID(0);
								o_dat.op_getNode(NID_k - 1)->getcoor(xk);
								for (int m = 1; m < numNodeFami; m++)
								{
									NID_m = o_dat.op_getFami(famID - 1)->getNodeID(m);
									o_dat.op_getNode(NID_m - 1)->getcoor(xm);

									crosCrack = segmentPlaneIntersection(xk, xm, crack[cidx]);
									if (crosCrack)
									{
										o_dat.op_getFami(famID - 1)->setbondstate(m, 0);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void pdsolve::setNoFailRegion(datModel& o_dat)
{
	int NID, famID;
	for (int i = 0; i < o_dat.ci_numNOFAILnode; i++)
	{
		NID = o_dat.cip_NOFailNode[i];
		famID = o_dat.op_getNode(NID - 1)->getFamID();
		o_dat.op_getFami(famID - 1)->cb_allowFail = false;
	}
	delete[] o_dat.cip_NOFailNode;
}

double pdsolve::inflFunc( double xi[], pdFamily* p_fami, datModel&o_dat)
{
	double delt = p_fami->gethorizon();
	/*double fac = o_dat.op_getGeomP()->getFactor();
	double delt =fac*o_dat.op_getGeomP()->getmaxDelta();*/
	/*if (ci_rank==0)
	{
		cout << "==" << p_fami->gethorizon() << "\t max=" << delt << endl;
	}*/
	double a = 1/4.0;
	double Aa = a * a;
	double omega = exp(-(xi[0] * xi[0] + xi[1] * xi[1]+ xi[2] * xi[2])
		/ Aa / delt / delt);
	return omega;
}


void pdsolve::shapTens2D(Matrix *A,  pdFamily* p_fami, datModel & o_dat)
{
	//for 2D A matrix with size 5*5;
	A->zero();
	double xi[3], omega, temp, xj[3], xk[3],dv;
	int Nid_m,Nid_k;
	int numNodeOFfami = p_fami->getNumNode();
	Nid_k = p_fami->getNodeID(0);
	o_dat.op_getNode(Nid_k-1)->getcoor(xk);
	for (int j = 1; j < numNodeOFfami; j++)
	{
		Nid_m = p_fami->getNodeID(j);
		o_dat.op_getNode(Nid_m - 1)->getcoor(xj);
		dv = o_dat.op_getNode(Nid_m - 1)->getvolume();
		for (int ii = 0; ii < 3; ii++)
		{
			xi[ii] = xj[ii] - xk[ii];
		}
		omega = inflFunc(xi, p_fami, o_dat);
		//row 0;
		temp = omega*xi[0] * xi[0] *dv;
		A->addCoeff(0, 0, temp);
		temp = omega*dv*xi[0] * xi[1];
		A->addCoeff(0, 1, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[0] * 0.5;
		A->addCoeff(0, 2, temp);
		temp = omega * dv*xi[0] * xi[1] * xi[1] * 0.5;
		A->addCoeff(0, 3, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[1];
		A->addCoeff(0, 4, temp);
		//row 1
		temp = omega * dv*xi[0] * xi[1];
		A->addCoeff(1, 0, temp);
		temp = omega * dv*xi[1] * xi[1];
		A->addCoeff(1, 1, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[1] * 0.5;
		A->addCoeff(1, 2, temp);
		temp = omega * dv*xi[1] * xi[1] * xi[1] * 0.5;
		A->addCoeff(1, 3, temp);
		temp = omega * dv*xi[0] * xi[1] * xi[1];
		A->addCoeff(1, 4, temp);
		// row 2
		temp = omega * dv*xi[0] * xi[0] * xi[0];
		A->addCoeff(2, 0, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[1];
		A->addCoeff(2, 1, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[0] * xi[0] * 0.5;
		A->addCoeff(2, 2, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[1] * xi[1] * 0.5;
		A->addCoeff(2, 3, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[0] * xi[1];
		A->addCoeff(2, 4, temp);
		// row 3
		temp = omega * dv*xi[0] * xi[1] * xi[1];
		A->addCoeff(3, 0, temp);
		temp = omega * dv*xi[1] * xi[1] * xi[1];
		A->addCoeff(3, 1, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[1] * xi[1] * 0.5;
		A->addCoeff(3, 2, temp);
		temp = omega * dv*xi[1] * xi[1] * xi[1] * xi[1] * 0.5;
		A->addCoeff(3, 3, temp);
		temp = omega * dv*xi[0] * xi[1] * xi[1] * xi[1];
		A->addCoeff(3, 4, temp);
		// row 4;
		temp = omega * dv*xi[0] * xi[0] * xi[1];
		A->addCoeff(4, 0, temp);
		temp = omega * dv*xi[0] * xi[1] * xi[1];
		A->addCoeff(4, 1, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[0] * xi[1] * 0.5;
		A->addCoeff(4, 2, temp);
		temp = omega * dv*xi[0] * xi[1] * xi[1] * xi[1] * 0.5;
		A->addCoeff(4, 3, temp);
		temp = omega * dv*xi[0] * xi[0] * xi[1] * xi[1];
		A->addCoeff(4, 4, temp);
	}
}

void pdsolve::shapTens3D(Matrix* A, pdFamily* p_fami, datModel& o_dat)
{
	// matrix A's size is 9*9;
	A->zero();
	double xi[3], omega, xj[3], xk[3], dv;
	int Nid_m, Nid_k;
	int numNodeOFfami = p_fami->getNumNode();
	Nid_k = p_fami->getNodeID(0);
	o_dat.op_getNode(Nid_k - 1)->getcoor(xk);
	for (int j = 1; j < numNodeOFfami; j++)
	{
		Nid_m = p_fami->getNodeID(j);
		o_dat.op_getNode(Nid_m - 1)->getcoor(xj);
		dv = o_dat.op_getNode(Nid_m - 1)->getvolume();
		for (int ii = 0; ii < 3; ii++)
		{
			xi[ii] = xj[ii] - xk[ii];
		}
		omega = inflFunc(xi, p_fami, o_dat);
		/*double temp[9][9] = { xi[0] * xi[0],xi[0] * xi[1],xi[0] * xi[2],0.5 * xi[0] * xi[0] * xi[0],0.5 * xi[0] * xi[1] * xi[1],0.5 * xi[0] * xi[2] * xi[2],xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[2],xi[0] * xi[1] * xi[2],
							xi[0] * xi[1],xi[1] * xi[1],xi[1] * xi[2],0.5 * xi[0] * xi[0] * xi[1],0.5 * xi[1] * xi[1] * xi[1],0.5 * xi[1] * xi[2] * xi[2],xi[0] * xi[1] * xi[1],xi[0] * xi[1] * xi[2],xi[1] * xi[1] * xi[2],
							xi[0] * xi[2],xi[1] * xi[2],xi[2] * xi[2],0.5 * xi[0] * xi[0] * xi[2],0.5 * xi[1] * xi[1] * xi[2],0.5 * xi[2] * xi[2] * xi[2],xi[0] * xi[1] * xi[2],xi[0] * xi[2] * xi[2],xi[1] * xi[2] * xi[2],
							xi[0] * xi[0] * xi[0],xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[2],0.5 * xi[0] * xi[0] * xi[0] * xi[0],0.5 * xi[0] * xi[0] * xi[1] * xi[1],0.5 * xi[0] * xi[0] * xi[2] * xi[2],xi[0] * xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[0] * xi[2],xi[0] * xi[0] * xi[1] * xi[2],
							xi[0] * xi[1] * xi[1],xi[1] * xi[1] * xi[1],xi[1] * xi[1] * xi[2],0.5 * xi[0] * xi[0] * xi[1] * xi[1],0.5 * xi[1] * xi[1] * xi[1] * xi[1],0.5 * xi[1] * xi[1] * xi[2] * xi[2],xi[0] * xi[1] * xi[1] * xi[1],xi[0] * xi[1] * xi[1] * xi[2],xi[1] * xi[1] * xi[1] * xi[2],
							xi[0] * xi[2] * xi[2],xi[1] * xi[2] * xi[2],xi[2] * xi[2] * xi[2],0.5 * xi[0] * xi[0] * xi[2] * xi[2],0.5 * xi[1] * xi[1] * xi[2] * xi[2],0.5 * xi[2] * xi[2] * xi[2] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],xi[0] * xi[2] * xi[2] * xi[2],xi[1] * xi[2] * xi[2] * xi[2],
							xi[0] * xi[0] * xi[1],xi[0] * xi[1] * xi[1],xi[0] * xi[1] * xi[2],0.5 * xi[0] * xi[0] * xi[0] * xi[1],0.5 * xi[0] * xi[1] * xi[1] * xi[1],0.5 * xi[0] * xi[1] * xi[2] * xi[2],xi[0] * xi[0] * xi[1] * xi[1],xi[0] * xi[0] * xi[1] * xi[2],xi[0] * xi[1] * xi[1] * xi[2],
							xi[0] * xi[0] * xi[2],xi[0] * xi[1] * xi[2],xi[0] * xi[2] * xi[2],0.5 * xi[0] * xi[0] * xi[0] * xi[2],0.5 * xi[0] * xi[1] * xi[1] * xi[2],0.5 * xi[0] * xi[2] * xi[2] * xi[2],xi[0] * xi[0] * xi[1] * xi[2],xi[0] * xi[0] * xi[2] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],
							xi[0] * xi[1] * xi[2],xi[1] * xi[1] * xi[2],xi[1] * xi[2] * xi[2],0.5 * xi[0] * xi[0] * xi[1] * xi[2],0.5 * xi[1] * xi[1] * xi[1] * xi[2],0.5 * xi[1] * xi[2] * xi[2] * xi[2],xi[0] * xi[1] * xi[1] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],xi[1] * xi[1] * xi[2] * xi[2]
		};*/

		double temp[9][9] = { xi[0] * xi[0],xi[0] * xi[1],xi[0] * xi[2] ,(xi[0] * xi[0] * xi[0]) * 0.5,(xi[0] * xi[1] * xi[1]) * 0.5,(xi[0] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[1],xi[0] * xi[1] * xi[2],xi[0] * xi[0] * xi[2] ,
		xi[0] * xi[1],xi[1] * xi[1],xi[1] * xi[2], (xi[0] * xi[0] * xi[1]) * 0.5,(xi[1] * xi[1] * xi[1]) * 0.5,(xi[1] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[1],xi[1] * xi[1] * xi[2],xi[0] * xi[1] * xi[2],
			xi[0] * xi[2],xi[1] * xi[2],xi[2] * xi[2], (xi[0] * xi[0] * xi[2]) * 0.5,(xi[1] * xi[1] * xi[2]) * 0.5,(xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[2],xi[1] * xi[2] * xi[2],xi[0] * xi[2] * xi[2],
			xi[0] * xi[0] * xi[0],xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[2],(xi[0] * xi[0] * xi[0] * xi[0]) * 0.5,(xi[0] * xi[0] * xi[1] * xi[1]) * 0.5,(xi[0] * xi[0] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[1] * xi[2],xi[0] * xi[0] * xi[0] * xi[2],
			xi[0] * xi[1] * xi[1],xi[1] * xi[1] * xi[1],xi[1] * xi[1] * xi[2],(xi[0] * xi[0] * xi[1] * xi[1]) * 0.5,(xi[1] * xi[1] * xi[1] * xi[1]) * 0.5,(xi[1] * xi[1] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[1] * xi[1],xi[1] * xi[1] * xi[1] * xi[2],xi[0] * xi[1] * xi[1] * xi[2],
			xi[0] * xi[2] * xi[2],xi[1] * xi[2] * xi[2],xi[2] * xi[2] * xi[2],(xi[0] * xi[0] * xi[2] * xi[2]) * 0.5,(xi[1] * xi[1] * xi[2] * xi[2]) * 0.5,(xi[2] * xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[2] * xi[2],xi[1] * xi[2] * xi[2] * xi[2],xi[0] * xi[2] * xi[2] * xi[2],
			xi[0] * xi[0] * xi[1],xi[0] * xi[1] * xi[1],xi[0] * xi[1] * xi[2],(xi[0] * xi[0] * xi[0] * xi[1]) * 0.5,(xi[0] * xi[1] * xi[1] * xi[1]) * 0.5,(xi[0] * xi[1] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[1] * xi[1],xi[0] * xi[1] * xi[1] * xi[2],xi[0] * xi[0] * xi[1] * xi[2],
			 xi[0] * xi[1] * xi[2],xi[1] * xi[1] * xi[2],xi[1] * xi[2] * xi[2],(xi[0] * xi[0] * xi[1] * xi[2]) * 0.5,(xi[1] * xi[1] * xi[1] * xi[2]) * 0.5,(xi[1] * xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[1] * xi[2],xi[1] * xi[1] * xi[2] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],
			xi[0] * xi[0] * xi[2],xi[0] * xi[1] * xi[2],xi[0] * xi[2] * xi[2],(xi[0] * xi[0] * xi[0] * xi[2]) * 0.5,(xi[0] * xi[1] * xi[1] * xi[2]) * 0.5,(xi[0] * xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[1] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],xi[0] * xi[0] * xi[2] * xi[2],
		};

		for (int row = 0; row < 9; row++)
		{
			for (int col = 0; col < 9; col++)
			{
				A->addCoeff(row, col, omega * dv * temp[row][col]);
			}
		}
	}
}

void pdsolve::shapTens3D3rd(Matrix* A, pdFamily* p_fami, datModel& o_dat)
{
	// matrix A's size is 9*9;
	A->zero();
	double xi[3], omega, xj[3], xk[3], dv;
	int Nid_m, Nid_k;
	int numNodeOFfami = p_fami->getNumNode();
	Nid_k = p_fami->getNodeID(0);
	o_dat.op_getNode(Nid_k - 1)->getcoor(xk);
	for (int j = 1; j < numNodeOFfami; j++)
	{
		Nid_m = p_fami->getNodeID(j);
		o_dat.op_getNode(Nid_m - 1)->getcoor(xj);
		dv = o_dat.op_getNode(Nid_m - 1)->getvolume();
		for (int ii = 0; ii < 3; ii++)
		{
			xi[ii] = xj[ii] - xk[ii];
		}
		omega = inflFunc(xi, p_fami, o_dat);
		double tempA1[9][9] = { xi[0] * xi[0],xi[0] * xi[1],xi[0] * xi[2] ,(xi[0] * xi[0] * xi[0]) * 0.5,(xi[0] * xi[1] * xi[1]) * 0.5,(xi[0] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[1],xi[0] * xi[1] * xi[2],xi[0] * xi[0] * xi[2] ,
		xi[0] * xi[1],xi[1] * xi[1],xi[1] * xi[2], (xi[0] * xi[0] * xi[1]) * 0.5,(xi[1] * xi[1] * xi[1]) * 0.5,(xi[1] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[1],xi[1] * xi[1] * xi[2],xi[0] * xi[1] * xi[2],
			xi[0] * xi[2],xi[1] * xi[2],xi[2] * xi[2], (xi[0] * xi[0] * xi[2]) * 0.5,(xi[1] * xi[1] * xi[2]) * 0.5,(xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[2],xi[1] * xi[2] * xi[2],xi[0] * xi[2] * xi[2],
			xi[0] * xi[0] * xi[0],xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[2],(xi[0] * xi[0] * xi[0] * xi[0]) * 0.5,(xi[0] * xi[0] * xi[1] * xi[1]) * 0.5,(xi[0] * xi[0] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[0] * xi[1],xi[0] * xi[0] * xi[1] * xi[2],xi[0] * xi[0] * xi[0] * xi[2],
			xi[0] * xi[1] * xi[1],xi[1] * xi[1] * xi[1],xi[1] * xi[1] * xi[2],(xi[0] * xi[0] * xi[1] * xi[1]) * 0.5,(xi[1] * xi[1] * xi[1] * xi[1]) * 0.5,(xi[1] * xi[1] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[1] * xi[1],xi[1] * xi[1] * xi[1] * xi[2],xi[0] * xi[1] * xi[1] * xi[2],
			xi[0] * xi[2] * xi[2],xi[1] * xi[2] * xi[2],xi[2] * xi[2] * xi[2],(xi[0] * xi[0] * xi[2] * xi[2]) * 0.5,(xi[1] * xi[1] * xi[2] * xi[2]) * 0.5,(xi[2] * xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[2] * xi[2],xi[1] * xi[2] * xi[2] * xi[2],xi[0] * xi[2] * xi[2] * xi[2],
			xi[0] * xi[0] * xi[1],xi[0] * xi[1] * xi[1],xi[0] * xi[1] * xi[2],(xi[0] * xi[0] * xi[0] * xi[1]) * 0.5,(xi[0] * xi[1] * xi[1] * xi[1]) * 0.5,(xi[0] * xi[1] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[1] * xi[1],xi[0] * xi[1] * xi[1] * xi[2],xi[0] * xi[0] * xi[1] * xi[2],
			 xi[0] * xi[1] * xi[2],xi[1] * xi[1] * xi[2],xi[1] * xi[2] * xi[2],(xi[0] * xi[0] * xi[1] * xi[2]) * 0.5,(xi[1] * xi[1] * xi[1] * xi[2]) * 0.5,(xi[1] * xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[1] * xi[1] * xi[2],xi[1] * xi[1] * xi[2] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],
			xi[0] * xi[0] * xi[2],xi[0] * xi[1] * xi[2],xi[0] * xi[2] * xi[2],(xi[0] * xi[0] * xi[0] * xi[2]) * 0.5,(xi[0] * xi[1] * xi[1] * xi[2]) * 0.5,(xi[0] * xi[2] * xi[2] * xi[2]) * 0.5,xi[0] * xi[0] * xi[1] * xi[2],xi[0] * xi[1] * xi[2] * xi[2],xi[0] * xi[0] * xi[2] * xi[2]
		};

		double tempA2[9][10] = { (xi[0]*xi[0]*xi[0]*xi[0]) / 6,(xi[0] * xi[1]*xi[1]*xi[1]) / 6,(xi[0] * xi[2]*xi[2]*xi[2]) / 6,(xi[0]*xi[0]*xi[0] * xi[1]) / 2,(xi[0]*xi[0] * xi[1]*xi[1]) / 2,(xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[0] * xi[1] * xi[2]*xi[2]) / 2,(xi[0]*xi[0] * xi[2]*xi[2]) / 2,(xi[0]*xi[0]*xi[0] * xi[2]) / 2,xi[0]*xi[0] * xi[1] * xi[2],
	(xi[0]*xi[0]*xi[0] * xi[1]) / 6,(xi[1]*xi[1]*xi[1]*xi[1]) / 6,(xi[1] * xi[2]*xi[2]*xi[2]) / 6,(xi[0]*xi[0] * xi[1]*xi[1]) / 2,(xi[0] * xi[1]*xi[1]*xi[1]) / 2,(xi[1]*xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1] * xi[2]*xi[2]) / 2,(xi[0] * xi[1] * xi[2]*xi[2]) / 2,(xi[0]*xi[0] * xi[1] * xi[2]) / 2,xi[0] * xi[1]*xi[1] * xi[2],
	(xi[0]*xi[0]*xi[0] * xi[2]) / 6,(xi[1]*xi[1]*xi[1] * xi[2]) / 6,(xi[2]*xi[2]*xi[2]*xi[2]) / 6,(xi[0]*xi[0] * xi[1] * xi[2]) / 2,(xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1] * xi[2]*xi[2]) / 2,(xi[1] * xi[2]*xi[2]*xi[2]) / 2,(xi[0] * xi[2]*xi[2]*xi[2]) / 2,(xi[0]*xi[0] * xi[2]*xi[2]) / 2,xi[0] * xi[1] * xi[2]*xi[2],
		(pow(xi[0],5)) / 6,(xi[0]*xi[0] * xi[1]*xi[1]*xi[1]) / 6,(xi[0]*xi[0] * xi[2]*xi[2]*xi[2]) / 6,(xi[0]*xi[0]*xi[0]*xi[0] * xi[1]) / 2,(xi[0]*xi[0]*xi[0] * xi[1]*xi[1]) / 2,(xi[0]*xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[0]*xi[0] * xi[1] * xi[2]*xi[2]) / 2,(xi[0]*xi[0]*xi[0] * xi[2]*xi[2]) / 2,(xi[0]*xi[0]*xi[0]*xi[0] * xi[2]) / 2,xi[0]*xi[0]*xi[0] * xi[1] * xi[2],
	(xi[0]*xi[0]*xi[0] * xi[1]*xi[1]) / 6,(pow(xi[1],5)) / 6,(xi[1]*xi[1] * xi[2]*xi[2]*xi[2]) / 6,(xi[0]*xi[0] * xi[1]*xi[1]*xi[1]) / 2,(xi[1]*xi[1]*xi[1]*xi[1] * xi[0]) / 2,(xi[1]*xi[1]*xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1]*xi[1] * xi[2]*xi[2]) / 2,(xi[1]*xi[1] * xi[0] * xi[2]*xi[2]) / 2,(xi[0]*xi[0] * xi[1]*xi[1] * xi[2]) / 2,xi[1]*xi[1]*xi[1] * xi[0] * xi[2],
	(xi[0]*xi[0]*xi[0] * xi[2]*xi[2]) / 6,(xi[1]*xi[1]*xi[1] * xi[2]*xi[2]) / 6,(pow(xi[2],5)) / 6,(xi[0]*xi[0] * xi[1] * xi[2]*xi[2]) / 2,(xi[1]*xi[1] * xi[0] * xi[2]*xi[2]) / 2,(xi[1]*xi[1] * xi[2]*xi[2]*xi[2]) / 2,(xi[2]*xi[2]*xi[2]*xi[2] * xi[1]) / 2,(xi[2]*xi[2]*xi[2]*xi[2] * xi[0]) / 2,(xi[0]*xi[0] * xi[2]*xi[2]*xi[2]) / 2,xi[2]*xi[2]*xi[2] * xi[0] * xi[1],
	(xi[0]*xi[0]*xi[0]*xi[0] * xi[1]) / 6,(xi[1]*xi[1]*xi[1]*xi[1] * xi[0]) / 6,(xi[2]*xi[2]*xi[2] * xi[0] * xi[1]) / 6,(xi[0]*xi[0]*xi[0] * xi[1]*xi[1]) / 2,(xi[0]*xi[0] * xi[1]*xi[1]*xi[1]) / 2,(xi[1]*xi[1]*xi[1] * xi[0] * xi[2]) / 2,(xi[1]*xi[1] * xi[0] * xi[2]*xi[2]) / 2,(xi[0]*xi[0] * xi[1] * xi[2]*xi[2]) / 2,(xi[0]*xi[0]*xi[0] * xi[1] * xi[2]) / 2,xi[0]*xi[0] * xi[1]*xi[1] * xi[2],
	(xi[0]*xi[0]*xi[0] * xi[1] * xi[2]) / 6,(xi[1]*xi[1]*xi[1]*xi[1] * xi[2]) / 6,(xi[2]*xi[2]*xi[2]*xi[2] * xi[1]) / 6,(xi[0]*xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1]*xi[1] * xi[0] * xi[2]) / 2,(xi[1]*xi[1]*xi[1] * xi[2]*xi[2]) / 2,(xi[1]*xi[1] * xi[2]*xi[2]*xi[2]) / 2,(xi[2]*xi[2]*xi[2] * xi[0] * xi[1]) / 2,(xi[0]*xi[0] * xi[1] * xi[2]*xi[2]) / 2,xi[1]*xi[1] * xi[0] * xi[2]*xi[2],
	(xi[0]*xi[0]*xi[0]*xi[0] * xi[2]) / 6,(xi[1]*xi[1]*xi[1] * xi[0] * xi[2]) / 6,(xi[2]*xi[2]*xi[2]*xi[2] * xi[0]) / 6,(xi[0]*xi[0]*xi[0] * xi[1] * xi[2]) / 2,(xi[0]*xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1] * xi[0] * xi[2]*xi[2]) / 2,(xi[2]*xi[2]*xi[2] * xi[0] * xi[1]) / 2,(xi[0]*xi[0] * xi[2]*xi[2]*xi[2]) / 2,(xi[0]*xi[0]*xi[0] * xi[2]*xi[2]) / 2,xi[0]*xi[0] * xi[1] * xi[2]*xi[2]
		};

		double A31[10][3] = {
			xi[0]*xi[0]*xi[0]*xi[0],xi[0]*xi[0]*xi[0] * xi[1],xi[0]*xi[0]*xi[0] * xi[2],
	xi[0] * xi[1]*xi[1]*xi[1],xi[1]*xi[1]*xi[1]*xi[1],xi[1]*xi[1]*xi[1] * xi[2],
	xi[0] * xi[2]*xi[2]*xi[2],xi[1] * xi[2]*xi[2]*xi[2],xi[2]*xi[2]*xi[2]*xi[2],
	xi[0]*xi[0]*xi[0] * xi[1],xi[0]*xi[0] * xi[1]*xi[1],xi[0]*xi[0] * xi[1] * xi[2],
	xi[0]*xi[0] * xi[1]*xi[1],xi[0] * xi[1]*xi[1]*xi[1],xi[0] * xi[1]*xi[1] * xi[2],
	xi[0] * xi[1]*xi[1] * xi[2],xi[1]*xi[1]*xi[1] * xi[2],xi[1]*xi[1] * xi[2]*xi[2],
	xi[0] * xi[1] * xi[2]*xi[2],xi[1]*xi[1] * xi[2]*xi[2],xi[1] * xi[2]*xi[2]*xi[2],
	xi[0]*xi[0] * xi[2]*xi[2],xi[0] * xi[1] * xi[2]*xi[2],xi[0] * xi[2]*xi[2]*xi[2],
	xi[0]*xi[0]*xi[0] * xi[2],xi[0]*xi[0] * xi[1] * xi[2],xi[0]*xi[0] * xi[2]*xi[2],
	xi[0]*xi[0] * xi[1] * xi[2],xi[0] * xi[1]*xi[1] * xi[2],xi[0] * xi[1] * xi[2]*xi[2]
		};
		
		double A32[10][6] = { (pow(xi[0], 5)) / 2,(xi[0]*xi[0]*xi[0] * xi[1]*xi[1]) / 2,(xi[0]*xi[0]*xi[0] * xi[2]*xi[2]) / 2,xi[0]*xi[0]*xi[0]*xi[0] * xi[1],xi[0]*xi[0]*xi[0] * xi[1] * xi[2],xi[0]*xi[0]*xi[0]*xi[0] * xi[2],
	(xi[0]*xi[0] * xi[1]*xi[1]*xi[1]) / 2,(pow(xi[1], 5)) / 2,(xi[1]*xi[1]*xi[1] * xi[2]*xi[2]) / 2,xi[1]*xi[1]*xi[1]*xi[1] * xi[0],xi[1]*xi[1]*xi[1]*xi[1] * xi[2],xi[1]*xi[1]*xi[1] * xi[0] * xi[2],
	(xi[0]*xi[0] * xi[2]*xi[2]*xi[2]) / 2,(xi[1]*xi[1] * xi[2]*xi[2]*xi[2]) / 2,(pow(xi[2], 5)) / 2,xi[2]*xi[2]*xi[2] * xi[0] * xi[1],xi[2]*xi[2]*xi[2]*xi[2] * xi[1],xi[2]*xi[2]*xi[2]*xi[2] * xi[0],
	(xi[0]*xi[0]*xi[0]*xi[0] * xi[1]) / 2,(xi[0]*xi[0] * xi[1]*xi[1]*xi[1]) / 2,(xi[0]*xi[0] * xi[1] * xi[2]*xi[2]) / 2,xi[0]*xi[0]*xi[0] * xi[1]*xi[1],xi[0]*xi[0] * xi[1]*xi[1] * xi[2],xi[0]*xi[0]*xi[0] * xi[1] * xi[2],
	(xi[0]*xi[0]*xi[0] * xi[1]*xi[1]) / 2,(xi[1]*xi[1]*xi[1]*xi[1] * xi[0]) / 2,(xi[1]*xi[1] * xi[0] * xi[2]*xi[2]) / 2,xi[0]*xi[0] * xi[1]*xi[1]*xi[1],xi[1]*xi[1]*xi[1] * xi[0] * xi[2],xi[0]*xi[0] * xi[1]*xi[1] * xi[2],
	(xi[0]*xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1]*xi[1]*xi[1] * xi[2]) / 2,(xi[1]*xi[1] * xi[2]*xi[2]*xi[2]) / 2,xi[1]*xi[1]*xi[1] * xi[0] * xi[2],xi[1]*xi[1]*xi[1] * xi[2]*xi[2],xi[1]*xi[1] * xi[0] * xi[2]*xi[2],
	(xi[0]*xi[0] * xi[1] * xi[2]*xi[2]) / 2,(xi[1]*xi[1]*xi[1] * xi[2]*xi[2]) / 2,(xi[2]*xi[2]*xi[2]*xi[2] * xi[1]) / 2,xi[1]*xi[1] * xi[0] * xi[2]*xi[2],xi[1]*xi[1] * xi[2]*xi[2]*xi[2],xi[2]*xi[2]*xi[2] * xi[0] * xi[1],
	(xi[0]*xi[0]*xi[0] * xi[2]*xi[2]) / 2,(xi[1]*xi[1] * xi[0] * xi[2]*xi[2]) / 2,(xi[2]*xi[2]*xi[2]*xi[2] * xi[0]) / 2,xi[0]*xi[0] * xi[1] * xi[2]*xi[2],xi[2]*xi[2]*xi[2] * xi[0] * xi[1],xi[0]*xi[0] * xi[2]*xi[2]*xi[2],
	(xi[0]*xi[0]*xi[0]*xi[0] * xi[2]) / 2,(xi[0]*xi[0] * xi[1]*xi[1] * xi[2]) / 2,(xi[0]*xi[0] * xi[2]*xi[2]*xi[2]) / 2,xi[0]*xi[0]*xi[0] * xi[1] * xi[2],xi[0]*xi[0] * xi[1] * xi[2]*xi[2],xi[0]*xi[0]*xi[0] * xi[2]*xi[2],
	(xi[0]*xi[0]*xi[0] * xi[1] * xi[2]) / 2,(xi[1]*xi[1]*xi[1] * xi[0] * xi[2]) / 2,(xi[2]*xi[2]*xi[2] * xi[0] * xi[1]) / 2,xi[0]*xi[0] * xi[1]*xi[1] * xi[2],xi[1]*xi[1] * xi[0] * xi[2]*xi[2],xi[0]*xi[0] * xi[1] * xi[2]*xi[2] 
		};

		double A33[10][10] = {
			(pow(xi[0],6)) / 6,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[1]) / 6,(xi[0] * xi[0] * xi[0] * xi[2] * xi[2] * xi[2]) / 6,(pow(xi[0],5) * xi[1]) / 2,(xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[1]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[0] * xi[2] * xi[2]) / 2,(pow(xi[0],5) * xi[2]) / 2,xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[2],
	(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[1]) / 6,(pow(xi[1],6)) / 6,(xi[1] * xi[1] * xi[1] * xi[2] * xi[2] * xi[2]) / 6,(xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[0]) / 2,(pow(xi[1],5) * xi[0]) / 2,(pow(xi[1],5) * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[2] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[0] * xi[2]) / 2,xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[2],
	(xi[0] * xi[0] * xi[0] * xi[2] * xi[2] * xi[2]) / 6,(xi[1] * xi[1] * xi[1] * xi[2] * xi[2] * xi[2]) / 6,(pow(xi[2],6)) / 6,(xi[2] * xi[2] * xi[2] * xi[0] * xi[0] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[1] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[2] * xi[1] * xi[1]) / 2,(pow(xi[2],5) * xi[1]) / 2,(pow(xi[2],5) * xi[0]) / 2,(xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[0]) / 2,xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[1],
	(pow(xi[0],5) * xi[1]) / 6,(xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[0]) / 6,(xi[2] * xi[2] * xi[2] * xi[0] * xi[0] * xi[1]) / 6,(xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[1]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[1]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[0] * xi[2]) / 2,(xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[2]) / 2,xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[2],
	(xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[1]) / 6,(pow(xi[1],5) * xi[0]) / 6,(xi[2] * xi[2] * xi[2] * xi[0] * xi[1] * xi[1]) / 6,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[1]) / 2,(xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[0]) / 2,(xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[2]) / 2,xi[1] * xi[1] * xi[1] * xi[0] * xi[0] * xi[2],
	(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[2]) / 6,(pow(xi[1],5) * xi[2]) / 6,(xi[2] * xi[2] * xi[2] * xi[2] * xi[1] * xi[1]) / 6,(xi[1] * xi[1] * xi[1] * xi[0] * xi[0] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[2] * xi[2] * xi[2]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[1] * xi[1]) / 2,(xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,xi[1] * xi[1] * xi[1] * xi[0] * xi[2] * xi[2],
	(xi[0] * xi[0] * xi[0] * xi[1] * xi[2] * xi[2]) / 6,(xi[1] * xi[1] * xi[1] * xi[1] * xi[2] * xi[2]) / 6,(pow(xi[2],5) * xi[1]) / 6,(xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[2] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[2] * xi[2] * xi[2]) / 2,(xi[2] * xi[2] * xi[2] * xi[2] * xi[1] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[0] * xi[1]) / 2,xi[2] * xi[2] * xi[2] * xi[0] * xi[1] * xi[1],
	(xi[0] * xi[0] * xi[0] * xi[0] * xi[2] * xi[2]) / 6,(xi[1] * xi[1] * xi[1] * xi[0] * xi[2] * xi[2]) / 6,(pow(xi[2],5) * xi[0]) / 6,(xi[0] * xi[0] * xi[0] * xi[1] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[1] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[0]) / 2,(xi[0] * xi[0] * xi[0] * xi[2] * xi[2] * xi[2]) / 2,xi[2] * xi[2] * xi[2] * xi[0] * xi[0] * xi[1],
	(pow(xi[0],5) * xi[2]) / 6,(xi[1] * xi[1] * xi[1] * xi[0] * xi[0] * xi[2]) / 6,(xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[0]) / 6,(xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[2]) / 2,(xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[0] * xi[1]) / 2,(xi[0] * xi[0] * xi[0] * xi[2] * xi[2] * xi[2]) / 2,(xi[0] * xi[0] * xi[0] * xi[0] * xi[2] * xi[2]) / 2,xi[0] * xi[0] * xi[0] * xi[1] * xi[2] * xi[2],
	(xi[0] * xi[0] * xi[0] * xi[0] * xi[1] * xi[2]) / 6,(xi[1] * xi[1] * xi[1] * xi[1] * xi[0] * xi[2]) / 6,(xi[2] * xi[2] * xi[2] * xi[2] * xi[0] * xi[1]) / 6,(xi[0] * xi[0] * xi[0] * xi[1] * xi[1] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[0] * xi[2]) / 2,(xi[1] * xi[1] * xi[1] * xi[0] * xi[2] * xi[2]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[1] * xi[1]) / 2,(xi[2] * xi[2] * xi[2] * xi[0] * xi[0] * xi[1]) / 2,(xi[0] * xi[0] * xi[0] * xi[1] * xi[2] * xi[2]) / 2,xi[0] * xi[0] * xi[1] * xi[1] * xi[2] * xi[2]
		};

		for (int row = 0; row < 9; row++)
		{
			for (int col = 0; col < 9; col++)
			{
				A->addCoeff(row, col, omega * dv * tempA1[row][col]);
			}
			for (int col = 9; col < 19; col++)
			{
				A->addCoeff(row, col, omega * dv * tempA2[row][col-9]);
			}
		}
		for (int row = 9; row < 19; row++)
		{
			for (int col = 0; col < 3; col++)
			{
				A->addCoeff(row, col, omega * dv * A31[row-9][col]);
			}
			for (int col = 3; col < 9; col++)
			{
				A->addCoeff(row, col, omega * dv * A32[row-9][col-3]);
			}
			for (int col = 9; col < 19; col++)
			{
				A->addCoeff(row, col, omega * dv * A33[row-9][col-9]);
			}
		}
	}
}

void pdsolve::vec_gd2D(double g[],double d[], Matrix *A, pdFamily* p_fami,double xi[], datModel & o_dat)
{
	// A(p,q) size 5*5;
	Vector *vec_xi, *vecg;
	vec_xi = new Vector(5);
	vecg = new Vector(5);
	
	vec_xi->setCoeff(0, xi[0]);
	vec_xi->setCoeff(1, xi[1]);
	vec_xi->setCoeff(2, xi[0] * xi[0]);
	vec_xi->setCoeff(3, xi[1] * xi[1]);
	vec_xi->setCoeff(4, xi[0] * xi[1]);

	matoperat.PLUSolve(A, vec_xi, vecg);

	g[0] = vecg->d_getCoeff(0);
	g[1] = vecg->d_getCoeff(1);
	d[0] = vecg->d_getCoeff(2);
	d[1] = vecg->d_getCoeff(3);
	d[2] = vecg->d_getCoeff(4);
	delete vec_xi;
	delete vecg;
	vec_xi = NULL;
	vecg = NULL;

}

void pdsolve::vec_gd3D(double g[], double d[], Matrix* A, pdFamily* p_fami, double xi[], datModel& o_dat)
{
	// A(p,q) size 9*9;
	Vector* vec_xi, * vecg;
	vec_xi = new Vector(9);
	vecg = new Vector(9);
	//=====
	vec_xi->setCoeff(0, xi[0]);
	vec_xi->setCoeff(1, xi[1]);
	vec_xi->setCoeff(2, xi[2]);
	vec_xi->setCoeff(3, xi[0] * xi[0]);
	vec_xi->setCoeff(4, xi[1] * xi[1]);
	vec_xi->setCoeff(5, xi[2] * xi[2]);
	vec_xi->setCoeff(6, xi[0] * xi[1]);
	vec_xi->setCoeff(7, xi[1] * xi[2]);
	vec_xi->setCoeff(8, xi[2] * xi[0]);
	//===solve==
	matoperat.PLUSolve(A, vec_xi, vecg);
	//==assign==
	g[0] = vecg->d_getCoeff(0);
	g[1] = vecg->d_getCoeff(1);
	g[2] = vecg->d_getCoeff(2);
	for (int i = 0; i < 6; i++)
	{
		d[i] = vecg->d_getCoeff(i + 3);
	}
	delete vec_xi;
	delete vecg;
	vec_xi = NULL;
	vecg = NULL;
}

void pdsolve::vec_gd3D3rd(double g[], double d[], Matrix* A, pdFamily* p_fami, double xi[], datModel& o_dat)
{
	// A(p,q) size 19*19;
	Vector* vec_xi, * vecg;
	vec_xi = new Vector(19);
	vecg = new Vector(19);
	//=====
	vec_xi->setCoeff(0, xi[0]);
	vec_xi->setCoeff(1, xi[1]);
	vec_xi->setCoeff(2, xi[2]);
	vec_xi->setCoeff(3, xi[0] * xi[0]);
	vec_xi->setCoeff(4, xi[1] * xi[1]);
	vec_xi->setCoeff(5, xi[2] * xi[2]);
	vec_xi->setCoeff(6, xi[0] * xi[1]);
	vec_xi->setCoeff(7, xi[1] * xi[2]);
	vec_xi->setCoeff(8, xi[2] * xi[0]);
	//==r3
	vec_xi->setCoeff(9, xi[0] * xi[0] * xi[0]);
	vec_xi->setCoeff(10, xi[1] * xi[1] * xi[1]);
	vec_xi->setCoeff(11, xi[2] * xi[2] * xi[2]);
	vec_xi->setCoeff(12, xi[0] * xi[0] * xi[1]);
	vec_xi->setCoeff(13, xi[0] * xi[1] * xi[1]);
	vec_xi->setCoeff(14, xi[1] * xi[1] * xi[2]);
	vec_xi->setCoeff(15, xi[1] * xi[2] * xi[2]);
	vec_xi->setCoeff(16, xi[0] * xi[2] * xi[2]);
	vec_xi->setCoeff(17, xi[0] * xi[0] * xi[2]);
	vec_xi->setCoeff(18, xi[0] * xi[1] * xi[2]);
	//===solve==
	matoperat.PLUSolve(A, vec_xi, vecg);
	//==assign==
	g[0] = vecg->d_getCoeff(0);
	g[1] = vecg->d_getCoeff(1);
	g[2] = vecg->d_getCoeff(2);
	for (int i = 0; i < 6; i++)
	{
		d[i] = vecg->d_getCoeff(i + 3);
	}
	delete vec_xi;
	delete vecg;
	vec_xi = NULL;
	vecg = NULL;
}

void pdsolve::matG2D(Matrix * G, Matrix * A, pdFamily* p_fami, int m, datModel & o_dat)
{
	double xk[3], xm[3], xi[3];
	int Nid_k, Nid_m;
	Nid_k = p_fami->getNodeID(0);
	Nid_m = p_fami->getNodeID(m);
	o_dat.op_getNode(Nid_k - 1)->getcoor(xk);
	o_dat.op_getNode(Nid_m - 1)->getcoor(xm);
	for (int i = 0; i < 3; i++)
	{
		xi[i] = xm[i] - xk[i];
	}

	int mu_km = p_fami->getbondstate(m);

	double g[2], d[3], omega, temp;
	vec_gd2D(g, d, A, p_fami, xi, o_dat);
	omega = inflFunc(xi, p_fami, o_dat);

	temp = mu_km* omega*((2 * cd_mu + cd_lambda)*d[0] + cd_mu * d[1]);
	G->setCoeff(0, 0, temp);
	temp = mu_km* omega*(cd_lambda + cd_mu)*d[2];
	G->setCoeff(0, 1, temp);
	G->setCoeff(1, 0, temp);
	temp = mu_km* omega * ((2 * cd_mu + cd_lambda)*d[1] + cd_mu * d[0]);
	G->setCoeff(1, 1, temp);
}

void pdsolve::matG3D(Matrix* G, Matrix* A, pdFamily* p_fami, int m, datModel& o_dat)
{
	double xk[3], xm[3], xi[3];
	int Nid_k, Nid_m;
	Nid_k = p_fami->getNodeID(0);
	Nid_m = p_fami->getNodeID(m);
	o_dat.op_getNode(Nid_k - 1)->getcoor(xk);
	o_dat.op_getNode(Nid_m - 1)->getcoor(xm);
	for (int i = 0; i < 3; i++)
	{
		xi[i] = xm[i] - xk[i];
	}

	int mu_km = p_fami->getbondstate(m);

	double g[3], d[6], omega, temp, trD;
	if (o_dat.ci_TESflag==2)
	{
		vec_gd3D(g, d, A, p_fami, xi, o_dat);
	}
	else if (o_dat.ci_TESflag==3)
	{
		vec_gd3D3rd(g, d, A, p_fami, xi, o_dat);
	}
	omega = inflFunc(xi, p_fami, o_dat);
	trD = d[0] + d[1] + d[2];
	//==assin values
	temp = mu_km * omega * ((cd_mu + cd_lambda) * d[0] + cd_mu * trD);
	G->setCoeff(0, 0, temp);
	temp = mu_km * omega * (cd_lambda + cd_mu) * d[3];
	G->setCoeff(0, 1, temp);
	G->setCoeff(1, 0, temp);
	temp = mu_km * omega * (cd_lambda + cd_mu) * d[5];
	G->setCoeff(0, 2, temp);
	G->setCoeff(2, 0, temp);
	temp = mu_km * omega * ((cd_mu + cd_lambda) * d[1] + cd_mu *trD);
	G->setCoeff(1, 1, temp);
	temp = mu_km * omega * (cd_lambda + cd_mu) * d[4];
	G->setCoeff(1, 2, temp);
	G->setCoeff(2, 1, temp);
	temp = mu_km * omega * ((cd_mu + cd_lambda) * d[2] + cd_mu * trD);
	G->setCoeff(2, 2, temp);
}

void pdsolve::matH2D(Matrix * H, pdFamily* p_fami, datModel & o_dat)
{
	H->zero();
	Matrix *A,*G;
	A = new Matrix(5, 5);
	shapTens2D(A, p_fami, o_dat);
	G = new Matrix(2, 2);
	double dv_m, temp;
	int  NID_m;
	int numNodeOfFam = p_fami->getNumNode();
	for (int m = 1; m < numNodeOfFam; m++)
	{
		NID_m = p_fami->getNodeID(m);
		dv_m = o_dat.op_getNode(NID_m - 1)->getvolume();
		matG2D(G, A, p_fami, m, o_dat);
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				temp = G->d_getCoeff(i, j)*dv_m;
				H->addCoeff(i, j, -temp);
				H->setCoeff(i, j + 2 * m, temp);
			}
		}

	}

	delete A, delete G;
	A = NULL; G = NULL;
}

void pdsolve::matH3D(Matrix* H, pdFamily* p_fami, datModel& o_dat)
{
	H->zero();
	Matrix* A, * G;
	if (o_dat.ci_TESflag==2)
	{
		A = new Matrix(9, 9);
		shapTens3D(A, p_fami, o_dat);
	}
	else if (o_dat.ci_TESflag == 3)
	{
		A = new Matrix(19, 19);
		shapTens3D3rd(A, p_fami, o_dat);
	}
	G = new Matrix(3, 3);
	double dv_m, temp;
	int  NID_m;
	int numNodeOfFam = p_fami->getNumNode();
	for (int m = 1; m < numNodeOfFam; m++)
	{
		NID_m = p_fami->getNodeID(m);
		dv_m = o_dat.op_getNode(NID_m - 1)->getvolume();
		matG3D(G, A, p_fami, m, o_dat);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				temp = G->d_getCoeff(i, j) * dv_m;
				H->addCoeff(i, j, -temp);
				H->setCoeff(i, j + 3 * m, temp);
			}
		}
	}
	delete A, delete G;
	A = NULL; G = NULL;
	
}

void pdsolve::assembleInterWorkPD(datModel & o_dat)
{
	// assemble H============= ;
	int numDime = o_dat.ci_Numdimen;
	pdFamily* temP_fami;
	if (numDime==2)
	{
		//==================2D=================
		Matrix* H;
		int numNodeOfFami, NID_k, NID_m, eqindex_row, eqindex_col;
		double dv_k, temp, tempu;
		
		for (int famkk = 0; famkk < o_dat.getTotnumFami(); famkk++)
		{
			temP_fami = o_dat.op_getFami(famkk);
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			dv_k = o_dat.op_getNode(NID_k - 1)->getvolume();
			//==get H==
			H = new Matrix(numDime, numDime * numNodeOfFami);
			matH2D(H, temP_fami, o_dat);
			//====assemble===
			for (int i = 0; i < numDime; i++)
			{
				eqindex_row = o_dat.op_getNode(NID_k - 1)->op_getDof(i)->i_getEqInde();
				if (eqindex_row != -1)
				{
					for (int m = 0; m < numNodeOfFami; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < numDime; j++)
						{
							eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
							temp = -(H->d_getCoeff(i, numDime * m + j) * dv_k);
							if (eqindex_col != -1)
							{
								cop_Ku->addCoeff(eqindex_row, eqindex_col, temp);
							}
							else
							{
								//fixIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getFixIndex();
								//cop_Kp->addCoeff(eqindex_row, fixIndex_col, -temp);
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cop_F->addCoeff(eqindex_row, -temp * tempu);
							}
						}

					}
				}

			}
			delete H;
			H = NULL;
		}
	}
	else if (numDime==3)
	{
		//==================3D=================
		Matrix* H;
		int numNodeOfFami, NID_k, NID_m, eqindex_row, eqindex_col;
		double dv_k, temp, tempu;
		for (int famkk = 0; famkk < o_dat.getTotnumFami(); famkk++)
		{
			temP_fami = o_dat.op_getFami(famkk);
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			dv_k = o_dat.op_getNode(NID_k - 1)->getvolume();
			//==get H==
			H = new Matrix(numDime, numDime * numNodeOfFami);
			matH3D(H, temP_fami, o_dat);
			//====assemble===
			for (int i = 0; i < numDime; i++)
			{
				eqindex_row = o_dat.op_getNode(NID_k - 1)->op_getDof(i)->i_getEqInde();
				if (eqindex_row != -1)
				{
					for (int m = 0; m < numNodeOfFami; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < numDime; j++)
						{
							eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
							temp = -(H->d_getCoeff(i, numDime * m + j) * dv_k);
							if (eqindex_col != -1)
							{
								cop_Ku->addCoeff(eqindex_row, eqindex_col, temp);
							}
							else
							{
								//fixIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getFixIndex();
								//cop_Kp->addCoeff(eqindex_row, fixIndex_col, -temp);
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cop_F->addCoeff(eqindex_row, -temp * tempu);
							}
						}

					}
				}

			}
			delete H;
			H = NULL;
		}
	}
}

void pdsolve::assembleInterWorkPD_CSRformat(datModel& o_dat)
{
	// assemble H============= ;
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	pdFamily* temP_fami;
	if (numDime == 2)
	{
		//==================2D=================
		Matrix* H;
		int numNodeOfFami, NID_k, NID_m, eqindex_row, eqindex_col;
		long long int i_temp;
		double dv_k, temp, tempu;
		for (int famkk = startP; famkk < endP; famkk++)
		{
			temP_fami = o_dat.op_getFami(famkk);
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			dv_k = o_dat.op_getNode(NID_k - 1)->getvolume();
			//==get H==
			H = new Matrix(numDime, numDime * numNodeOfFami);
			matH2D(H, temP_fami, o_dat);
			//====assemble===
			for (int i = 0; i < numDime; i++)
			{
				eqindex_row = o_dat.op_getNode(NID_k - 1)->op_getDof(i)->i_getEqInde();
				if (eqindex_row != -1)
				{
					for (int m = 0; m < numNodeOfFami; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < numDime; j++)
						{
							
							temp = -(H->d_getCoeff(i, numDime * m + j) * dv_k);
							if (cb_InteralForce == false)
							{
								// non-dynamic solver;
								eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								if (eqindex_col != -1)
								{
									i_temp = findCSRIndexOfMat(eqindex_row, eqindex_col);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								else
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								}
							}
							else
							{
								// dynamic solver;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
							}
							
						}

					}
				}

			}
			delete H;
			H = NULL;
		}
	}
	else if (numDime == 3)
	{
		//==================3D=================
		Matrix* H;
		int numNodeOfFami, NID_k, NID_m, eqindex_row, eqindex_col;
		long long int i_temp;
		double dv_k, temp, tempu;
		for (int famkk = startP; famkk < endP; famkk++)
		{
			temP_fami = o_dat.op_getFami(famkk);
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			dv_k = o_dat.op_getNode(NID_k - 1)->getvolume();
			//==get H==
			H = new Matrix(numDime, numDime * numNodeOfFami);
			matH3D(H, temP_fami, o_dat);
			//====assemble===
			for (int i = 0; i < numDime; i++)
			{
				eqindex_row = o_dat.op_getNode(NID_k - 1)->op_getDof(i)->i_getEqInde();
				if (eqindex_row != -1)
				{
					for (int m = 0; m < numNodeOfFami; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < numDime; j++)
						{
							temp = -(H->d_getCoeff(i, numDime * m + j) * dv_k);
							if (cb_InteralForce == false)
							{
								// non-dynamic solver;
								eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								if (eqindex_col != -1)
								{
									i_temp = findCSRIndexOfMat(eqindex_row, eqindex_col);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								else
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								}
							}
							else
							{
								// dynamic solver;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
							}
						}
					}
				}

			}
			delete H;
			H = NULL;
		}
	}
}

void pdsolve::matC2D(Matrix * C, pdFamily* p_fami, datModel & o_dat)
{
	C->zero();
	Matrix *A;
	A = new Matrix(5, 5);
	shapTens2D(A, p_fami, o_dat);
	double g[2], d[3];
	double xk[3], xm[3], xi[3], dv_m, temp,omega;
	int NID_k, NID_m, mu_km;
	//1st node;
	NID_k = p_fami->getNodeID(0);
	o_dat.op_getNode(NID_k - 1)->getcoor(xk);
	// integration
	int numNodeOfFam = p_fami->getNumNode();
	for (int m = 1; m < numNodeOfFam; m++)
	{
		//m-th node
		NID_m = p_fami->getNodeID(m);
		o_dat.op_getNode(NID_m - 1)->getcoor(xm);
		dv_m = o_dat.op_getNode(NID_m - 1)->getvolume();
		for (int i = 0; i < 3; i++)
		{
			xi[i] = xm[i] - xk[i];
		}
		mu_km = p_fami->getbondstate(m);
		omega = inflFunc(xi, p_fami, o_dat);
		vec_gd2D(g, d, A, p_fami, xi, o_dat);
		temp = mu_km* omega * g[0] * dv_m;
		C->setCoeff(0, 2 * m, temp);
		C->setCoeff(2, 2 * m + 1, temp);
		C->addCoeff(0, 0, -temp);
		C->addCoeff(2, 1, -temp);
		temp = mu_km* omega * g[1] * dv_m;
		C->setCoeff(1, 2 * m + 1, temp);
		C->setCoeff(2, 2 * m, temp);
		C->addCoeff(1, 1, -temp);
		C->addCoeff(2, 0, -temp);
	}
	delete A;
	A = NULL;
}

void pdsolve::matC3D(Matrix* C, pdFamily* p_fami, datModel& o_dat)
{
	C->zero();
	Matrix* A;
	if (o_dat.ci_TESflag == 2)
	{
		A = new Matrix(9, 9);
		shapTens3D(A, p_fami, o_dat);
	}
	else if (o_dat.ci_TESflag == 3)
	{
		A = new Matrix(19, 19);
		shapTens3D3rd(A, p_fami, o_dat);
	}
	else
	{
		printf("ERROR: TES flag is wrong\n");
		exit(0);
	}
	double g[3], d[6];
	double xk[3], xm[3], xi[3], dv_m, omega;
	int NID_k, NID_m, mu_km;
	//1st node;
	NID_k = p_fami->getNodeID(0);
	o_dat.op_getNode(NID_k - 1)->getcoor(xk);
	// integration
	int numNodeOfFam = p_fami->getNumNode();
	double d_c[3];
	for (int m = 1; m < numNodeOfFam; m++)
	{
		//m-th node
		NID_m = p_fami->getNodeID(m);
		o_dat.op_getNode(NID_m - 1)->getcoor(xm);
		dv_m = o_dat.op_getNode(NID_m - 1)->getvolume();
		for (int i = 0; i < 3; i++)
		{
			xi[i] = xm[i] - xk[i];
		}
		mu_km = p_fami->getbondstate(m);
		omega = inflFunc(xi, p_fami, o_dat);
		if (o_dat.ci_TESflag == 2)
		{
			vec_gd3D(g, d, A, p_fami, xi, o_dat);
		}
		else if (o_dat.ci_TESflag == 3)
		{
			vec_gd3D3rd(g, d, A, p_fami, xi, o_dat);
		}
		for (int i = 0; i < 3; i++)
		{
			d_c[i]= mu_km * omega * g[i] * dv_m;
		}
		C->setCoeff(0, 3 * m, d_c[0]);
		C->setCoeff(1, 3 * m + 1, d_c[1]);
		C->setCoeff(2, 3 * m + 2, d_c[2]);
		C->setCoeff(3, 3 * m, d_c[1]);
		C->setCoeff(3, 3 * m + 1, d_c[0]);
		C->setCoeff(4, 3 * m + 1, d_c[2]);
		C->setCoeff(4, 3 * m + 2, d_c[1]);
		C->setCoeff(5, 3 * m, d_c[2]);
		C->setCoeff(5, 3 * m + 2, d_c[0]);
		//===first 3 column;
		C->addCoeff(0, 0, -d_c[0]);
		C->addCoeff(1, 1, -d_c[1]);
		C->addCoeff(2, 2, -d_c[2]);
		C->addCoeff(3, 0, -d_c[1]);
		C->addCoeff(3, 1, -d_c[0]);
		C->addCoeff(4, 1, -d_c[2]);
		C->addCoeff(4, 2, -d_c[1]);
		C->addCoeff(5, 0, -d_c[2]);
		C->addCoeff(5, 2, -d_c[0]);
	}
	delete A;
	A = NULL;
}



void pdsolve::assemblePDBEwork(datModel & o_dat)
{
	int numDime = o_dat.ci_Numdimen;
	pdFamily* temP_fami;
	if (numDime==2)
	{
		//===========2D=======================================
		Matrix* C, * N, * DC, * NDC;
		N = new Matrix(2, 3);
		N->zero();
		int bLinNID[2], NID_m;
		double xL1[2], xL2[2], nx, ny, temp, tempu;
		bool bL1_isPDNODE, bL2_isPDNODE;
		int numNodeFam, famkk;
		int eqIndex_row[2], eqIndex_col[2];
		for (int bL = 0; bL < o_dat.getTotnumPDBEs(); bL++)
		{
			o_dat.op_getPDBE(bL)->getNodeID(bLinNID);
			bL1_isPDNODE = o_dat.op_getNode(bLinNID[0] - 1)->getNodeType();
			bL2_isPDNODE = o_dat.op_getNode(bLinNID[1] - 1)->getNodeType();
			if (bL1_isPDNODE == false || bL2_isPDNODE == false)
			{
				cout << "Warning: This node of PD domain boundary is not PD node!" << endl;
				exit(0);
			}
			o_dat.op_getNode(bLinNID[0] - 1)->getcoor(xL1);
			o_dat.op_getNode(bLinNID[1] - 1)->getcoor(xL2);
			nx = xL2[1] - xL1[1];
			ny = xL1[0] - xL2[0];
			N->setCoeff(0, 0, nx);
			N->setCoeff(0, 2, ny);
			N->setCoeff(1, 1, ny);
			N->setCoeff(1, 2, nx);
			//=============t1=================
			famkk = o_dat.op_getNode(bLinNID[0] - 1)->getFamID() - 1;
			temP_fami = o_dat.op_getFami(famkk);
			numNodeFam = temP_fami->getNumNode();
			C = new Matrix(3, 2 * numNodeFam);
			DC = new Matrix(3, 2 * numNodeFam);
			NDC = new Matrix(2, 2 * numNodeFam);
			matC2D(C, temP_fami, o_dat);
			matoperat.matMultiply(cop_D, C, DC);
			matoperat.matMultiply(N, DC, NDC);
			//u1;
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[0] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
							temp = 1.0 / 3 * NDC->d_getCoeff(i, 2 * m + j);
							if (eqIndex_col[j] != -1)
							{
								cop_Ku->addCoeff(eqIndex_row[i], eqIndex_col[j], temp);
							}
							else
							{
								//fixIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getFixIndex();
								//cop_Kp->addCoeff(eqIndex_row[i], fixIndex_col[j], -temp);
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cop_F->addCoeff(eqIndex_row[i], -temp * tempu);
							}
						}
					}
				}

			}
			//u2; node 2
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[1] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
							temp = 1.0 / 6 * NDC->d_getCoeff(i, 2 * m + j);
							if (eqIndex_col[j] != -1)
							{
								cop_Ku->addCoeff(eqIndex_row[i], eqIndex_col[j], temp);
							}
							else
							{
								//fixIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getFixIndex();
								//cop_Kp->addCoeff(eqIndex_row[i], fixIndex_col[j], -temp);
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cop_F->addCoeff(eqIndex_row[i], -temp * tempu);

							}
						}
					}
				}

			}
			delete C, delete DC, delete NDC;
			C = NULL; DC = NULL; NDC = NULL;
			//===========t2**** node 2========
			famkk = o_dat.op_getNode(bLinNID[1] - 1)->getFamID() - 1;
			temP_fami = o_dat.op_getFami(famkk);
			numNodeFam = temP_fami->getNumNode();
			C = new Matrix(3, 2 * numNodeFam);
			DC = new Matrix(3, 2 * numNodeFam);
			NDC = new Matrix(2, 2 * numNodeFam);
			matC2D(C, temP_fami, o_dat);
			matoperat.matMultiply(cop_D, C, DC);
			matoperat.matMultiply(N, DC, NDC);
			//u1; node 1
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[0] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
							temp = 1.0 / 6 * NDC->d_getCoeff(i, 2 * m + j);
							if (eqIndex_col[j] != -1)
							{
								cop_Ku->addCoeff(eqIndex_row[i], eqIndex_col[j], temp);
							}
							else
							{
								//fixIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getFixIndex();
								//cop_Kp->addCoeff(eqIndex_row[i], fixIndex_col[j], -temp);
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cop_F->addCoeff(eqIndex_row[i], -temp * tempu);
							}
						}
					}
				}

			}
			//u2; node 2===========
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[1] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
							temp = 1.0 / 3 * NDC->d_getCoeff(i, 2 * m + j);
							if (eqIndex_col[j] != -1)
							{
								cop_Ku->addCoeff(eqIndex_row[i], eqIndex_col[j], temp);
							}
							else
							{
								//fixIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getFixIndex();
								//cop_Kp->addCoeff(eqIndex_row[i], fixIndex_col[j], -temp);
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cop_F->addCoeff(eqIndex_row[i], -temp * tempu);
							}
						}
					}
				}

			}
			delete C, delete DC, delete NDC;
			C = NULL; DC = NULL; NDC = NULL;
		}
		delete N;
		N = NULL;
	}
	else if (numDime==3)
	{	
		//====================3D===================
		int nG, conNID[4];
		int famiID, numNodeFam, eqIndex_row, eqIndex_col;
		int NID_m;
		double p, q, wp, wq, fac, xN[4][3], temp, tempu;
		double N[4];//N[4] are shape functions;
		Matrix* mNt, * Nmat, * C, * f_mNt, * f_mNt_Nm, * f_mNt_Nm_D, * finMat;
		mNt = new Matrix(12, 3);
		Nmat = new Matrix(3, 6);
		f_mNt = new Matrix(12, 3);
		f_mNt_Nm = new Matrix(12, 6);
		f_mNt_Nm_D = new Matrix(12, 6);
		nG = o_globGP.i_getNumPts();
		//===========3D=======================================
		for (int pdbe = 0; pdbe <o_dat.getTotnumPDBEs(); pdbe++)
		{
			//==get PDBEs==
			o_dat.op_getPDBE(pdbe)->getNodeID(conNID);
			for (int i = 0; i < 4; i++)
			{
				o_dat.op_getNode(conNID[i] - 1)->getcoor(xN[i]);
			}
			//===Gauss integration===
			for (int mp = 0; mp < nG; mp++)
			{
				p = o_globGP.d_getGaussPt(mp);
				wp = o_globGP.d_getWeight(mp);
				for (int mq = 0; mq < nG; mq++)
				{
					q = o_globGP.d_getGaussPt(mq);
					wq = o_globGP.d_getWeight(mq);
					shapeFunctionQuad4N(N, p, q);
					matMathcalNt(mNt, p, q);
					for (int ei = 0; ei < 4; ei++)
					{
						//====get the final matrix for assembling;
						fac = N[ei] * wp * wq;
						matoperat.matMultiply(mNt, fac, f_mNt);
						matN_trans(Nmat, xN, p, q);
						famiID = o_dat.op_getNode(conNID[ei] - 1)->getFamID();
						temP_fami = o_dat.op_getFami(famiID - 1);
						numNodeFam= temP_fami->getNumNode();
						C = new Matrix(6, 3 * numNodeFam);
						finMat = new Matrix(12, 3 * numNodeFam);
						matC3D(C, temP_fami, o_dat);
						matoperat.matMultiply(f_mNt, Nmat, f_mNt_Nm);
						matoperat.matMultiply(f_mNt_Nm, cop_D, f_mNt_Nm_D);
						matoperat.matMultiply(f_mNt_Nm_D, C, finMat);
						//=======assembling======
						for (int i = 0; i < 4; i++)
						{
							for (int j = 0; j < 3; j++)
							{
								eqIndex_row = o_dat.op_getNode(conNID[i] - 1)->op_getDof(j)->i_getEqInde();
								if (eqIndex_row!=-1)
								{
									for (int m = 0; m < numNodeFam; m++)
									{
										NID_m = temP_fami->getNodeID(m);
										for (int jj = 0; jj < 3; jj++)
										{
											eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
											temp = finMat->d_getCoeff(i * 3 + j, 3 * m + jj);
											if (eqIndex_col!=-1)
											{
												cop_Ku->addCoeff(eqIndex_row, eqIndex_col, temp);
											}
											else
											{
												//fixIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getFixIndex();
												//cop_Kp->addCoeff(eqIndex_row, fixIndex_col, -temp);
												tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
												cop_F->addCoeff(eqIndex_row, -temp * tempu);
											}
										}
									}
								}
							}
						}
						//===finish assembling for this ei-th node;
						// delete;
						delete C, delete finMat;
						C = NULL, finMat = NULL;
					}
				}
			}
		}
		delete mNt, delete Nmat, delete f_mNt, delete f_mNt_Nm, delete f_mNt_Nm_D;
		mNt = NULL, Nmat = NULL, f_mNt = NULL, f_mNt_Nm = NULL, f_mNt_Nm_D = NULL;
	}
}

void pdsolve::assemblePDBEwork_CSRformat(datModel& o_dat)
{

	int numPDBEs, startP, endP;
	numPDBEs = o_dat.getTotnumPDBEs();
	startP = ci_rank * numPDBEs / ci_numProce;
	endP = (ci_rank + 1) * numPDBEs / ci_numProce;
	int numDime = o_dat.ci_Numdimen;
	pdFamily* temP_fami;
	if (numDime == 2)
	{
		//===========2D=======================================
		Matrix* C, * N, * DC, * NDC;
		N = new Matrix(2, 3);
		N->zero();
		int bLinNID[2], NID_m;
		long long int i_temp;
		double xL1[2], xL2[2], nx, ny, temp, tempu;
		bool bL1_isPDNODE, bL2_isPDNODE;
		int numNodeFam, famkk;
		int eqIndex_row[2], eqIndex_col[2];
		for (int bL = startP; bL < endP; bL++)
		{
			o_dat.op_getPDBE(bL)->getNodeID(bLinNID);
			bL1_isPDNODE = o_dat.op_getNode(bLinNID[0] - 1)->getNodeType();
			bL2_isPDNODE = o_dat.op_getNode(bLinNID[1] - 1)->getNodeType();
			if (bL1_isPDNODE == false || bL2_isPDNODE == false)
			{
				cout << "Warning: This node of PD domain boundary is not PD node!" << endl;
				exit(0);
			}
			o_dat.op_getNode(bLinNID[0] - 1)->getcoor(xL1);
			o_dat.op_getNode(bLinNID[1] - 1)->getcoor(xL2);
			nx = xL2[1] - xL1[1];
			ny = xL1[0] - xL2[0];
			N->setCoeff(0, 0, nx);
			N->setCoeff(0, 2, ny);
			N->setCoeff(1, 1, ny);
			N->setCoeff(1, 2, nx);
			//=============t1=================
			famkk = o_dat.op_getNode(bLinNID[0] - 1)->getFamID() - 1;
			temP_fami = o_dat.op_getFami(famkk);
			numNodeFam = temP_fami->getNumNode();
			C = new Matrix(3, 2 * numNodeFam);
			DC = new Matrix(3, 2 * numNodeFam);
			NDC = new Matrix(2, 2 * numNodeFam);
			matC2D(C, temP_fami, o_dat);
			matoperat.matMultiply(cop_D, C, DC);
			matoperat.matMultiply(N, DC, NDC);
			//u1;
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[0] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							temp = 1.0 / 3 * NDC->d_getCoeff(i, 2 * m + j);
							if (cb_InteralForce==false)
							{
								//non-dynamical solver
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								else
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
							}
							else
							{
								//dynamic solver;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
							}
							
						}
					}
				}

			}
			//u2; node 2
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[1] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							temp = 1.0 / 6 * NDC->d_getCoeff(i, 2 * m + j);
							if (cb_InteralForce == false)
							{
								//non-dynamic solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								else
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
							}
							else
							{
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
							}
							
						}
					}
				}

			}
			delete C, delete DC, delete NDC;
			C = NULL; DC = NULL; NDC = NULL;
			//===========t2**** node 2========
			famkk = o_dat.op_getNode(bLinNID[1] - 1)->getFamID() - 1;
			temP_fami = o_dat.op_getFami(famkk);
			numNodeFam = temP_fami->getNumNode();
			C = new Matrix(3, 2 * numNodeFam);
			DC = new Matrix(3, 2 * numNodeFam);
			NDC = new Matrix(2, 2 * numNodeFam);
			matC2D(C, temP_fami, o_dat);
			matoperat.matMultiply(cop_D, C, DC);
			matoperat.matMultiply(N, DC, NDC);
			//u1; node 1
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[0] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							
							temp = 1.0 / 6 * NDC->d_getCoeff(i, 2 * m + j);
							if (cb_InteralForce == false)
							{
								//non-dynamic solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								else
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
							}
							else
							{
								//dynamic solver;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
							}
							
						}
					}
				}

			}
			//u2; node 2===========
			for (int i = 0; i < 2; i++)
			{
				eqIndex_row[i] = o_dat.op_getNode(bLinNID[1] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex_row[i] != -1)
				{
					for (int m = 0; m < numNodeFam; m++)
					{
						NID_m = temP_fami->getNodeID(m);
						for (int j = 0; j < 2; j++)
						{
							temp = 1.0 / 3 * NDC->d_getCoeff(i, 2 * m + j);
							if (cb_InteralForce == false)
							{
								//non-dynamic solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								else
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
							}
							else
							{
								// dynamic solver;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
							}
						}
					}
				}

			}
			delete C, delete DC, delete NDC;
			C = NULL; DC = NULL; NDC = NULL;
		}
		delete N;
		N = NULL;
	}
	else if (numDime == 3)
	{
		if (o_dat.ci_eleType==12)
		{
			assemblePDBEworkQuad_CSRformat(o_dat);
		}
		else if (o_dat.ci_eleType==10)
		{
			assemblePDBEworkTetrahe_CSRformat(o_dat);
		}
		else
		{
			printf("ERROR: element Type for 3D must be 10 or 12.\n");
			exit(0);
		}
	}
}

void pdsolve::assemblePDBEworkQuad_CSRformat(datModel& o_dat)
{
	int numPDBEs, startP, endP;
	numPDBEs = o_dat.getTotnumPDBEs();
	startP = ci_rank * numPDBEs / ci_numProce;
	endP = (ci_rank + 1) * numPDBEs / ci_numProce;
	int numDime = o_dat.ci_Numdimen;
	pdFamily* temP_fami;
	
	//====================3D===================
	int nG, conNID[4];
	int famiID, numNodeFam, eqIndex_row, eqIndex_col;
	int NID_m;
	long long int i_temp;
	double p, q, wp, wq, fac, xN[4][3], temp, tempu;
	double N[4];//N[4] are shape functions;
	Matrix* mNt, * Nmat, * C, * f_mNt, * f_mNt_Nm, * f_mNt_Nm_D, * finMat;
	mNt = new Matrix(12, 3);
	Nmat = new Matrix(3, 6);
	f_mNt = new Matrix(12, 3);
	f_mNt_Nm = new Matrix(12, 6);
	f_mNt_Nm_D = new Matrix(12, 6);
	nG = o_globGP.i_getNumPts();
	//===========3D=======================================
	for (int pdbe = startP; pdbe < endP; pdbe++)
	{
		//==get PDBEs==
		o_dat.op_getPDBE(pdbe)->getNodeID(conNID);
		for (int i = 0; i < 4; i++)
		{
			o_dat.op_getNode(conNID[i] - 1)->getcoor(xN[i]);
		}
		//===Gauss integration===
		for (int mp = 0; mp < nG; mp++)
		{
			p = o_globGP.d_getGaussPt(mp);
			wp = o_globGP.d_getWeight(mp);
			for (int mq = 0; mq < nG; mq++)
			{
				q = o_globGP.d_getGaussPt(mq);
				wq = o_globGP.d_getWeight(mq);
				shapeFunctionQuad4N(N, p, q);
				matMathcalNt(mNt, p, q);
				for (int ei = 0; ei < 4; ei++)
				{
					//====get the final matrix for assembling;
					fac = N[ei] * wp * wq;
					matoperat.matMultiply(mNt, fac, f_mNt);
					matN_trans(Nmat, xN, p, q);
					famiID = o_dat.op_getNode(conNID[ei] - 1)->getFamID();
					temP_fami = o_dat.op_getFami(famiID - 1);
					numNodeFam = temP_fami->getNumNode();
					C = new Matrix(6, 3 * numNodeFam);
					finMat = new Matrix(12, 3 * numNodeFam);
					matC3D(C, temP_fami, o_dat);
					matoperat.matMultiply(f_mNt, Nmat, f_mNt_Nm);
					matoperat.matMultiply(f_mNt_Nm, cop_D, f_mNt_Nm_D);
					matoperat.matMultiply(f_mNt_Nm_D, C, finMat);
					//=======assembling======
					for (int i = 0; i < 4; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							eqIndex_row = o_dat.op_getNode(conNID[i] - 1)->op_getDof(j)->i_getEqInde();
							if (eqIndex_row != -1)
							{
								for (int m = 0; m < numNodeFam; m++)
								{
									NID_m = temP_fami->getNodeID(m);
									for (int jj = 0; jj < 3; jj++)
									{

										temp = finMat->d_getCoeff(i * 3 + j, 3 * m + jj);
										if (cb_InteralForce == false)
										{
											//non-dynamic solver
											eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
											if (eqIndex_col != -1)
											{
												i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
												cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
											}
											else
											{
												tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
												cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
											}
										}
										else
										{
											//dynamic solver
											tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
											cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
										}

									}
								}
							}
						}
					}
					//===finish assembling for this ei-th node;
					// delete;
					delete C, delete finMat;
					C = NULL, finMat = NULL;
				}
			}
		}
	}
	delete mNt, delete Nmat, delete f_mNt, delete f_mNt_Nm, delete f_mNt_Nm_D;
	mNt = NULL, Nmat = NULL, f_mNt = NULL, f_mNt_Nm = NULL, f_mNt_Nm_D = NULL;
	
}

void pdsolve::assemblePDBEworkTetrahe_CSRformat(datModel& o_dat)
{
	int numPDBEs, startP, endP;
	numPDBEs = o_dat.getTotnumPDBEs();
	startP = ci_rank * numPDBEs / ci_numProce;
	endP = (ci_rank + 1) * numPDBEs / ci_numProce;
	int numDime = o_dat.ci_Numdimen;
	pdFamily* temP_fami;
	
	//====================3D===================
	int conNID[4];
	int famiID, numNodeFam, eqIndex_row, eqIndex_col;
	int NID_m;
	long long int i_temp;
	double  xN[4][3], temp, tempu;
	Matrix* mNt[3], * mNtN[3], * C, * mNt_N_D, * finMat;
	for (int i = 0; i < 3; i++)
	{
		mNt[i] = new Matrix(9, 3);
		mNtN[i] = new Matrix(9, 6);
	}
	matMathcalNt_trian(mNt);
	mNt_N_D = new Matrix(9, 6);
	//===========3D=======================================
	for (int pdbe = startP; pdbe < endP; pdbe++)
	{
		//==get PDBEs==
		o_dat.op_getPDBE(pdbe)->getNodeID(conNID);
		for (int i = 0; i < 3; i++)
		{
			o_dat.op_getNode(conNID[i] - 1)->getcoor(xN[i]);
		}
		//=== integration===
		matMathcalNtN(mNtN, mNt, xN);
		for (int ei = 0; ei < 3; ei++)
		{
			matoperat.matMultiply(mNtN[ei], cop_D, mNt_N_D);
			famiID = o_dat.op_getNode(conNID[ei] - 1)->getFamID();
			temP_fami = o_dat.op_getFami(famiID - 1);
			numNodeFam = temP_fami->getNumNode();
			C = new Matrix(6, 3 * numNodeFam);
			finMat = new Matrix(9, 3 * numNodeFam);
			matC3D(C, temP_fami, o_dat);
			matoperat.matMultiply(mNt_N_D, C, finMat);
			//=======assembling======
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					eqIndex_row = o_dat.op_getNode(conNID[i] - 1)->op_getDof(j)->i_getEqInde();
					if (eqIndex_row != -1)
					{
						for (int m = 0; m < numNodeFam; m++)
						{
							NID_m = temP_fami->getNodeID(m);
							for (int jj = 0; jj < 3; jj++)
							{

								temp = finMat->d_getCoeff(i * 3 + j, 3 * m + jj);
								if (cb_InteralForce == false)
								{
									//non-dynamic solver
									eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
									if (eqIndex_col != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
									else
									{
										tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
										cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
									}
								}
								else
								{
									//dynamic solver
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
									cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
								}

							}
						}
					}
				}
			}
			//===finish assembling for this ei-th node;
			// delete;
			delete C, delete finMat;
			C = NULL, finMat = NULL;
		}
	}
	for (int i = 0; i < 3; i++)
	{
		delete mNt[i], delete mNtN[i];
		mNt[i]=NULL, mNtN[i]=NULL;
	}
	delete  mNt_N_D;
	mNt_N_D = NULL;
}

void pdsolve::assembleMassMatPD_CSRformat(datModel& o_dat)
{
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;

	int numDime = o_dat.ci_Numdimen;
	int Nid_k, eqIndex;
	double rho, dv;
	long long int i_temp;
	rho = o_dat.op_getmaterial()->getrho();
	for (int famkk = startP; famkk < endP; famkk++)
	{
		Nid_k = o_dat.op_getFami(famkk)->getNodeID(0);
		dv = o_dat.op_getNode(Nid_k - 1)->getvolume();
		for (int i = 0; i < numDime; i++)
		{
			eqIndex = o_dat.op_getNode(Nid_k - 1)->op_getDof(i)->i_getEqInde();
			if (eqIndex != -1)
			{
				i_temp = findCSRIndexOfMat(eqIndex, eqIndex);
				cdp_M[i_temp] = cdp_M[i_temp] + dv * rho;
			}
		}
	}
}

void pdsolve::assembleMassMatPD(datModel & o_dat)
{
	int numDime = o_dat.ci_Numdimen;
	int Nid_k, eqIndex;
	double rho, dv;
	rho = o_dat.op_getmaterial()->getrho();
	for (int famk = 0; famk < o_dat.getTotnumFami(); famk++)
	{
		Nid_k = o_dat.op_getFami(famk)->getNodeID(0);
		dv = o_dat.op_getNode(Nid_k - 1)->getvolume();
		for (int i = 0; i < numDime; i++)
		{
			eqIndex = o_dat.op_getNode(Nid_k - 1)->op_getDof(i)->i_getEqInde();
			if (eqIndex!=-1)
			{
				cop_M->addCoeff(eqIndex, eqIndex, dv*rho);
			}
		}
	}
}

void pdsolve::assembleLumpedMass(datModel& o_dat, int numEqua)
{
	cdp_M = new double[numEqua];
	cdp_MGlo = new double[numEqua];
	for (long long int i = 0; i < numEqua; i++)
	{
		cdp_M[i] = 0;
		cdp_MGlo[i] = 0;
	}
	//=======
	int totNumEle, startP, endP;
	totNumEle = o_dat.getTotnumEle();
	startP = ci_rank * totNumEle / ci_numProce;
	endP = (ci_rank + 1) * totNumEle / ci_numProce;
	//==
	int numDime, eqIndex, numNodeELE, *conNID;
	double(*xN)[3],rho,dv, VolEle = 0;
	numDime = o_dat.ci_Numdimen;
	rho = o_dat.op_getmaterial()->getrho();
	for (int ele = startP; ele < endP; ele++)
	{
		numNodeELE = o_dat.op_getEles(ele)->ci_numNodes;
		conNID = new int[numNodeELE];
		xN = new double[numNodeELE][3];
		o_dat.op_getEles(ele)->getConNid(conNID);
		for (int nd = 0; nd < numNodeELE; nd++)
		{
			o_dat.op_getNode(conNID[nd] - 1)->getcoor(xN[nd]);
		}
		//Gauss integration to get element volumn
		VolEle = o_dat.op_getEles(ele)->eleVolume(xN);
		dv = VolEle / numNodeELE;
		for (int nd = 0; nd < numNodeELE; nd++)
		{
			for (int i = 0; i < numDime; i++)
			{
				eqIndex = o_dat.op_getNode(conNID[nd] - 1)->op_getDof(i)->i_getEqInde();
				if (eqIndex != -1)
				{
					cdp_M[eqIndex] = cdp_M[eqIndex] + dv * rho;
				}
			}
		}
		delete[] conNID, delete[] xN;
		conNID = NULL, xN = NULL;

	}
	//===
	MPI_Allreduce(cdp_M, cdp_MGlo, numEqua, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	delete[] cdp_M;
	cdp_M = NULL;
}

void pdsolve::matMathcalNt(Matrix* mNt, double p, double q)
{
	// this function is to get the shape function matrix transpos for quad 4N PDBE;
	// the matrix is named as mathcal{N};
	double N[4];
	shapeFunctionQuad4N(N, p, q);
	mNt->zero();
	for (int k = 0; k < 4; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			mNt->setCoeff(i + 3 * k, i, N[k]);
		}
	}
	
}

void pdsolve::matN_trans(Matrix* Nmat, double xN[][3], double p, double q)
{
	// For each quad 4N PDBE;
	// This function is to get the matrix for transform ds in global cooridinate..
	// to local coordi dpdq;
	double dNdp[4], dNdq[4];
	dNdp[0] = 0.25 * (-1 + q);
	dNdp[1] = 0.25 * (1 - q);
	dNdp[2] = 0.25 * (1 + q);
	dNdp[3] = 0.25 * (-1 - q);
	dNdq[0] = 0.25 * (-1 + p);
	dNdq[1] = 0.25 * (-1 - p);
	dNdq[2] = 0.25 * (1 + p);
	dNdq[3] = 0.25 * (1 - p);
	double dXdp[3] = { 0 }, dXdq[3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			dXdp[i] = dXdp[i] + dNdp[j] * xN[j][i];
			dXdq[i] = dXdq[i] + dNdq[j] * xN[j][i];
		}
	}
	double N1, N2, N3;
	N1 = dXdp[1] * dXdq[2] - dXdp[2] * dXdq[1];
	N2 = dXdq[0] * dXdp[2] - dXdp[0] * dXdq[2];
	N3 = dXdp[0] * dXdq[1] - dXdq[0] * dXdp[1];
	Nmat->zero();
	Nmat->setCoeff(0, 0, N1);
	Nmat->setCoeff(1, 1, N2);
	Nmat->setCoeff(2, 2, N3);
	Nmat->setCoeff(0, 3, N2);
	Nmat->setCoeff(0, 5, N3);
	Nmat->setCoeff(1, 3, N1);
	Nmat->setCoeff(1, 4, N3);
	Nmat->setCoeff(2, 4, N2);
	Nmat->setCoeff(2, 5, N1);
}

void pdsolve::matMathcalNtN(Matrix* mNtN[], Matrix* mNt[], double xN[][3])
{
	Matrix* Nmat = new Matrix(3, 6);
	double v1[3], v2[3];
	for (int i = 0; i < 3; i++)
	{
		v1[i] = xN[1][i] - xN[0][i];
		v2[i] = xN[2][i] - xN[0][i];
	}
	double N1, N2, N3;
	N1 = v1[1] * v2[2] - v1[2] * v2[1];
	N2 = v2[0] * v1[2] - v1[0] * v2[2];
	N3 = v1[0] * v2[1] - v2[0] * v1[1];
	Nmat->zero();
	Nmat->setCoeff(0, 0, N1);
	Nmat->setCoeff(1, 1, N2);
	Nmat->setCoeff(2, 2, N3);
	Nmat->setCoeff(0, 3, N2);
	Nmat->setCoeff(0, 5, N3);
	Nmat->setCoeff(1, 3, N1);
	Nmat->setCoeff(1, 4, N3);
	Nmat->setCoeff(2, 4, N2);
	Nmat->setCoeff(2, 5, N1);
	for (int i = 0; i < 3; i++)
	{
		matoperat.matMultiply(mNt[i], Nmat, mNtN[i]);
	}
	delete Nmat;
	Nmat = NULL;
}

void pdsolve::matMathcalNt_trian(Matrix* mNt[])
{
	for (int i = 0; i < 3; i++)
	{
		mNt[i]->zero();
	}
	for (int i = 0; i < 3; i++)
	{
		//==
		mNt[0]->setCoeff(i, i, 1.0 / 12);
		mNt[0]->setCoeff(i + 3, i, 1.0 / 24);
		mNt[0]->setCoeff(i + 6, i, 1.0 / 24);
		//==
		mNt[1]->setCoeff(i, i, 1.0 / 24);
		mNt[1]->setCoeff(i + 3, i, 1.0 / 12);
		mNt[1]->setCoeff(i + 6, i, 1.0 / 24);
		//==
		mNt[2]->setCoeff(i, i, 1.0 / 24);
		mNt[2]->setCoeff(i + 3, i, 1.0 / 24);
		mNt[2]->setCoeff(i + 6, i, 1.0 / 12);
	}
}

void pdsolve::assembleSEDbyFEM(datModel & o_dat)
{
	// By traditional FEM
	Matrix *Ke;
	int totNumEle, algoType, numDime;
	totNumEle = o_dat.getTotnumEle();
	numDime = o_dat.ci_Numdimen;

	int numNode;
	int *conNid, NID_i, NID_j, eqInde_i, eqInde_j;
	double temp, (*xN)[3], tempu;
	
	for (int ele = 0; ele < totNumEle; ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		if (algoType == 2)
		{
			numNode = o_dat.op_getEles(ele)->getNumNodes();
			conNid = new int[numNode];
			xN = new double[numNode][3];
			Ke = new Matrix(numDime * numNode, numDime * numNode);
			o_dat.op_getEles(ele)->getConNid(conNid);
			for (int  nd = 0; nd < numNode; nd++)
			{
				o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			}
			o_dat.op_getEles(ele)->eleStiffMatFEM(Ke, cop_D, xN);
			for (int i = 0; i < numNode; i++)
			{
				NID_i = conNid[i];
				for (int ii = 0; ii < numDime; ii++)
				{
					eqInde_i = o_dat.op_getNode(NID_i - 1)->op_getDof(ii)->i_getEqInde();
					if (eqInde_i != -1)
					{
						for (int j = 0; j < numNode; j++)
						{
							NID_j = conNid[j];
							for (int jj = 0; jj < numDime; jj++)
							{
								temp = Ke->d_getCoeff(numDime * i + ii, numDime * j + jj);
								eqInde_j = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->i_getEqInde();
								if (eqInde_j != -1)
								{
									cop_Ku->addCoeff(eqInde_i, eqInde_j, temp);
								}
								else
								{
									tempu = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->d_getValue();
									cop_F->addCoeff(eqInde_i, -temp * tempu);
								}
							}
						}
					}
				}
			}
			delete[] conNid, delete[] xN;
			conNid = NULL, xN = NULL;
			delete Ke;
			Ke = NULL;
		}

	}
}

void pdsolve::assembleSEDbyFEM_CSRformat(datModel& o_dat)
{
	// By traditional FEM
		
	int totNumFE, startP, endP;
	totNumFE = o_dat.civ_feIDX.size();
	startP = ci_rank * totNumFE / ci_numProce;
	endP = (ci_rank + 1) * totNumFE / ci_numProce;

	Matrix* Ke;
	int algoType, numDime;
	numDime = o_dat.ci_Numdimen;

	int numNode, ele;
	long long int i_temp;
	int* conNid, NID_i, NID_j, eqInde_i, eqInde_j;
	double temp, (*xN)[3], tempu;

	
	for (int fe = startP; fe < endP; fe++)
	{
		ele = o_dat.civ_feIDX[fe];
		//algoType = o_dat.op_getEles(ele)->getAlgoType();
		//if (algoType == 2)
		//{
			numNode = o_dat.op_getEles(ele)->getNumNodes();
			conNid = new int[numNode];
			xN = new double[numNode][3];
			Ke = new Matrix(numDime * numNode, numDime * numNode);
			o_dat.op_getEles(ele)->getConNid(conNid);
			for (int nd = 0; nd < numNode; nd++)
			{
				o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			}
			o_dat.op_getEles(ele)->eleStiffMatFEM(Ke, cop_D, xN);
			for (int i = 0; i < numNode; i++)
			{
				NID_i = conNid[i];
				for (int ii = 0; ii < numDime; ii++)
				{
					eqInde_i = o_dat.op_getNode(NID_i - 1)->op_getDof(ii)->i_getEqInde();
					if (eqInde_i != -1)
					{
						for (int j = 0; j < numNode; j++)
						{
							NID_j = conNid[j];
							for (int jj = 0; jj < numDime; jj++)
							{
								
								temp = Ke->d_getCoeff(numDime * i + ii, numDime * j + jj);
								if (cb_InteralForce == false)
								{
									//cout << "== " << ci_solvFlag << endl;
									//non-dynamic solver;
									eqInde_j = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->i_getEqInde();
									if (eqInde_j != -1)
									{
										i_temp = findCSRIndexOfMat(eqInde_i, eqInde_j);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
									else
									{
										tempu = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->d_getValue();
										cdp_F[eqInde_i] = cdp_F[eqInde_i] - temp * tempu;
									}
								}
								else
								{
									//dynamic solver;
									tempu = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->d_getValue();
									cdp_F[eqInde_i] = cdp_F[eqInde_i] - temp * tempu;
								}
								
							}
						}
					}
				}
			}
			delete[] conNid, delete[] xN;
			conNid = NULL, xN = NULL;
			delete Ke;
			Ke = NULL;
		//}

	}
}

void pdsolve::calExternalForce_CSRformat(datModel& o_dat)
{
	//====equivalent force from essential BCs
	//==============cal nodal force from point BC;
	int numPBCs, startP, endP;
	numPBCs = o_dat.getTotnumPointBCs();
	startP = ci_rank * numPBCs / ci_numProce;
	endP = (ci_rank + 1) * numPBCs / ci_numProce;

	int NID, fDof, eqIndex, numDime = o_dat.ci_Numdimen;
	double vForce, temp;
	for (int i = startP; i < endP; i++)
	{
		NID = o_dat.op_getPointBC(i)->i_getNID();
		fDof = o_dat.op_getPointBC(i)->i_getfDOF();
		vForce = o_dat.op_getPointBC(i)->d_getValue();
		eqIndex = o_dat.op_getNode(NID - 1)->op_getDof(fDof)->i_getEqInde();
		cdp_F[eqIndex] = cdp_F[eqIndex] + vForce;
	}

	//==================cal Nodal force from Nature force;
	//=====3D
	int numTracBC = o_dat.getTotnumNaturalBCs();
	
	double tVal, (*xN)[3];
	int numEle, numNoEle, * conNid;
	Vector* Fe;
	for (int nBC = 0; nBC < numTracBC; nBC++)
	{
		tVal = o_dat.op_getNaturalBC(nBC)->getValue();
		numEle = o_dat.op_getNaturalBC(nBC)->getNumEle();
		startP = ci_rank * numEle / ci_numProce;
		endP = (ci_rank + 1) * numEle / ci_numProce;
		for (int ele = startP; ele < endP; ele++)
		{
			numNoEle = o_dat.op_getNaturalBC(nBC)->op_getNBCsEle(ele)->getNumNodes();
			xN = new double[numNoEle][3];
			conNid = new int[numNoEle];
			Fe = new Vector(numDime * numNoEle);
			o_dat.op_getNaturalBC(nBC)->op_getNBCsEle(ele)->getConNid(conNid);
			for (int nd = 0; nd < numNoEle; nd++)
			{
				o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			}
			o_dat.op_getNaturalBC(nBC)->op_getNBCsEle(ele)->eleEquivNodalForce(Fe, tVal, xN);
			for (int nd = 0; nd < numNoEle; nd++)
			{
				for (int ii = 0; ii < numDime; ii++)
				{
					eqIndex = o_dat.op_getNode(conNid[nd] - 1)->op_getDof(ii)->i_getEqInde();
					if (eqIndex != -1)
					{
						temp = Fe->d_getCoeff(numDime * nd + ii);
						cdp_F[eqIndex] = cdp_F[eqIndex] + temp;
					}
				}
			}
			delete[] xN, delete[] conNid;
			xN = NULL, conNid = NULL;
			delete Fe; Fe = NULL;
		}
	}
}

void pdsolve::pdfemSolver(datModel& o_dat, fioFiles& o_files)
{
	if (o_dat.ci_solvFlag==1)
	{
		//static solver
		pdfemStaticSolver_CSRformat(o_dat, o_files);
	}
	else if (o_dat.ci_solvFlag==0)
	{
		//dynamic solver;
		pdfemDynamicSolver_CSRformat(o_dat, o_files);
	}
	else if (o_dat.ci_solvFlag == 2)
	{
		//quasi static solver;
		pdfemQuasiStaticSolver_CSRformat(o_dat, o_files);
	}
	else
	{
		printf("ERROR: solver flag is wrong\n");
		exit(0);
	}
}

void pdsolve::calAcceleration(int numEq, double* dp_A)
{
	for (int i = 0; i < numEq; i++)
	{
		dp_A[i] = cdp_FGlo[i] / cdp_MGlo[i];
	}
}

void pdsolve::timeIntegration(datModel& o_dat,Vector* Vu_n, Vector* Vu_nm1, Vector* Vu_np1, int numEq)
{
	int  startP, endP;
	startP = ci_rank * numEq / ci_numProce;
	endP = (ci_rank + 1) * numEq / ci_numProce;
	double tempA;
	double dt = o_dat.getTstep();
	for (int i = startP; i < endP; i++)
	{
		tempA = dt * dt * cdp_FGlo[i] / cdp_MGlo[i];
		Vu_np1->cdp_vecCoeff[i] = Vu_n->cdp_vecCoeff[i] * 2.0 + tempA - Vu_nm1->cdp_vecCoeff[i];
	}
	MPI_Bcast(&(Vu_np1->cdp_vecCoeff[startP]), endP - startP, MPI_DOUBLE, ci_rank, MPI_COMM_WORLD);
}

void pdsolve::assembleElemassMatFEM_CSRformat(datModel& o_dat)
{
	int totNumFE, startP, endP;
	totNumFE = o_dat.civ_feIDX.size();
	startP = ci_rank * totNumFE / ci_numProce;
	endP = (ci_rank + 1) * totNumFE / ci_numProce;

	int numNode, numDime = o_dat.ci_Numdimen;
	Matrix* Me;
	int * conNid, eqIndex_row, eqIndex_col, ele;
	double temp, (*xN)[3], rho;
	rho = o_dat.op_getmaterial()->getrho();
	long long int i_temp;
	for (int fe = startP; fe < endP; fe++)
	{
		ele = o_dat.civ_feIDX[fe];
		//if (algoType == 2)
		//{
			numNode = o_dat.op_getEles(ele)->getNumNodes();
			conNid = new int[numNode];
			xN = new double[numNode][3];
			Me = new Matrix(numNode * numDime, numNode * numDime);
			o_dat.op_getEles(ele)->getConNid(conNid);
			for (int nd = 0; nd < numNode; nd++)
			{
				o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			}
			o_dat.op_getEles(ele)->eleMassMat(Me, rho, xN);
			for (int i = 0; i < numNode; i++)
			{
				for (int ii = 0; ii < numDime; ii++)
				{
					eqIndex_row = o_dat.op_getNode(conNid[i] - 1)->op_getDof(ii)->i_getEqInde();
					if (eqIndex_row != -1)
					{
						for (int j = 0; j < numNode; j++)
						{
							for (int jj = 0; jj < numDime; jj++)
							{
								eqIndex_col = o_dat.op_getNode(conNid[j] - 1)->op_getDof(jj)->i_getEqInde();
								if (eqIndex_col != -1)
								{
									temp = Me->d_getCoeff(numDime * i + ii, numDime * j + jj);
									i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
									cdp_M[i_temp] = cdp_M[i_temp] + temp;
									//cop_M->addCoeff(eqIndex_row, eqIndex_col, temp);
								}
							}
						}
					}
				}
			}
			delete[] xN, delete[] conNid; xN = NULL; conNid = NULL;
			delete Me; Me = NULL;
		//}
	}
}

void pdsolve::assembleElemassMatFEM(datModel & o_dat, ofstream & test)
{
	int numNode, numDime = o_dat.ci_Numdimen;
	Matrix *Me;
	//Me = new Matrix(8, 8);
	int algoType, *conNid, eqIndex_row, eqIndex_col;
	double temp, (*xN)[3], rho;
	rho = o_dat.op_getmaterial()->getrho();
	for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		if (algoType==2)
		{
			numNode = o_dat.op_getEles(ele)->getNumNodes();
			conNid = new int[numNode];
			xN = new double[numNode][3];
			Me = new Matrix(numNode * numDime, numNode * numDime);
			o_dat.op_getEles(ele)->getConNid(conNid);
			for (int nd = 0; nd < numNode; nd++)
			{
				o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			}
			o_dat.op_getEles(ele)->eleMassMat(Me, rho, xN);
			for (int i = 0; i < numNode; i++)
			{
				for (int ii = 0; ii < numDime; ii++)
				{
					eqIndex_row = o_dat.op_getNode(conNid[i] - 1)->op_getDof(ii)->i_getEqInde();
					if (eqIndex_row!=-1)
					{
						for (int j = 0; j < numNode; j++)
						{
							for (int jj = 0; jj < numDime; jj++)
							{
								eqIndex_col = o_dat.op_getNode(conNid[j] - 1)->op_getDof(jj)->i_getEqInde();
								if (eqIndex_col!=-1)
								{
									temp = Me->d_getCoeff(numDime * i + ii, numDime * j + jj);
									cop_M->addCoeff(eqIndex_row, eqIndex_col, temp);
								}
							}
						}
					}
				}
			}
			delete[] xN, delete[] conNid; xN = NULL; conNid = NULL;
			delete Me; Me = NULL;
		}
	}
}

void pdsolve::setVaryEssentialBC(datModel & o_dat)
{
	// displacement increasing boundary condition, for dynamic or quasi-static solver;
	double dt, Uval;
	dt = o_dat.getTstep();
	int numEBCs, numNODE, i_dof, Nid;
	numEBCs = o_dat.getTotnumEssentialBCs();
	for (int i = 0; i < numEBCs; i++)
	{
		if (o_dat.op_getEssenBC(i)->cb_varing)
		{
			Uval = o_dat.op_getEssenBC(i)->cd_value + dt * o_dat.op_getEssenBC(i)->cd_velocity;
			o_dat.op_getEssenBC(i)->cd_value = Uval;
			numNODE = o_dat.op_getEssenBC(i)->getNumNODE();
			i_dof = o_dat.op_getEssenBC(i)->getDOF();
			for (int j = 0; j < numNODE; j++)
			{
				Nid = o_dat.op_getEssenBC(i)->cip_NID[j];
				o_dat.op_getNode(Nid-1)->op_getDof(i_dof)->setValue(Uval);
			}
		}
	}
}

void pdsolve::setVaryNaturalBC(datModel& o_dat)
{
	//  increasing natural boundary condition, for dynamic or quasi-static solver;

	double fval, dt;
	dt = o_dat.getTstep();
	int numNBCs = o_dat.getTotnumNaturalBCs();
	for (int i = 0; i < numNBCs; i++)
	{
		if (o_dat.op_getNaturalBC(i)->cb_varing)
		{
			fval = o_dat.op_getNaturalBC(i)->cd_value + dt * o_dat.op_getNaturalBC(i)->cd_velocity;
			o_dat.op_getNaturalBC(i)->cd_value = fval;
		}
	}
}


void pdsolve::resetDispBC(datModel& o_dat,double Multip)
{
	int Nid, direc;
	double xc[2], tempu;
	for (int k = 0; k < o_dat.getTotnumVaryEssenBC(); k++)
	{
		Nid = o_dat.op_getVaryEssenBC(k)->getNid();
		direc = o_dat.op_getVaryEssenBC(k)->getDirec();
		o_dat.op_getNode(Nid - 1)->getcoor(xc);
		tempu = o_dat.op_getNode(Nid - 1)->op_getDof(direc)->d_getValue();
		tempu = tempu * Multip;
		o_dat.op_getNode(Nid - 1)->op_getDof(direc)->setValue(tempu);
	}
}


void pdsolve::calExternalForce( datModel & o_dat)
{
	//====equivalent force from essential BCs===
	//************************************
	//==============cal nodal force from point BC;
	int NID, fDof, eqIndex, numDime = o_dat.ci_Numdimen;
	double vForce, temp;
	for (int i = 0; i < o_dat.getTotnumPointBCs(); i++)
	{
		NID = o_dat.op_getPointBC(i)->i_getNID();
		fDof = o_dat.op_getPointBC(i)->i_getfDOF();
		vForce = o_dat.op_getPointBC(i)->d_getValue();
		eqIndex = o_dat.op_getNode(NID - 1)->op_getDof(fDof)->i_getEqInde();
		cop_F->addCoeff(eqIndex, vForce);
	}

	//==================cal Nodal force from Nature force;
	//=====3D
	int numTracBC = o_dat.getTotnumNaturalBCs();
	double tVal, (*xN)[3];
	int numEle, numNoEle, * conNid;
	Vector* Fe;
	for (int nBC = 0; nBC < numTracBC; nBC++)
	{
		tVal = o_dat.op_getNaturalBC(nBC)->getValue();
		numEle = o_dat.op_getNaturalBC(nBC)->getNumEle();
		for (int ele = 0; ele < numEle; ele++)
		{
			numNoEle = o_dat.op_getNaturalBC(nBC)->op_getNBCsEle(ele)->getNumNodes();
			xN = new double[numNoEle][3];
			conNid = new int[numNoEle];
			Fe = new Vector(numDime * numNoEle);
			o_dat.op_getNaturalBC(nBC)->op_getNBCsEle(ele)->getConNid(conNid);
			for (int nd = 0; nd < numNoEle; nd++)
			{
				o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			}
			o_dat.op_getNaturalBC(nBC)->op_getNBCsEle(ele)->eleEquivNodalForce(Fe, tVal, xN);
			for (int nd = 0; nd < numNoEle; nd++)
			{
				for (int ii = 0; ii < numDime; ii++)
				{
					eqIndex = o_dat.op_getNode(conNid[nd] - 1)->op_getDof(ii)->i_getEqInde();
					if (eqIndex != -1)
					{
						temp = Fe->d_getCoeff(numDime * nd + ii);
						cop_F->addCoeff(eqIndex, temp);
					}
				}
			}
			delete[] xN, delete[] conNid;
			xN = NULL, conNid = NULL;
			delete Fe; Fe = NULL;
		}
	}
}

void pdsolve::pdfemAssembleEquaSys(datModel& o_dat)
{

	printf("TEST==\n");
	cop_F->zero();

	cop_Ku->zero();
	cop_F->zero();
	assembleInterWorkPD(o_dat);
	assemblePDBEwork(o_dat);
	assembleSEDbyFEM(o_dat);

	
	calExternalForce(o_dat);
}

void pdsolve::pdfemAssembleEquasSys_CSRformat(datModel& o_dat, int numEq)
{
	cdp_F = new double[numEq];
	cdp_FGlo = new double[numEq];
	for (int i = 0; i < numEq; i++)
	{
		cdp_F[i] = 0;
		cdp_FGlo[i] = 0;
	}
	cb_InteralForce = false;
	setCSRIndexes_gloStiffMat(o_dat);
	assembleInterWorkPD_CSRformat(o_dat);
	assemblePDBEwork_CSRformat(o_dat);
	assembleSEDbyFEM_CSRformat(o_dat);
	calExternalForce_CSRformat(o_dat);
	MPI_Allreduce(cdp_Ku, cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce(cdp_Ku, cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	delete[] cdp_F, delete[] cdp_Ku;
	cdp_F = NULL, cdp_Ku = NULL;
	//MPI_Bcast(cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Bcast(cdp_FGlo, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void pdsolve::pdfemStaticSolver(datModel& o_dat)
{
	//Don't need block data any more for static solver;
	o_dat.deleteBLOCK();
	//================
	o_dat.ci_solvFlag = 1;// remove later;
	if (o_dat.ci_solvFlag!=1)
	{
		printf("ERROR: This is not static solver\n");
		printf("Please reset the solver flag\n");
		exit(0);
	}
	//================
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode, numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = totNumNode * numDime - numPreDof;
	cop_Ku = new Matrix(numEq, numEq);
	cop_F = new Vector(numEq);
	pdfemAssembleEquaSys(o_dat);
	ofstream fout("test11.out");
	cop_F->print(fout);
	fout.close();
	Vector *uu_n = new Vector(numEq);
	matoperat.PARDISOsolveSparse(cop_Ku, cop_F, uu_n);
	//store displacememts results;
	double tempu;
	int eqIndex;
	for (int k = 0; k < totNumNode; k++)
	{
		for (int i = 0; i < numDime; i++)
		{
			eqIndex = o_dat.op_getNode(k)->op_getDof(i)->i_getEqInde();
			if (eqIndex != -1)
			{
				tempu = uu_n->d_getCoeff(eqIndex);
				o_dat.op_getNode(k)->op_getDof(i)->setValue(tempu);
			}
		}
	}
	delete uu_n, delete cop_Ku, delete cop_F;
	uu_n=NULL, cop_Ku = NULL, cop_F = NULL;
	//==stresses
	calGlobalNodeStresses(o_dat);
}

void pdsolve::pdfemStaticSolver_CSRformat(datModel& o_dat, fioFiles& o_files)
{
	
	//Don't need block data any more for static solver;
	//o_dat.deleteBLOCK();
	//================
	o_dat.ci_solvFlag = 1;// remove later;
	if (o_dat.ci_solvFlag != 1)
	{
		printf("ERROR: This is not static solver\n");
		printf("Please reset the solver flag\n");
		exit(0);
	}
	//====
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode;
	int numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = (totNumNode * numDime) - (numPreDof);
	long long int numEq_long = numEq;
	if (ci_rank == 0)
	{
		printf("Assembling Equations.....\n");
	}
	pdfemAssembleEquasSys_CSRformat(o_dat,numEq);
	if (ci_rank == 0)
	{
		printf("Finished assembling Equations.\n");
	}
	
	int comm = MPI_Comm_c2f(MPI_COMM_WORLD);
	Vector* Ug = new Vector(numEq);
	if (ci_rank == 0)
	{
		printf("cluster_PARDISO static solving.....\n");
	}
	matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_KuGlo,
				cdp_FGlo, Ug->cdp_vecCoeff, &comm);
	if (ci_rank == 0)
	{
		printf("Finshed cluster_PARDISO static solving.\n");
	}
	/*store results========;*/
	MPI_Bcast(Ug->cdp_vecCoeff, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	storeDisplacementResult(o_dat, Ug);
	delete Ug; Ug = NULL;
	//sequential solve omp solver
	//if (ci_rank==0)
	//{
	//	//cdp_Ug = new double[numEq];
	//	Vector* Ug = new Vector(numEq);
	//	//mkl_set_dynamic(true);
	//	//mkl_set_num_threads(4);
	//	//mkl_domain_set_num_threads(4, MKL_DOMAIN_PARDISO);
	//	//omp_set_num_threads(1);
	//	matoperat.PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_KuGlo, cdp_FGlo,Ug->cdp_vecCoeff);
	//	/*int comm = MPI_Comm_c2f(MPI_COMM_WORLD);
	//	matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_KuGlo, 
	//		cdp_FGlo, Ug->cdp_vecCoeff, &comm);*/
	//	//store displacememts results;
	//	storeDisplacementResult(o_dat, Ug);
	//	delete Ug; Ug = NULL;
	//}
	delete[] cip_ia, delete[] cip_ja, delete[] cdp_FGlo, delete[] cdp_KuGlo;
	cip_ia = NULL, cip_ja = NULL, cdp_FGlo = NULL, cdp_KuGlo = NULL;
	/*==stresses==========*/
	calGlobalNodeStresses(o_dat);
	if (ci_rank == 0)
	{
		//===========post processing=====================
		printf("writing results.......\n");
		o_files.writeResults(o_dat);
	}
}

void pdsolve::storeDisplacementResult(datModel& o_dat, Vector* U)
{
	//store;
	int totNumNode = o_dat.getTotnumNode();
	int eqIndex, numDime = o_dat.ci_Numdimen;
	for (int k = 0; k < totNumNode; k++)
	{
		for (int i = 0; i < numDime; i++)
		{
			eqIndex = o_dat.op_getNode(k)->op_getDof(i)->i_getEqInde();
			if (eqIndex!= -1)
			{
				o_dat.op_getNode(k)->op_getDof(i)->setValue(U->d_getCoeff(eqIndex));
			}
		}
	}
}

void pdsolve::setCSRIndexes_gloMassMat(datModel& o_dat)
{
	/*This function is set indexes for global mass matrix in CSR format;
	output is ia and ja;*/
	
	//allocate ia zero-based format;
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode, numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = totNumNode * numDime - numPreDof;
	cip_ia = new long long int[numEq + 1];
	for (int i = 0; i < numEq + 1; i++)
	{
		//zero-based, ia[0]=0 always.
		cip_ia[i] = 0;
	}
	//============ find interactions of each node=========;
	vector<int>* inteNode;
	inteNode = new vector<int>[totNumNode];
	//====PD family===;
	int NID_k;
	for (int famk = 0; famk < o_dat.getTotnumFami(); famk++)
	{
		NID_k = o_dat.op_getFami(famk)->getNodeID(0);
		inteNode[NID_k - 1].push_back(NID_k); // for pd node, only itself has contribution to mass;
	}
	//===fem interactions;
	int algoType, numNele, * conNID;
	for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
	{
		algoType = o_dat.op_getEles(ele)->getAlgoType();
		if (algoType == 2)
		{
			numNele = o_dat.op_getEles(ele)->getNumNodes();
			conNID = new int[numNele];
			o_dat.op_getEles(ele)->getConNid(conNID);
			for (int nd = 0; nd < numNele; nd++)
			{
				for (int i = 0; i < numNele; i++)
				{
					inteNode[conNID[nd] - 1].push_back(conNID[i]);
				}
			}
			delete[]conNID;
			conNID = NULL;
		}
	}
	//=====delete the duplicated nodes and sort them in ascending order;
	unordered_set<int> inteSet;
	for (int nd = 0; nd < totNumNode; nd++)
	{
		//delete duplicted nodes
		for (int NODEID : inteNode[nd])
		{
			inteSet.insert(NODEID);
		}
		inteNode[nd].assign(inteSet.begin(), inteSet.end());
		//sort;
		sort(inteNode[nd].begin(), inteNode[nd].end());
		inteSet.clear();
	}
	//============get ia, ja==========================================
	int eqInd_row, eqInd_col;
	vector<int>vec_ja;
	for (int nd = 0; nd < totNumNode; nd++)
	{
		if (o_dat.op_getNode(nd)->getNodeType()!=2)
		{
			//if not pure pd node;
			for (int i = 0; i < numDime; i++)
			{
				eqInd_row = o_dat.op_getNode(nd)->op_getDof(i)->i_getEqInde();
				if (eqInd_row != -1)
				{
					for (int itN = 0; itN < int(inteNode[nd].size()); itN++)
					{
						for (int j = 0; j < numDime; j++)
						{
							eqInd_col = o_dat.op_getNode(inteNode[nd][itN] - 1)->op_getDof(j)->i_getEqInde();
							if (eqInd_col != -1)
							{
								vec_ja.push_back(eqInd_col);//zero-based;
								// cip_ia[0]=0, always for zero-based, so eqInd_row+1;
								//count the num of non-zero coeff of each row
								// will be revised later;
								cip_ia[eqInd_row + 1] = cip_ia[eqInd_row + 1] + 1;
							}
						}
					}
				}
			}

		}
		else
		{
			//pure PD node;
			for (int i = 0; i < numDime; i++)
			{
				eqInd_row = o_dat.op_getNode(nd)->op_getDof(i)->i_getEqInde();
				if (eqInd_row != -1)
				{
					vec_ja.push_back(eqInd_row);//zero-based; eqInd_row=eqInd_col
					// cip_ia[0]=0, always for zero-based, so eqInd_row+1;
					//count the num of non-zero coeff of each row
					// will be revised later;
					cip_ia[eqInd_row + 1] = cip_ia[eqInd_row + 1] + 1;
				}
			}

		}
	}
	for (int nd = 0; nd < totNumNode; nd++)
	{
		inteNode[nd].clear();
	}
	delete[] inteNode;
	inteNode = NULL;
	//====get real ia=== caution: zero-based;
	for (int i = 0; i < numEq; i++)
	{
		cip_ia[i + 1] = cip_ia[i] + cip_ia[i + 1];
	}
	//====store ja in 1D array===;
	long long int size_ja = vec_ja.size();
	cip_ja = new long long int[size_ja];
	for (int i = 0; i < size_ja; i++)
	{
		cip_ja[i] = vec_ja[i];
	}
	vec_ja.clear();
	//===initializa cdp_M;
	cdp_M = new double[size_ja];
	cdp_MGlo = new double[size_ja];
	for (long long int i = 0; i < size_ja; i++)
	{
		cdp_M[i] = 0;
		cdp_MGlo[i] = 0;
	}
}

void pdsolve::calinternalForce_CSRformat(datModel& o_dat, int numEq)
{
	if (!cdp_F)
	{
		cdp_F = new double[numEq];
	}
	if (!cdp_FGlo)
	{
		cdp_FGlo = new double[numEq];
	}
	for (int i = 0; i < numEq; i++)
	{
		cdp_F[i] = 0;
		cdp_FGlo[i] = 0;
	}
	cb_InteralForce = true;
	assembleInterWorkPD_CSRformat(o_dat);
	assemblePDBEwork_CSRformat(o_dat);
	assembleSEDbyFEM_CSRformat(o_dat);
	//calExternalForce_CSRformat(o_dat);
	//MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void pdsolve::pdfemDynamicSolver_CSRformat(datModel& o_dat, fioFiles& o_files)
{
	////===========static solver to set the start ebc of dynamic solver;
	//o_dat.ci_solvFlag = 1;
	//pdfemStaticSolver_CSRformat(o_dat,o_files); // stress calculated;
	//double maxSig = failureCriterion_stress(o_dat);
	//double ratio = o_dat.op_getmaterial()->getSigult() / maxSig*1.2;
	//double Uval;
	//int numEBCs, numNODE, i_dof, Nid;
	//numEBCs = o_dat.getTotnumEssentialBCs();
	////set start ebc;
	//for (int i = 0; i < numEBCs; i++)
	//{
	//	if (o_dat.op_getEssenBC(i)->cb_varing)
	//	{
	//		Uval = o_dat.op_getEssenBC(i)->cd_value * ratio;
	//		o_dat.op_getEssenBC(i)->cd_value = Uval;
	//		numNODE = o_dat.op_getEssenBC(i)->getNumNODE();
	//		i_dof = o_dat.op_getEssenBC(i)->getDOF();
	//		for (int j = 0; j < numNODE; j++)
	//		{
	//			Nid = o_dat.op_getEssenBC(i)->cip_NID[j];
	//			o_dat.op_getNode(Nid - 1)->op_getDof(i_dof)->setValue(Uval);
	//		}
	//	}
	//}
	////static solver for start ebc;
	//pdfemStaticSolver_CSRformat(o_dat, o_files); // stress calculated;
	//// updat bond-state
	//o_dat.ci_solvFlag = 0;
	//failureCriterion_stress(o_dat);
	//Don't need block data any more for dynamic solver;
	o_dat.deleteBLOCK();
	//================================================================
	o_dat.ci_solvFlag = 0;// remove later;
	if (o_dat.ci_solvFlag != 0)
	{
		printf("ERROR: This is not dynamic solver\n");
		printf("Please reset the solver flag\n");
		exit(0);
	}
	//============number of equations===================
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode;
	int numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = totNumNode * numDime - numPreDof;
	long long int numEq_long = numEq;
	//========================================================================
	//=============inital some variables=====
	/*Vu_nm1--Vector: displacement of (n-1)th step;
	 Vu_n--Vector: displacement of (n)th step;
	 Vu_np1--Vector: displacement of (n+1)th step;
	 Va_n--Vector: accerations of (n)th step;*/
	int saveFren = o_dat.getSaveFreq();
	int numTstep = o_dat.getnumTstep();
	double* dp_A = new double[numEq];// accerations, Va_n point to this pointer;
	Vector* Vu_nm1, * Vu_n, * Vu_np1, * Va_n;
	Vu_nm1 = new Vector(numEq);
	Vu_n = new Vector(numEq);
	Vu_np1 = new Vector(numEq);
	Va_n = new Vector();
	Va_n->setNumRows(numEq);
	Va_n->setCoeff(dp_A);

	//=====mass===
	if (o_dat.cb_lumpedMass)
	{
		assembleLumpedMass(o_dat, numEq);
	}
	else
	{
		//==========mass matrix CSR format======
		setCSRIndexes_gloMassMat(o_dat);
		assembleElemassMatFEM_CSRformat(o_dat);
		assembleMassMatPD_CSRformat(o_dat);
		MPI_Allreduce(cdp_M, cdp_MGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		delete[]cdp_M;
		cdp_M = NULL;
	}
	
	//===================
	//===============start dynamic solving==========
	double dt = o_dat.getTstep();
	if (ci_rank==0)
	{
		printf("start to Dynamic solving......\n");
		printf("Time step is %e\n", dt);
	}
	//============get 1th step displacement;
	//initial Vu_n========
	int eqIndex;
	for (int k = 0; k < totNumNode; k++)
	{
		for (int i = 0; i < numDime; i++)
		{
			eqIndex = o_dat.op_getNode(k)->op_getDof(i)->i_getEqInde();
			if (eqIndex != -1)
			{
				Vu_n->setCoeff(eqIndex, o_dat.op_getNode(k)->op_getDof(i)->d_getValue());
			}
		}
	}
	//==set varied prescribed displacemets;
	setVaryEssentialBC(o_dat);
	//====resultant force;
	calinternalForce_CSRformat(o_dat, numEq);
	calExternalForce_CSRformat(o_dat);
	MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//===solve for acceleration
	int comm = MPI_Comm_c2f(MPI_COMM_WORLD);
	if (o_dat.cb_lumpedMass)
	{
		for (int i = 0; i < numEq; i++)
		{
			dp_A[i] = cdp_FGlo[i] / cdp_MGlo[i];
		}
	}
	else
	{
		matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_MGlo,
			cdp_FGlo, dp_A, &comm);
		MPI_Bcast(dp_A, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	//===updat displacement and store;
	matoperat.matMultiply(Va_n, 0.5 * dt * dt, Vu_np1);
	matoperat.matAdd(Vu_np1, Vu_n, Vu_np1);
	storeDisplacementResult(o_dat, Vu_np1);
	//calculate stress==
	calGlobalNodeStresses(o_dat);
	if (o_dat.ci_failFlag==0)
	{
		//===updat bond-state==
		failureCriterion_stress(o_dat);
	}
	else if (o_dat.ci_failFlag == 1)
	{
		failureCriterion_stretch(o_dat);
	}
	
	
	//===
	//===start loop=============;
	for (int n = 1; n < numTstep+1; n++)
	{
		if (ci_rank==0)
		{
			cout << "time step n= " << n << endl;
		}
		//initial Vu_n, Vu_nm1;
		for (int i = 0; i < numEq; i++)
		{
			Vu_nm1->setCoeff(i, Vu_n->d_getCoeff(i));
			Vu_n->setCoeff(i, Vu_np1->d_getCoeff(i));
		}

		//==set varied prescribed displacemets;
		setVaryEssentialBC(o_dat);
		//====resultant force;
		calinternalForce_CSRformat(o_dat, numEq);
		calExternalForce_CSRformat(o_dat);
		MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//===time integration
		if (o_dat.cb_lumpedMass)
		{
			timeIntegration(o_dat, Vu_n, Vu_nm1, Vu_np1, numEq);
		}
		else
		{
			//===solve for acceleration
			matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_MGlo,
				cdp_FGlo, dp_A, &comm);
			MPI_Bcast(dp_A, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			////===updat displacement and store the results;
			matoperat.matMultiply(Va_n, dt* dt, Vu_np1); // a*t^2;
			matoperat.matAdd(Vu_np1, 2.0, Vu_n, Vu_np1);// 2*U_n+a*t^2;
			matoperat.matMinus(Vu_np1, Vu_nm1, Vu_np1);//2*U_n+a*t^2-U_(n-1);
		}
		//==store displacement
		storeDisplacementResult(o_dat, Vu_np1);
		//calculate stress==
		calGlobalNodeStresses(o_dat);
		if (o_dat.ci_failFlag == 0)
		{
			//===updat bond-state==
			failureCriterion_stress(o_dat);
		}
		else if (o_dat.ci_failFlag == 1)
		{
			failureCriterion_stretch(o_dat);
		}
		////==calculta reaction force if needed;

		//write results===
		if (n%saveFren==0||n== numTstep)
		{
			calLocalDamage(o_dat);
			if (ci_rank == 0)
			{
				printf("Writing results at step of %d ......\n", n);
				o_files.writeReslutsTOTAL_vtk(o_dat, n); //**set phi
			}
		}
	}

	// release memory;
	delete Vu_nm1, delete Vu_n, delete Vu_np1, delete Va_n; // no need to delete dp_A, if Va_n is deleted;
	Vu_nm1 = NULL, Vu_n = NULL, Vu_np1 = NULL, Va_n = NULL;
	delete[] cdp_F, delete[] cdp_FGlo, delete[] cdp_MGlo;
	cdp_F = NULL, cdp_FGlo = NULL, cdp_MGlo = NULL;
	if (cip_ia)
	{
		delete[] cip_ia, delete[] cip_ja;
		cip_ia = NULL, cip_ja = NULL;
	}
}

void pdsolve::pdfemDynamicNewmarkSolver_CSRformat(datModel& o_dat, fioFiles& o_files)
{
	//Don't need block data any more for dynamic solver;
	o_dat.deleteBLOCK();
	//================================================================
	o_dat.ci_solvFlag = 0;// remove later;
	if (o_dat.ci_solvFlag != 0)
	{
		printf("ERROR: This is not dynamic solver\n");
		printf("Please reset the solver flag\n");
		exit(0);
	}
	//============number of equations===================
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode;
	int numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = totNumNode * numDime - numPreDof;
	long long int numEq_long = numEq;
	//========================================================================
	//=============inital some variables=====
	/*Vu_nm1--Vector: displacement of (n-1)th step;
	 Vu_n--Vector: displacement of (n)th step;
	 Vu_np1--Vector: displacement of (n+1)th step;
	 Va_n--Vector: accerations of (n)th step;*/
	Vector* Vu_n, * Vu_np1, * Va_n, * Va_np1, * Vv_n, * Vv_np1;

	//=====mass==========
	if (o_dat.cb_lumpedMass)
	{
		assembleLumpedMass(o_dat, numEq);
	}
	else
	{
		//==========mass matrix CSR format======
		setCSRIndexes_gloMassMat(o_dat);
		assembleElemassMatFEM_CSRformat(o_dat);
		assembleMassMatPD_CSRformat(o_dat);
		MPI_Allreduce(cdp_M, cdp_MGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		delete[]cdp_M;
		cdp_M = NULL;
	}

	//===================
	//===============start dynamic solving==========
	double dt = o_dat.getTstep();
	if (ci_rank == 0)
	{
		printf("start to Newmark's Dynamic solving......\n");
		printf("Time step is %e\n", dt);
	}
}

void pdsolve::pdfemQuasiStaticAssembleEquasSys_CSRformat(datModel& o_dat, int numEq)
{
	o_dat.ci_solvFlag = 2;
	//initial K matrix
	for (long long int i = 0; i < cip_ia[numEq]; i++)
	{
		cdp_Ku[i] = 0;
		cdp_KuGlo[i] = 0;
	}
	//====internal force;
	calinternalForce_CSRformat(o_dat, numEq);
	//==set varied NBC (NBC, not prescribed displacemets); 
	setVaryNaturalBC(o_dat);
	//==external force;
	calExternalForce_CSRformat(o_dat);
	//assembling equations
	cb_InteralForce = false;//not for calculate internal force
	assembleInterWorkPD_CSRformat(o_dat);
	assemblePDBEwork_CSRformat(o_dat);
	assembleSEDbyFEM_CSRformat(o_dat);
	//==all reduce sum;
	MPI_Allreduce(cdp_Ku, cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void pdsolve::pdfemQuasiStaticSolver_CSRformat(datModel& o_dat, fioFiles& o_files)
{
	int saveFren = o_dat.getSaveFreq();
	int numTstep = o_dat.getnumTstep();
	//============number of equations===================
	int numDime = o_dat.ci_Numdimen;
	int numPreDof = 0;
	for (int i = 0; i < o_dat.getTotnumEssentialBCs(); i++)
	{
		numPreDof = numPreDof + o_dat.op_getEssenBC(i)->getNumNODE();
	}
	int totNumNode;
	int numEq;
	totNumNode = o_dat.getTotnumNode();
	numEq = totNumNode * numDime - numPreDof;
	long long int numEq_long = numEq;
	//=============inital some variables=====
	Vector* Vdu, * Vu_n;
	Vu_n = new Vector(numEq);
	Vdu = new Vector(numEq);
	int comm = MPI_Comm_c2f(MPI_COMM_WORLD);
	//===============start dynamic solving==========
	if (ci_rank == 0)
	{
		printf("start to quasi-static solving......\n");
	}
	//============get initial step displacement;
	//initial Vu_n========
	int eqIndex;
	for (int k = 0; k < totNumNode; k++)
	{
		for (int i = 0; i < numDime; i++)
		{
			eqIndex = o_dat.op_getNode(k)->op_getDof(i)->i_getEqInde();
			if (eqIndex != -1)
			{
				Vu_n->setCoeff(eqIndex, o_dat.op_getNode(k)->op_getDof(i)->d_getValue());
			}
		}
	}
	//===set global stiff K as CSR format;
	setCSRIndexes_gloStiffMat(o_dat);
	//===start loop=============;
	for (int n = 1; n < numTstep + 1; n++)
	{
		if (ci_rank == 0)
		{
			cout << "load step n= " << n << endl;
		}
		//assenble Euqation system;
		pdfemQuasiStaticAssembleEquasSys_CSRformat(o_dat, numEq);

		// solver du;
		matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_KuGlo,
			cdp_FGlo, Vdu->cdp_vecCoeff, &comm);

		/*store results========;*/
		MPI_Bcast(Vdu->cdp_vecCoeff, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		matoperat.matAdd(Vu_n, Vdu, Vu_n);
		storeDisplacementResult(o_dat, Vu_n);
		/*==stresses==========*/
		calGlobalNodeStresses(o_dat);
		//===updat bond-state==
		if (o_dat.ci_failFlag == 0)
		{
			failureCriterion_stress(o_dat);
		}
		else if (o_dat.ci_failFlag == 1)
		{
			failureCriterion_stretch(o_dat);
		}
		//==calculta reaction force if needed;

		//write results===
		if (n % saveFren == 0 || n == numTstep)
		{
			calLocalDamage(o_dat);
			if (ci_rank == 0)
			{
				printf("Writing results at step of %d ......\n", n);
				o_files.writeReslutsTOTAL_vtk(o_dat, n); //**set phi
			}
		}
	}
	//release memory;
	delete Vu_n, delete Vdu;
	Vu_n = NULL, Vdu = NULL;
	delete[] cip_ia, delete[] cip_ja, delete[] cdp_F, delete[] cdp_FGlo;
	delete[] cdp_Ku, delete[]cdp_KuGlo;
	cip_ia = NULL; cip_ja = NULL, cdp_F = NULL, cdp_FGlo = NULL;
	cdp_Ku = NULL, cdp_KuGlo = NULL;
}

double pdsolve::failureCriterion_stretch(datModel& o_dat)
{
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	pdFamily* temP_fami;
	vector<int>mp_fami(numFami);
	vector<int>mp_fami_glo(numFami);
	int numNodeOfFami, NID_k, NID_m;
	double xk[3] = { 0 }, xm[3] = { 0 }, uk[3] = { 0 }, um[3] = { 0 }, xi[3] = { 0 }, eta[3] = { 0 };
	double mag_ref,mag, s,s_max=0, sc;//get latter
	double  KIC, G0, kappa, E, nu, delta_k, delta_m, delta, mu_km;
	double fac = o_dat.op_getGeomP()->getFactor();
	E = o_dat.op_getmaterial()->getE();
	nu = o_dat.op_getmaterial()->getnu();
	KIC = o_dat.op_getmaterial()->getKIc();
	G0 = KIC * KIC / E;
	kappa = E / (3.0 * (1 - 2 * nu));
	for (int famkk = startP; famkk < endP; famkk++)
	{
		temP_fami = o_dat.op_getFami(famkk);
		if (temP_fami->cb_allowFail)
		{
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			o_dat.op_getNode(NID_k - 1)->getcoor(xk);
			delta_k = temP_fami->gethorizon();
			for (int i = 0; i < numDime; i++)
			{
				uk[i] = o_dat.op_getNode(NID_k - 1)->op_getDof(i)->d_getValue();
			}
			for (int m = 1; m < numNodeOfFami; m++)
			{
				mu_km = temP_fami->getbondstate(m);
				if (mu_km == 1)
				{
					NID_m = temP_fami->getNodeID(m);
					o_dat.op_getNode(NID_m - 1)->getcoor(xm);
					for (int i = 0; i < numDime; i++)
					{
						um[i] = o_dat.op_getNode(NID_m - 1)->op_getDof(i)->d_getValue();
						xi[i] = xm[i] - xk[i];
						eta[i] = um[i] - uk[i];
					}
					mag_ref = sqrt(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2]);
					mag = sqrt((xi[0] + eta[0]) * (xi[0] + eta[0]) + (xi[1] + eta[1]) * (xi[1] + eta[1]) +
						(xi[2] + eta[2]) * (xi[2] + eta[2]));
					s = mag / mag_ref - 1.0;

					if (o_dat.ci_solvFlag != 1)//non-static solver
					{
						delta_m = fac * pow((o_dat.op_getNode(NID_m - 1)->getvolume()), 1.0 / numDime);
						delta = 0.5 * (delta_k + delta_m);
						sc = sqrt(5 * G0 / (9 * kappa * delta));
						if (s >= sc)
						{
							temP_fami->setbondstate(m, 0);//only non-static solver set mu=0;
							mp_fami[famkk] = 1;
						}

					}
					else if (o_dat.ci_solvFlag == 1)//static solver
					{
						if (s > s_max)
						{
							s_max = s;
						}
					}
				}

			}
		}
		
	}
	//==== sent and receive data;
	double maxS = 0;
	if (o_dat.ci_solvFlag !=1&&numFami>0)//non -static solver
	{
		MPI_Allreduce(&(mp_fami[0]), &(mp_fami_glo[0]), numFami, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		int LocSP, LocEP;
		for (int rank = 0; rank < ci_numProce; rank++)
		{
			LocSP = rank * numFami / ci_numProce;
			LocEP = (rank + 1) * numFami / ci_numProce;
			for (int famk = LocSP; famk < LocEP; famk++)
			{
				if (mp_fami_glo[famk] == 1)
				{
					numNodeOfFami = o_dat.op_getFami(famk)->getNumNode();
					MPI_Bcast(&(o_dat.op_getFami(famk)->cip_bondState[0]), numNodeOfFami,
						MPI_INT, rank, MPI_COMM_WORLD);
				}
			}
		}
	}
	else if (o_dat.ci_solvFlag==1)//static solver
	{
		MPI_Allreduce(&s_max, &maxS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	return maxS;
}

double pdsolve::failureCriterion_stress(datModel& o_dat)
{
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	//startP = 0;
	//endP = numFami;

	pdFamily* temP_fami;
	vector<int>mp_fami(numFami);
	vector<int>mp_fami_glo(numFami);
	int numNodeOfFami, NID_k, NID_m;
	double sig[6], sigma_k[6], sigma_m[6], sigma1, sigult, sigma1_max = 0, mu_km;
	sigult = o_dat.op_getmaterial()->getSigult();
	//Matrix* sigma_kg = new Matrix(3, 3);
	//Matrix* simga_priaxis = new Matrix(3, 3);// dummy here.
	//Vector* sigma_pri = new Vector(3);
	for (int famkk = startP; famkk < endP; famkk++)
	{
		temP_fami = o_dat.op_getFami(famkk);
		if (temP_fami->cb_allowFail)
		{
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			o_dat.op_getNode(NID_k - 1)->getStress(sigma_k);
			for (int m = 1; m < numNodeOfFami; m++)
			{
				mu_km = temP_fami->getbondstate(m);
				if (mu_km == 1)
				{
					NID_m = temP_fami->getNodeID(m);
					o_dat.op_getNode(NID_m - 1)->getStress(sigma_m);
					// symmetry matrix, store by upper triangle
				/*	for (int ii = 0; ii < 3; ii++)
					{
						sigma_kg->setCoeff(ii, ii, 0.5 * (sigma_k[ii] + sigma_m[ii]));
					}
					sigma_kg->setCoeff(0, 1, 0.5 * (sigma_k[3] + sigma_m[3]));
					sigma_kg->setCoeff(0, 2, 0.5 * (sigma_k[5] + sigma_m[5]));
					sigma_kg->setCoeff(1, 2, 0.5 * (sigma_k[4] + sigma_m[4]));
					matoperat.dSymeEigenV('N', sigma_kg, sigma_pri, simga_priaxis);
					sigma1 = *max_element(sigma_pri->cdp_vecCoeff, sigma_pri->cdp_vecCoeff + 3);*/
					for (int ii = 0; ii < 6; ii++)
					{
						sig[ii] = 0.5 * (sigma_k[ii] + sigma_m[ii]);
					}
					sigma1 = sqrt(0.5 * ((sig[0] - sig[1]) * (sig[0] - sig[1]) +
						(sig[0] - sig[2]) * (sig[0] - sig[2]) + (sig[1] - sig[2]) * (sig[1] - sig[2])) +
						3 * (sig[3] * sig[3] + sig[4] * sig[4] + sig[5] * sig[5]));
					if (o_dat.ci_solvFlag !=1)//non-static solver
					{
						if (sigma1 > sigult)
						{
							temP_fami->setbondstate(m, 0);//only dynamic solver set mu=0;
							mp_fami[famkk] = 1;
						}
					}
					else if (o_dat.ci_solvFlag == 1)//static solver
					{
						if (sigma1 > sigma1_max)
						{
							sigma1_max = sigma1;
						}
					}
				}

			}
		}
		
	}

	//==== sent and receive data;
	double maxSig = 0;
	if (o_dat.ci_solvFlag!=1&&numFami>0)//non-static solver
	{
		MPI_Allreduce(&(mp_fami[0]), &(mp_fami_glo[0]), numFami, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		int LocSP, LocEP;
		for (int rank = 0; rank < ci_numProce; rank++)
		{
			LocSP = rank * numFami / ci_numProce;
			LocEP = (rank + 1) * numFami / ci_numProce;
			for (int famk = LocSP; famk < LocEP; famk++)
			{
				if (mp_fami_glo[famk] == 1)
				{
					numNodeOfFami = o_dat.op_getFami(famk)->getNumNode();
					MPI_Bcast(&(o_dat.op_getFami(famk)->cip_bondState[0]), numNodeOfFami,
						MPI_INT, rank, MPI_COMM_WORLD);
				}
			}
		}
	}
	else if (o_dat.ci_solvFlag==1)//static solver;
	{
		MPI_Allreduce(&sigma1_max, &maxSig, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	/*delete sigma_kg, delete simga_priaxis, delete sigma_pri;
	sigma_kg = NULL, simga_priaxis = NULL, sigma_pri = NULL;*/
	return maxSig;
	//return sigma1_max;
}

double pdsolve::failureCriterion_stressCut(datModel& o_dat)
{
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	pdFamily* temP_fami;
	vector<int>mp_fami(numFami);
	vector<int>mp_fami_glo(numFami);
	int numNodeOfFami, NID_k, NID_m;
	double sig[6], sigma_k[6], sigma_m[6], sigma1, sigult, sigma1_max = 0, mu_km;
	sigult = o_dat.op_getmaterial()->getSigult();
	Matrix* sigma_kg = new Matrix(3, 3);
	Matrix* simga_priaxis = new Matrix(3, 3);// dummy here.
	Vector* sigma_pri = new Vector(3);
	int maxIDX;
	double eigV[3];
	for (int famkk = startP; famkk < endP; famkk++)
	{
		temP_fami = o_dat.op_getFami(famkk);
		if (temP_fami->cb_allowFail)
		{
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			o_dat.op_getNode(NID_k - 1)->getStress(sigma_k);
			for (int m = 1; m < numNodeOfFami; m++)
			{
				mu_km = temP_fami->getbondstate(m);
				o_dat.op_getNode(NID_m - 1)->getStress(sigma_m);
				// symmetry matrix, store by upper triangle
				for (int ii = 0; ii < 3; ii++)
				{
					sigma_kg->setCoeff(ii, ii, 0.5 * (sigma_k[ii] + sigma_m[ii]));
				}
				sigma_kg->setCoeff(0, 1, 0.5 * (sigma_k[3] + sigma_m[3]));
				sigma_kg->setCoeff(0, 2, 0.5 * (sigma_k[5] + sigma_m[5]));
				sigma_kg->setCoeff(1, 2, 0.5 * (sigma_k[4] + sigma_m[4]));
				matoperat.dSymeEigenV('V', sigma_kg, sigma_pri, simga_priaxis);
				sigma1 = *max_element(sigma_pri->cdp_vecCoeff, sigma_pri->cdp_vecCoeff + 3);
				maxIDX = max_element(sigma_pri->cdp_vecCoeff, sigma_pri->cdp_vecCoeff + 3)- sigma_pri->cdp_vecCoeff;
				for (int ii = 0; ii < 3; ii++)
				{
					eigV[ii] = simga_priaxis->d_getCoeff(ii, maxIDX);
				}
			}
		}
	}

	return 0.0;
}

void pdsolve::calLocalDamage(datModel& o_dat)
{
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	pdFamily* temP_fami;
	int numNodeOfFami;
	double tempDam;
	vector<double>locaDama(numFami);
	vector<double>glo_locaDama(numFami);
	for (int famkk = startP; famkk < endP; famkk++)
	{
		tempDam = 0;
		temP_fami = o_dat.op_getFami(famkk);
		numNodeOfFami = temP_fami->getNumNode();
		for (int m = 1; m < numNodeOfFami; m++)
		{
			tempDam = tempDam + temP_fami->getbondstate(m);
		}
		locaDama[famkk] = 1.0 - tempDam / (numNodeOfFami-1.0);
	}
	if (numFami>0)
	{
		MPI_Reduce(&(locaDama[0]), &(glo_locaDama[0]), numFami, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	if (ci_rank==0)
	{
		int NID_k;
		for (int famk = 0; famk < numFami; famk++)
		{
			NID_k = o_dat.op_getFami(famk)->getNodeID(0);
			o_dat.op_getNode(NID_k - 1)->setLocalDamage(glo_locaDama[famk]);
		}
	}
}

void pdsolve::calBondstate(double x1[], double x2[], datModel & o_dat)
{
	//calculate bond state;
	//assuming crack between point x1 and x2;
	int numFam = o_dat.getTotnumFami();
	int NID_k, NID_m, numNodeOfFam;
	double x_k[2], x_m[2];
	int bondstate;
	
	for (int k = 0; k < numFam; k++)
	{
		numNodeOfFam = o_dat.op_getFami(k)->getNumNode();
		NID_k = o_dat.op_getFami(k)->getNodeID(0);
		
		o_dat.op_getNode(NID_k - 1)->getcoor(x_k);
		for (int m = 1; m < numNodeOfFam; m++)
		{
			bondstate = o_dat.op_getFami(k)->getbondstate(m);
			if (bondstate==1)
			{
				NID_m= o_dat.op_getFami(k)->getNodeID(m);
				o_dat.op_getNode(NID_m - 1)->getcoor(x_m);
				if (intersection(x_k, x_m, x1, x2))
				{
					o_dat.op_getFami(k)->setbondstate(m, 0);
				}
			}
		}

	}
}

bool pdsolve::intersection(double L1X1[], double L1X2[], double L2X1[], double L2X2[])
{
	//exclude test
	if ((L1X1[0] > L1X2[0] ? L1X1[0] : L1X2[0]) < (L2X1[0] < L2X2[0] ? L2X1[0] : L2X2[0]) ||
		(L1X1[1] > L1X2[1] ? L1X1[1] : L1X2[1]) < (L2X1[1] < L2X2[1] ? L2X1[1] : L2X2[1]) ||
		(L2X1[0] > L2X2[0] ? L2X1[0] : L2X2[0]) < (L1X1[0] < L1X2[0] ? L1X1[0] : L1X2[0]) ||
		(L2X1[1] > L2X2[1] ? L2X1[1] : L2X2[1]) < (L1X1[1] < L1X2[1] ? L1X1[1] : L1X2[1]))
	{
		return false;
	}
	//pass test

	if ((((L1X1[0] - L2X1[0])*(L2X2[1] - L2X1[1]) - (L1X1[1] - L2X1[1])*(L2X2[0] - L2X1[0]))*
		((L1X2[0] - L2X1[0])*(L2X2[1] - L2X1[1]) - (L1X2[1] - L2X1[1])*(L2X2[0] - L2X1[0]))) > 0 ||
		(((L2X1[0] - L1X1[0])*(L1X2[1] - L1X1[1]) - (L2X1[1] - L1X1[1])*(L1X2[0] - L1X1[0]))*
		((L2X2[0] - L1X1[0])*(L1X2[1] - L1X1[1]) - (L2X2[1] - L1X1[1])*(L1X2[0] - L1X1[0]))) > 0)
	{
		return false;
	}
	return true;
}

void pdsolve::calPDNodeStresses(datModel & o_dat, int* count)
{
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	pdFamily* temP_fami;
	if (numDime==2)
	{
		Matrix* C;
		Vector* epsilon, * sigma, * uk;
		epsilon = new Vector(3);
		sigma = new Vector(3);
		int numNodeOfFam, Nid_m, Nid_k;
		double tempu;
		for (int famkk = startP; famkk < endP; famkk++)
		{
			temP_fami = o_dat.op_getFami(famkk);
			numNodeOfFam = temP_fami->getNumNode();
			Nid_k = temP_fami->getNodeID(0);
			(count[Nid_k - 1])++;
			C = new Matrix(3, 2 * numNodeOfFam);
			uk = new Vector(2 * numNodeOfFam);
			matC2D(C, temP_fami, o_dat);
			for (int m = 0; m < numNodeOfFam; m++)
			{
				Nid_m = temP_fami->getNodeID(m);
				for (int i = 0; i < 2; i++)
				{
					tempu = o_dat.op_getNode(Nid_m - 1)->op_getDof(i)->d_getValue();
					uk->setCoeff(2 * m + i, tempu);
				}
			}
			matoperat.matMultiply(C, uk, epsilon);
			matoperat.matMultiply(cop_D, epsilon, sigma);
			for (int ii = 0; ii < 2; ii++)
			{
				o_dat.op_getNode(Nid_k - 1)->setStress(ii, sigma->d_getCoeff(ii));
			}
			o_dat.op_getNode(Nid_k - 1)->setStress(3, sigma->d_getCoeff(2));
			delete C, delete uk;
			C = NULL; uk = NULL;
		}
		delete sigma, delete epsilon;
		sigma = NULL; epsilon = NULL;
	}
	else if (numDime==3)
	{
		Matrix* C;
		Vector* epsilon, * sigma, * uk;
		epsilon = new Vector(6);
		sigma = new Vector(6);
		int numNodeOfFam, Nid_m, Nid_k;
		double tempu;
		for (int famkk = startP; famkk < endP; famkk++)
		{
			temP_fami = o_dat.op_getFami(famkk);
			numNodeOfFam = temP_fami->getNumNode();
			Nid_k = temP_fami->getNodeID(0);
			(count[Nid_k - 1]) = count[Nid_k - 1] + 1;
			C = new Matrix(6, numDime * numNodeOfFam);
			uk = new Vector(numDime * numNodeOfFam);
			matC3D(C, temP_fami, o_dat);
			for (int m = 0; m < numNodeOfFam; m++)
			{
				Nid_m = temP_fami->getNodeID(m);
				for (int i = 0; i < numDime; i++)
				{
					tempu = o_dat.op_getNode(Nid_m - 1)->op_getDof(i)->d_getValue();
					uk->setCoeff(numDime * m + i, tempu);
				}
			}
			matoperat.matMultiply(C, uk, epsilon);
			matoperat.matMultiply(cop_D, epsilon, sigma);
			for (int ii = 0; ii < 6; ii++)
			{
				o_dat.op_getNode(Nid_k - 1)->setStress(ii, sigma->d_getCoeff(ii));
			}
			delete C, delete uk;
			C = NULL; uk = NULL;
		}
		delete sigma, delete epsilon;
		sigma = NULL; epsilon = NULL;
	}
}

void pdsolve::calFEMNodeStresses_EXP(datModel& o_dat, int* count)
{
	int totNumFE, startP, endP;
	totNumFE = o_dat.civ_feIDX.size();
	startP = ci_rank * totNumFE / ci_numProce;
	endP = (ci_rank + 1) * totNumFE / ci_numProce;
	////====extrapolation method===
	int numDime = o_dat.ci_Numdimen;
	if (numDime == 3)
	{
		//=====initial matrix L;
		Matrix* L;
		L = new Matrix(8, 8);
		// represent Gauss points coordinates in local;
		double p_I[8] = { -1.0,1.,1.,-1.,-1.,1.,1.,-1. };
		double q_I[8] = { -1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0 };
		double r_I[8] = { -1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0 };
		double p[8], q[8], r[8], temp;
		for (int i = 0; i < 8; i++)
		{
			//represents node coordinate in local 
			p[i] = sqrt(3.0) * p_I[i];
			q[i] = sqrt(3.0) * q_I[i];
			r[i] = sqrt(3.0) * r_I[i];
		}
		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 8; j++)
			{
				temp = 1.0 / 8 * (1 + p_I[j] * p[i]) * (1 + q_I[j] * q[i]) * (1 + r_I[j] * r[i]);
				L->setCoeff(i, j, temp);
			}
		}
		//===get stresses of each fem elements
		int algoType, numNodeEle, * conNID, ele;
		Vector* Ue, * Nsig[6];
		for (int i = 0; i < 6; i++)
		{
			Nsig[i] = new Vector(8);
		}
		double(*xN)[3], tempu;
		for (int fe = startP; fe < endP; fe++)
		{
			ele = o_dat.civ_feIDX[fe];
			//algoType = o_dat.op_getEles(k)->getAlgoType();
			//if (algoType == 2)
			{
				numNodeEle = o_dat.op_getEles(ele)->getNumNodes();
				Ue = new Vector(numDime * numNodeEle);
				xN = new double[numNodeEle][3];
				conNID = new int[numNodeEle];
				o_dat.op_getEles(ele)->getConNid(conNID);
				for (int nd = 0; nd < numNodeEle; nd++)
				{
					o_dat.op_getNode(conNID[nd] - 1)->getcoor(xN[nd]);
					for (int j = 0; j < numDime; j++)
					{
						tempu = o_dat.op_getNode(conNID[nd] - 1)->op_getDof(j)->d_getValue();
						Ue->setCoeff(numDime * nd + j, tempu);
					}
					count[conNID[nd] - 1] = count[conNID[nd] - 1] + 1;
				}
				o_dat.op_getEles(ele)->eleFitStresses(1, Nsig, cop_D, L, Ue, xN);

				for (int nd = 0; nd < numNodeEle; nd++)//caution: only calculat the corner values
				{
					for (int si = 0; si < 6; si++)
					{
						o_dat.op_getNode(conNID[nd] - 1)->addStress(si, Nsig[si]->d_getCoeff(nd));
					}
				}
				delete[] xN, delete[] conNID;
				xN = NULL, conNID = NULL;
				delete Ue;
				Ue = NULL;
			}
		}
		for (int i = 0; i < 6; i++)
		{
			delete Nsig[i];
			Nsig[i] = NULL;
		}
		delete L;
		L = NULL;
	}
	else if (numDime == 2)
	{
		//=====initial matrix L;
		Matrix* L;
		L = new Matrix(4, 4);
		double a, b, c;
		a = 1 + 0.5 * sqrt(3.0);
		b = -0.5;
		c = 1 - 0.5 * sqrt(3.0);
		double row0[4] = { a,b,c,b };
		for (int i = 0; i < 4; i++)
		{
			L->setCoeff(0, i, row0[i]);
		}
		double row1[3] = { a,b,c };
		for (int i = 0; i < 3; i++)
		{
			L->setCoeff(1, i + 1, row1[i]);
		}
		L->setCoeff(2, 2, a);
		L->setCoeff(2, 3, b);
		L->setCoeff(3, 3, a);
		for (int i = 1; i < 4; i++)
		{
			for (int j = 0; j < i; j++)
			{
				L->setCoeff(i, j, L->d_getCoeff(j, i));
			}
		}
		//===get stresses of each fem elements
		int algoType, numNodeEle, * conNID, ele;
		Vector* Ue, * Nsig[3];
		for (int i = 0; i < 3; i++)
		{
			Nsig[i] = new Vector(4);
		}
		double(*xN)[3], tempu;
		for (int fe = startP; fe < endP; fe++)
		{
			ele = o_dat.civ_feIDX[fe];
			//algoType = o_dat.op_getEles(k)->getAlgoType();
			//if (algoType == 2)
			{
				numNodeEle = o_dat.op_getEles(ele)->getNumNodes();
				Ue = new Vector(numDime * numNodeEle);
				xN = new double[numNodeEle][3];
				conNID = new int[numNodeEle];
				o_dat.op_getEles(ele)->getConNid(conNID);
				for (int nd = 0; nd < numNodeEle; nd++)
				{
					o_dat.op_getNode(conNID[nd] - 1)->getcoor(xN[nd]);
					for (int j = 0; j < numDime; j++)
					{
						tempu = o_dat.op_getNode(conNID[nd] - 1)->op_getDof(j)->d_getValue();
						Ue->setCoeff(numDime * nd + j, tempu);
					}
					count[conNID[nd] - 1] = count[conNID[nd] - 1] + 1;
				}
				o_dat.op_getEles(ele)->eleFitStresses(1, Nsig, cop_D, L, Ue, xN);
				for (int nd = 0; nd < numNodeEle; nd++)// Caution: only calculat the corner values
				{
					for (int si = 0; si < 2; si++)
					{
						o_dat.op_getNode(conNID[nd] - 1)->addStress(si, Nsig[si]->d_getCoeff(nd));
					}
					o_dat.op_getNode(conNID[nd] - 1)->addStress(3, Nsig[2]->d_getCoeff(nd));
				}
				delete[] xN, delete[] conNID;
				xN = NULL, conNID = NULL;
				delete Ue;
				Ue = NULL;
			}
		}
		for (int i = 0; i < 3; i++)
		{
			delete Nsig[i];
			Nsig[i] = NULL;
		}
		delete L;
		L = NULL;
	}

}

void pdsolve::calFEMNodeStresses_LSM(datModel& o_dat, int* count)
{
	int totNumFE, startP, endP;
	totNumFE = o_dat.civ_feIDX.size();
	startP = ci_rank * totNumFE / ci_numProce;
	endP = (ci_rank + 1) * totNumFE / ci_numProce;
	//=============LSM method==========================;
	int numDime = o_dat.ci_Numdimen;
	if (numDime==2)
	{
		int nG;
		Matrix* fitA;
		fitA = new Matrix(3, 3);
		fitA->zero();
		//==get fitA;
		nG = o_globGP.i_getNumPts();
		double* pp, * qq, p, q;
		pp = new double[nG * nG];
		qq = new double[nG * nG];
		for (int np = 0; np < nG; np++)
		{
			p = o_globGP.d_getGaussPt(np);
			for (int nq = 0; nq < nG; nq++)
			{
				q = o_globGP.d_getGaussPt(nq);
				pp[np * nG + nq] = p;
				qq[np * nG + nq] = q;
			}
		}
		for (int i = 0; i < nG * nG; i++)
		{
			fitA->addCoeff(0, 0, 1.0);
			fitA->addCoeff(0, 1, pp[i]);
			fitA->addCoeff(0, 2, qq[i]);
			fitA->addCoeff(1, 0, pp[i]);
			fitA->addCoeff(1, 1, pp[i] * pp[i]);
			fitA->addCoeff(1, 2, pp[i] * qq[i]);
			fitA->addCoeff(2, 0, qq[i]);
			fitA->addCoeff(2, 1, pp[i] * qq[i]);
			fitA->addCoeff(2, 2, qq[i] * qq[i]);
		}
		delete[] pp; delete[] qq;
		pp = NULL; qq = NULL;
		//===get stresses of each fem elements
		int algoType, numNodeEle, * conNID, ele;
		Vector* Ue, * Nsig[3];
		double (*xN)[3], tempu;
		for (int fe = startP; fe < endP; fe++)
		{
			ele = o_dat.civ_feIDX[fe];
			/*algoType = o_dat.op_getEles(k)->getAlgoType();
			if (algoType == 2)*/
			/*{*/
				numNodeEle = o_dat.op_getEles(ele)->getNumNodes();
				Ue = new Vector(numDime * numNodeEle);
				xN = new double[numNodeEle][3];
				for (int ii = 0; ii < 3; ii++)
				{
					Nsig[ii] = new Vector(numNodeEle);
				}
				conNID = new int[numNodeEle];
				o_dat.op_getEles(ele)->getConNid(conNID);
				for (int nd = 0; nd < numNodeEle; nd++)
				{
					o_dat.op_getNode(conNID[nd] - 1)->getcoor(xN[nd]);
					for (int j = 0; j < numDime; j++)
					{
						tempu = o_dat.op_getNode(conNID[nd] - 1)->op_getDof(j)->d_getValue();
						Ue->setCoeff(numDime * nd + j, tempu);
					}
					count[conNID[nd] - 1] = count[conNID[nd] - 1] + 1;
				}
				o_dat.op_getEles(ele)->eleFitStresses(2, Nsig, cop_D, fitA, Ue, xN);
				for (int nd = 0; nd < numNodeEle; nd++)//only calculat the corner values
				{
					o_dat.op_getNode(conNID[nd] - 1)->addStress(0, Nsig[0]->d_getCoeff(nd));
					o_dat.op_getNode(conNID[nd] - 1)->addStress(1, Nsig[1]->d_getCoeff(nd));
					o_dat.op_getNode(conNID[nd] - 1)->addStress(3, Nsig[2]->d_getCoeff(nd));//sig_xy, caution
				}
				delete[] xN, delete[] conNID;
				xN = NULL, conNID = NULL;
				delete Ue;
				Ue = NULL;
				for (int ii = 0; ii < 3; ii++)
				{
					delete Nsig[ii];
					Nsig[ii] = NULL;
				}
			//}
		}
		delete fitA;
		fitA = NULL;
	}
	else if (numDime==3)
	{
		int nG;
		Matrix* fitA;
		fitA = new Matrix(4, 4);
		fitA->zero();
		//==get fitA;
		nG = o_globGP.i_getNumPts();
		double* pp, * qq, * rr, p, q, r, temp;
		pp = new double[nG * nG * nG];
		qq = new double[nG * nG * nG];
		rr = new double[nG * nG * nG];
		for (int np = 0; np < nG; np++)
		{
			p = o_globGP.d_getGaussPt(np);
			for (int nq = 0; nq < nG; nq++)
			{
				q = o_globGP.d_getGaussPt(nq);
				for (int nr = 0; nr < nG; nr++)
				{
					r= o_globGP.d_getGaussPt(nr);
					pp[np * nG * nG + nq * nG + nr] = p;
					qq[np * nG * nG + nq * nG + nr] = q;
					rr[np * nG * nG + nq * nG + nr] = r;
				}
			}
		}
		for (int i = 0; i < nG * nG * nG; i++)
		{
			fitA->addCoeff(0, 0, 1.0);
			fitA->addCoeff(0, 1, pp[i]);
			fitA->addCoeff(0, 2, qq[i]);
			fitA->addCoeff(0, 3, rr[i]);
			fitA->addCoeff(1, 1, pp[i] * pp[i]);
			fitA->addCoeff(1, 2, pp[i] * qq[i]);
			fitA->addCoeff(1, 3, pp[i] * rr[i]);
			fitA->addCoeff(2, 2, qq[i] * qq[i]);
			fitA->addCoeff(2, 3, qq[i] * rr[i]);
			fitA->addCoeff(3, 3, rr[i] * rr[i]);
		}
		temp = fitA->d_getCoeff(0, 1);
		fitA->setCoeff(1, 0, temp);
		temp = fitA->d_getCoeff(0, 2);
		fitA->setCoeff(2, 0, temp);
		temp = fitA->d_getCoeff(1, 2);
		fitA->setCoeff(2, 1, temp);
		temp = fitA->d_getCoeff(0, 3);
		fitA->setCoeff(3, 0, temp);
		temp = fitA->d_getCoeff(1, 3);
		fitA->setCoeff(3, 1, temp);
		temp = fitA->d_getCoeff(2, 3);
		fitA->setCoeff(3, 2, temp);
		delete[] pp, delete[] qq, delete[] rr;
		pp = NULL; qq = NULL; rr = NULL;
		//===get stresses of each fem elements
		int algoType, numNodeEle, * conNID, ele;
		Vector* Ue, * Nsig[6];
		double (*xN)[3], tempu;
		for (int fe = startP; fe < endP; fe++)
		{
			ele = o_dat.civ_feIDX[fe];
			/*algoType = o_dat.op_getEles(k)->getAlgoType();
			if (algoType == 2)*/
			{
				numNodeEle = o_dat.op_getEles(ele)->getNumNodes();
				Ue = new Vector(numDime * numNodeEle);
				for (int ii = 0; ii < 6; ii++)
				{
					Nsig[ii] = new Vector(numNodeEle);
				}
				xN = new double[numNodeEle][3];
				conNID = new int[numNodeEle];
				o_dat.op_getEles(ele)->getConNid(conNID);
				for (int nd = 0; nd < numNodeEle; nd++)
				{
					o_dat.op_getNode(conNID[nd] - 1)->getcoor(xN[nd]);
					for (int j = 0; j < numDime; j++)
					{
						tempu = o_dat.op_getNode(conNID[nd] - 1)->op_getDof(j)->d_getValue();
						Ue->setCoeff(numDime * nd + j, tempu);
					}
					count[conNID[nd] - 1] = count[conNID[nd] - 1] + 1;
				}
				o_dat.op_getEles(ele)->eleFitStresses(2, Nsig, cop_D, fitA, Ue, xN);
				for (int nd = 0; nd < numNodeEle; nd++)
				{
					for (int si = 0; si < 6; si++)
					{
						o_dat.op_getNode(conNID[nd] - 1)->addStress(si, Nsig[si]->d_getCoeff(nd));
					}
				}
				delete[] xN, delete[] conNID;
				xN = NULL, conNID = NULL;
				delete Ue;
				Ue = NULL;
				for (int ii = 0; ii < 6; ii++)
				{
					delete Nsig[ii];
					Nsig[ii] = NULL;
				}
			}
		}
		delete fitA;
		fitA = NULL;
	}



}

void pdsolve::calGlobalNodeStresses(datModel & o_dat)
{
	//initial stresses and the count;
	int numDime = o_dat.ci_Numdimen;
	int totNumNODE = o_dat.getTotnumNode();
	int *count = new int[totNumNODE];
	int* gloCount = new int[totNumNODE];
	double* glosig = new double[6 * totNumNODE];
	for (int i = 0; i < totNumNODE; i++)
	{
		gloCount[i] = 0;
		count[i] = 0;
		/*if (o_dat.op_getNode(i)->getNodeType())
		{
			count[i] = 1;
		}
		else
		{
			count[i] = 0;
		}*/
		for (int j = 0; j < 6; j++)
		{
			o_dat.op_getNode(i)->setStress(j, 0);
			glosig[i * 6 + j] = 0;
		}
	}
	//=== cal PD node stress by PD algorithem;
	calPDNodeStresses(o_dat,count); // be caution, always do PD node first;
	/*=====================================================================
	=======================================================================
	==================FEM nodal stress=====================================
	=======================================================================
	=======================================================================*/
	/*====extrapolation method===*/
	calFEMNodeStresses_EXP(o_dat, count);
	/*====LSM method===*/
	//calFEMNodeStresses_LSM(o_dat, count);

	//=============
	MPI_Allreduce(count, gloCount, totNumNODE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(o_dat.cop_datLev2->cdp_sigma, glosig, 6 * totNumNODE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce(count, gloCount, totNumNODE, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(o_dat.cop_datLev2->cdp_sigma, glosig, 6 * totNumNODE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	double* tempSig = o_dat.cop_datLev2->cdp_sigma;
	o_dat.cop_datLev2->cdp_sigma = glosig;
	delete[] tempSig; tempSig = NULL; /*Caution here;*/
	/*====get average stress=============;*/
	//if (ci_rank==0)
	{
		for (int i = 0; i < totNumNODE; i++)
		{
			o_dat.op_getNode(i)->calAverageStress(gloCount[i]);//parallel later;
		}
	}
	delete[] count, delete[] gloCount;
	count = NULL, gloCount = NULL;
}

void pdsolve::postProcessing(datModel & o_dat, ofstream & test)
{
	//calGlobalNodeStresses(o_dat, test);
}
