/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdsolve.h"
#include<unordered_set>
#include<algorithm>
extern calMatrixOperations matoperat;
extern pdGaussPt o_globGP;
const double PI = acos(-1.);
pdsolve::pdsolve(datModel & o_dat,int rank,int numProce)
{
	ci_rank = rank;
	ci_numProce = numProce;
	cd_blockFac = 2.5;
	//===initial the pointers============
	cop_M = NULL;
	cop_Ku = NULL; cop_F = NULL;
	cip_ia = NULL, cip_ja = NULL, cdp_F = NULL, cdp_Ku = NULL, cdp_Ug = NULL, cdp_M = NULL;
	cdp_FGlo = NULL, cdp_KuGlo = NULL;
	//initial parameters for Newmark's method
	cd_gamma = o_dat.cd_gamma;
	cd_beta = o_dat.cd_beta;
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
	int temp = 0;
	MPI_Bcast(&temp, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
	/*	if (ci_rank==0)
		{
			printf("finished findDomainDimen......\n");
		}*/
		calVolumeOfNode(o_dat);//mpi
		/*if (ci_rank == 0)
		{
			printf("finished calVolumeOfNode......\n");
		}*/
		setDeltaMaxMin(o_dat);//mpi
		/*if (ci_rank == 0)
		{
			printf("finished setDeltaMaxMin......\n");
		}*/
		setBlockAndFami(o_dat);//mpi
		/*if (ci_rank == 0)
		{
			printf("finished setBlockAndFami......\n");
		}*/
		initialBondState(o_dat);
		/*if (ci_rank == 0)
		{
			printf("finished initialBondState......\n");
		}*/
		setNoFailRegion(o_dat);
		/*if (ci_rank == 0)
		{
			printf("finished setNoFailRegion......\n");
		}*/
	}
	
	if (ci_rank == 0)
	{
		printf("Finished setting data model.\n");
	}
}

void pdsolve::findDomainDimen(datModel& o_dat)
{
	int numD = o_dat.ci_Numdimen;
	int numNode, startP, endP;
	numNode = o_dat.getTotnumNode();
	startP = ci_rank * numNode / ci_numProce;
	endP = (ci_rank + 1) * numNode / ci_numProce;

	double lbc[3] = { 0 }, rtc[3]={ 0 }, d_x[3] = { 0 }, lbcGlo[3] = { 0 }, rtcGlo[3] = { 0 };
	o_dat.op_getNode(startP)->getcoor(lbc);
	o_dat.op_getNode(startP)->getcoor(rtc);
	o_dat.op_getNode(startP)->getcoor(lbcGlo);
	o_dat.op_getNode(startP)->getcoor(rtcGlo);
	for (int i = startP+ 1; i < endP; i++)
	{
		o_dat.op_getNode(i)->getcoor(d_x);
		for (int j = 0; j < numD; j++)
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
		for (int i = 0; i < numD; i++)
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
	//for PD element
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
	//for FE
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
				/*if (o_dat.op_getNode(k)->op_getDof(i)->i_getEqInde()==-1)
				{
					printf("ERROR: The DOF of node %d is repeatly defined!\n", k + 1);
					exit(0);
				}*/
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
	if (o_dat.cb_FENSF==true)
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
	MPI_Allreduce(loc_V, glo_V, numNode, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
	blockSize = cd_blockFac * fac * maxDelta;
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
	if (o_dat.cb_FENSF==false)
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
	else
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
	//==allocate element into block; only for 2D quasi-static solver;
	int AlgoType, * conNID, numNele;
	double eleCen[3], (*xN)[3], *N;
	if (o_dat.ci_solvFlag==2&&o_dat.ci_Numdimen==2)
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
				//N = new double[numNele];
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
				delete[]conNID, delete[] xN;
				conNID = NULL, xN = NULL;
			}
		}
	}

	//===============set families===============================
	/*===parallize===*/
	int totNumPDN, startP, endP, nd;
	totNumPDN = o_dat.civ_pdNodeIDX.size();
	startP = ci_rank * totNumPDN / ci_numProce;
	endP = (ci_rank + 1) * totNumPDN / ci_numProce;
	/*startP = 0;
	endP = totNumPDN;*/
	//=====
	double delta_k, delta_j, maxDel;
	double MPx[3], cx[3];//MPx, cx,are coordinate of  center point of node j and k	 respectively
	double dist;//distance between MPx, cx
	int numNodeoB, NodeIDoB, numDimen; 
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
							dist = ((MPx[0] - cx[0]) * (MPx[0] - cx[0]) + (MPx[1] - cx[1]) * (MPx[1] - cx[1])
								+ (MPx[2] - cx[2]) * (MPx[2] - cx[2]));
							if (dist > maxDel* maxDel * 1.0e-16 &&( dist < (delta_k* delta_k + maxDel* maxDel * 1.0e-16)|| dist < (delta_j* delta_j + maxDel* maxDel * 1.0e-16)))
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
	if (cip_ia)
	{
		delete[] cip_ia;
		cip_ia = NULL;
	}
	cip_ia = new long long int[numEq + 1];
	for (int i = 0; i < numEq+1; i++)
	{
		//zero-based, ia[0]=0 always.
		cip_ia[i] = 0;
	}
	//============ find interactions of each node=========;
	vector<int>* inteNode = new vector<int> [totNumNode];
	//====PD family===;
	int NID_k;
	vector<int> famV;
	for (int famk = 0; famk < o_dat.getTotnumFami(); famk++)
	{
		NID_k = o_dat.op_getFami(famk)->getNodeID(0);
		//famV = o_dat.op_getFami(famk)->GetvecNid();
		inteNode[NID_k - 1].assign(o_dat.op_getFami(famk)->civ_NID.begin(), o_dat.op_getFami(famk)->civ_NID.end());
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
			for (int i = 0; i < numNele; i++)
			{
				if (i != nd)
				{
					inteNode[conNID[i] - 1].insert(inteNode[conNID[i] - 1].end(), o_dat.op_getFami(famID - 1)->civ_NID.begin(), o_dat.op_getFami(famID - 1)->civ_NID.end());
				}
			}
			/*for (int m = 0; m < numNodeOfFam; m++)
			{
				NID_m = o_dat.op_getFami(famID - 1)->getNodeID(m);
				for (int i = 0; i < numNele; i++)
				{
					if (i!=nd)
					{
						inteNode[conNID[i] - 1].push_back(NID_m);
					}
				}
			}*/
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
					inteNode[(conNID[nd] - 1)].push_back(conNID[i]);
				}
			}
			delete[] conNID;
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
	int eqInd_row, eqInd_col = -10, sz;
	vector<int>vec_ja;
	for (int nd = 0; nd < totNumNode; nd++)
	{
		sz= inteNode[nd].size();
		for (int i = 0; i < numDime; i++)
		{
			eqInd_row = o_dat.op_getNode(nd)->op_getDof(i)->i_getEqInde();
			if (eqInd_row !=-1)
			{
				for (int itN = 0; itN < sz; itN++)
				{
					for (int j = 0; j < numDime; j++)
					{
						eqInd_col = o_dat.op_getNode(inteNode[nd].at(itN) - 1)->op_getDof(j)->i_getEqInde();
						if (eqInd_col!=-1)
						{
							vec_ja.push_back(eqInd_col);//zero-based;
							// cip_ia[0]=0, always for zero-based, so eqInd_row+1;
							//count the num of non-zero coeff of each row
							// will be revised later;
							cip_ia[eqInd_row + 1] ++;
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
	if (cip_ja)
	{
		delete[] cip_ja;
		cip_ja = NULL;
	}
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
	if (temp1*temp2>0)// if p1 and p2 are at the same side of the crack plane
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
			Xisc[i] = xp1[i] + d * tempV2[i];// the intersecting point of p1p2 with the crack plane
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
	//int numDime = o_dat.ci_Numdimen;
	int numCrac;
	numCrac = o_dat.ci_numCrack;
	//double(*crack)[3][3] = o_dat.cdp_crack;
	//double Xmin[3], Xmax[3], xk[3], xm[3];
	////===block and family
	//double blockSize, lbc[3], temp[3];
	//blockSize = o_dat.op_getGeomP()->getBlockSize();
	//o_dat.op_getGeomP()->getlbc(lbc);
	//int i_bIndexMin[3], i_bIndexMax[3], numBlocks[3], blockIndex, numNodeoB, NodeIDoB;
	//int famID, numNodeFami, NID_k, NID_m;
	//int fn = 0;
	//if (numDime==3)
	//{
	//	fn = 3;
	//}
	//else if (numDime==2)
	//{
	//	fn = 2;
	//}
	//bool crosCrack;
	//o_dat.getNumOfBLock(numBlocks);
	//if (numDime==3)
	{
		for (int cidx = 0; cidx < numCrac; cidx++)
		{

			updateBondstate(o_dat.cdp_crack[cidx], o_dat);
			/*for (int kk = 0; kk < 3; kk++)
			{
				temp[0] = crack[cidx][0][kk];
				temp[1] = crack[cidx][1][kk];
				temp[2] = crack[cidx][2][kk];
				Xmin[kk] = *min_element(temp, temp + fn);
				Xmax[kk] = *max_element(temp, temp + fn);
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
									if (numDime==3)
									{
										crosCrack = segmentPlaneIntersection(xk, xm, crack[cidx]);
									}
									else if (numDime == 2)
									{
										crosCrack=intersection(xk, xm, crack[cidx][0], crack[cidx][1]);
									}
									if (crosCrack)
									{
										o_dat.op_getFami(famID - 1)->setbondstate(m, 0);
									}
								}
							}
						}
					}
				}
			}*/
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
	double Aa = (o_dat.cd_NLF) * (o_dat.cd_NLF);
	return (exp(-(xi[0] * xi[0] + xi[1] * xi[1] + xi[2] * xi[2])
		/ Aa / delt / delt));
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
	//====xy, yz, zx======
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

void pdsolve::assembleInterWorkPD_CSRformat(datModel& o_dat, double* U_N)
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
							if (o_dat.ci_solvFlag == 1)
							{
								// static solver;
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								//===Ku and Ku*du;
								if (eqindex_col != -1)
								{
									i_temp = findCSRIndexOfMat(eqindex_row, eqindex_col);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + j];// caution: here U_N is du;
								//	cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								//}
							}
							else
							{
								// dynamic solver;
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								}
								else
								{
									//New mark method
									//Residule force
									tempu = U_N[(NID_m - 1) * numDime + j];
									cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
									//for KN
									eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
									if (eqindex_col != -1)
									{
										i_temp = findCSRIndexOfMat(eqindex_row, eqindex_col);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
								}
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
							if (o_dat.ci_solvFlag == 1)
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								//===Ku and Ku*du;
								if (eqindex_col != -1)
								{
									i_temp = findCSRIndexOfMat(eqindex_row, eqindex_col);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + j];// caution: here U_N is du;
								//	cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								//}
							}
							else
							{
								// dynamic solver;
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
								}
								else
								{
									//newmark method
									//for residule force
									tempu = U_N[(NID_m - 1) * numDime + j];
									cdp_F[eqindex_row] = cdp_F[eqindex_row] - temp * tempu;
									//for KN
									eqindex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
									if (eqindex_col != -1)
									{
										i_temp = findCSRIndexOfMat(eqindex_row, eqindex_col);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
								}
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

void pdsolve::assemblePDBEwork3D_CSRformat(datModel& o_dat, double* U_N)
{
	int numPDBEs, startP, endP;
	numPDBEs = o_dat.getTotnumPDBEs();
	startP = ci_rank * numPDBEs / ci_numProce;
	endP = (ci_rank + 1) * numPDBEs / ci_numProce;
	
	//====================3D quad 4N PDBE===================
	int conNID[4];
	double xN[4][3], N[4];//N[4] are shape functions;
	Matrix* mNt, * Nmat, * f_mNt, * f_mNt_Nm, * f_mNt_Nm_D;
	mNt = new Matrix(12, 3);
	Nmat = new Matrix(3, 6);
	f_mNt = new Matrix(12, 3);
	f_mNt_Nm = new Matrix(12, 6);
	f_mNt_Nm_D = new Matrix(12, 6);
	//====================3D tri 3N PDBE===================
	Matrix* tri_mNt[3], * tri_mNtN[3],  *tri_mNt_N_D;
	for (int i = 0; i < 3; i++)
	{
		tri_mNt[i] = new Matrix(9, 3);
		tri_mNtN[i] = new Matrix(9, 6);
	}
	matMathcalNt_trian(tri_mNt);
	tri_mNt_N_D = new Matrix(9, 6);
	//===========3D=======================================
	for (int pdbe = startP; pdbe < endP; pdbe++)
	{
		//==get PDBEs==
		o_dat.op_getPDBE(pdbe)->getNodeID(conNID);
		if (conNID[3]!=conNID[2])
		{
			//====================3D quad 4N PDBE===================
			assemblePDBEworkQuad4NElement_CSRformat(o_dat, U_N, conNID,xN, N, mNt, Nmat, f_mNt, f_mNt_Nm, f_mNt_Nm_D);
		}
		else
		{
			//====================3D tri 3N PDBE===================
			assemblePDBEworkTri3NElement_CSRformat(o_dat, U_N, conNID, xN, tri_mNt, tri_mNtN, tri_mNt_N_D);
		}

	}
	//====================3D tri 3N PDBE===================
	delete mNt, delete Nmat, delete f_mNt, delete f_mNt_Nm, delete f_mNt_Nm_D;
	mNt = NULL, Nmat = NULL, f_mNt = NULL, f_mNt_Nm = NULL, f_mNt_Nm_D = NULL;
	//====================3D tri 3N PDBE===================
	for (int i = 0; i < 3; i++)
	{
		delete tri_mNt[i], delete tri_mNtN[i];
		tri_mNt[i] = NULL, tri_mNtN[i] = NULL;
	}
	delete  tri_mNt_N_D;
	tri_mNt_N_D = NULL;
}

void pdsolve::assemblePDBEworkQuad4NElement_CSRformat(datModel& o_dat, double* U_N,int conNID[],double xN[][3], double N[],Matrix *mNt,Matrix * Nmat, Matrix* f_mNt, Matrix* f_mNt_Nm, Matrix* f_mNt_Nm_D)
{
	pdFamily* temP_fami;
	int nG = o_globGP.i_getNumPts();
	int famiID, numNodeFam, eqIndex_row, eqIndex_col, NID_m;
	int numDime = o_dat.ci_Numdimen;
	long long int i_temp;
	double p, q, wp, wq, fac, temp, tempu;
	Matrix* C, *finMat;
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
									if (o_dat.ci_solvFlag == 1)
									{
										//static solver
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
									else if (o_dat.ci_solvFlag == 2)
									{
										//quasi-static solver;
										eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
										//==internal force==may modified later;
										tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
										cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
										//===Ku and Ku*du;
										if (eqIndex_col != -1)
										{
											i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
											cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
										}
										//else
										//{
										//	//no need this part, as du is added into value of dof;
										//	//= may modified later;
										//	tempu = U_N[(NID_m - 1) * numDime + jj];// caution: here U_N is du;
										//	cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
										//}
									}
									else
									{
										//dynamic solver
										if (!o_dat.cb_Newmark)
										{
											tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
											cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
										}
										else
										{
											//Newmark method
											//tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
											// for residule force
											tempu = U_N[(NID_m - 1) * numDime + jj];
											cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
											//for KN
											eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
											if (eqIndex_col != -1)
											{
												i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
												cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
											}
										}
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

void pdsolve::assemblePDBEworkTri3NElement_CSRformat(datModel& o_dat, double* U_N, int conNID[], double xN[][3], Matrix* mNt[], Matrix* mNtN[], Matrix* mNt_N_D)
{
	Matrix* C, * finMat;
	pdFamily* temP_fami;
	int famiID, numNodeFam, eqIndex_row, eqIndex_col;
	int NID_m;
	int numDime = o_dat.ci_Numdimen;
	long long int i_temp;
	double  temp, tempu;
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
							if (o_dat.ci_solvFlag == 1)
							{
								//static solver
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
								cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
								//===Ku and Ku*du;
								if (eqIndex_col != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + jj];// caution: here U_N is du;
								//	cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
								//}
							}
							else
							{
								//dynamic solver
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
									cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
								}
								else
								{
									//Newmark method;
									//tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
									//for residule force;
									tempu = U_N[(NID_m - 1) * numDime + jj];
									cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
									//for KN
									eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
									if (eqIndex_col != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
								}
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

void pdsolve::DispGrad(double DG[], int NID, datModel& o_dat)
{
	// this function is to get displacement gradient;
	//DG[0] ~[3] =pUxpx,pUypy, pUxpy,pUypx;
	int famk = o_dat.op_getNode(NID - 1)->getFamID() - 1;
	pdFamily* p_fam = o_dat.op_getFami(famk);
	int numNodeFam = p_fam->getNumNode();
	Matrix* matDG = new Matrix(4, 2 * numNodeFam);
	matDG->zero();

	Matrix* A;
	A = new Matrix(5, 5);
	//shapTens(A, famk, o_dat, test);
	shapTens2D(A, p_fam, o_dat);

	double g[2], d[3];

	double xk[3], xm[3], xi[3], dv_m, temp, omega;
	int NID_k, NID_m, mu_km;
	NID_k = p_fam->getNodeID(0);
	o_dat.op_getNode(NID_k - 1)->getcoor(xk);


	for (int m = 1; m < numNodeFam; m++)
	{
		NID_m = p_fam->getNodeID(m);
		o_dat.op_getNode(NID_m - 1)->getcoor(xm);
		dv_m = o_dat.op_getNode(NID_m - 1)->getvolume();
		for (int i = 0; i < 3; i++)
		{
			xi[i] = xm[i] - xk[i];
		}
		mu_km = p_fam->getbondstate(m);
		omega = inflFunc(xi, p_fam, o_dat);
		//vec_gd(g, d, A, famk, xi, o_dat, test);
		vec_gd2D(g, d, A, p_fam, xi, o_dat);
		temp = mu_km * omega * g[0] * dv_m;
		matDG->setCoeff(0, 2 * m, temp);
		matDG->setCoeff(3, 2 * m + 1, temp);
		matDG->addCoeff(0, 0, -temp);
		matDG->addCoeff(3, 1, -temp);

		temp = mu_km * omega * g[1] * dv_m;
		matDG->setCoeff(1, 2 * m + 1, temp);
		matDG->setCoeff(2, 2 * m, temp);
		matDG->addCoeff(1, 1, -temp);
		matDG->addCoeff(2, 0, -temp);
	}

	Vector* v_DG, * Uf;
	v_DG = new Vector(4);
	Uf = new Vector(2 * numNodeFam);
	for (int m = 0; m < numNodeFam; m++)
	{
		NID_m = o_dat.op_getFami(famk)->getNodeID(m);
		for (int i = 0; i < 2; i++)
		{
			temp = o_dat.op_getNode(NID_m - 1)->op_getDof(i)->d_getValue();
			Uf->setCoeff(2 * m + i, temp);
		}
	}
	matoperat.matMultiply(matDG, Uf, v_DG);
	for (int i = 0; i < 4; i++)
	{
		DG[i] = v_DG->d_getCoeff(i);
	}
	delete A, delete v_DG, delete Uf, delete matDG;
	A = NULL, v_DG = NULL, Uf = NULL, matDG = NULL;
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
		double xL1[3], xL2[3], nx, ny, temp, tempu;
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

void pdsolve::assemblePDBEwork_CSRformat(datModel& o_dat,double* U_N)
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
		double xL1[3], xL2[3], nx, ny, temp, tempu;
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
							if (o_dat.ci_solvFlag == 1)
							{
								//static solver
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//===Ku and Ku*du;
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + j];// caution: here U_N is du;
								//	cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//}
							}
							else
							{
								//dynamic solver;
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
								else
								{
									//Newmark method
									//for residule force
									tempu = U_N[(NID_m - 1) * numDime + j];
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
									//for KN
									eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
									if (eqIndex_col[j] != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
								}
								
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
							if (o_dat.ci_solvFlag == 1)
							{
								//static solver;
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//===Ku and Ku*du;
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + j];// caution: here U_N is du;
								//	cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//}
							}
							else
							{
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
								else
								{
									//Newmark method
									//for residule force
									tempu = U_N[(NID_m - 1) * numDime + j];
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
									//for KN
									eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
									if (eqIndex_col[j] != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}

								}
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
							if (o_dat.ci_solvFlag == 1)
							{
								//static solver;
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//===Ku and Ku*du;
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + j];// caution: here U_N is du;
								//	cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//}
							}
							else
							{
								//dynamic solver;
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
								else
								{
									//Newmark method
									//for residule force
									tempu = U_N[(NID_m - 1) * numDime + j];
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
									//for KN
									eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
									if (eqIndex_col[j] != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}

									
								}
								/*tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;*/
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
							if (o_dat.ci_solvFlag == 1)
							{
								//static solver;
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
							else if (o_dat.ci_solvFlag == 2)
							{
								//quasi-static solver;
								eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
								//==internal force==may modified later;
								tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//===Ku and Ku*du;
								if (eqIndex_col[j] != -1)
								{
									i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
									cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
								}
								//else
								//{
								//	//no need this part, as du is added into value of dof;
								//	//= may modified later;
								//	tempu = U_N[(NID_m - 1) * numDime + j];// caution: here U_N is du;
								//	cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								//}
							}
							else
							{
								// dynamic solver;
								if (!o_dat.cb_Newmark)
								{
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
								}
								else
								{
									//Newmark method
									//for residule force
									tempu = U_N[(NID_m - 1) * numDime + j];
									cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;
									//for KN
									eqIndex_col[j] = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->i_getEqInde();
									if (eqIndex_col[j] != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row[i], eqIndex_col[j]);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}

									
								}
								/*tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(j)->d_getValue();
								cdp_F[eqIndex_row[i]] = cdp_F[eqIndex_row[i]] - temp * tempu;*/
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
		//assemblePDBEworkQuad_CSRformat(o_dat, U_N);
		//assemblePDBEworkTetrahe_CSRformat(o_dat, U_N);
		/*if (o_dat.ci_eleType==12)
		{
			assemblePDBEworkQuad_CSRformat(o_dat,U_N);
		}
		else if (o_dat.ci_eleType==10)
		{
			assemblePDBEworkTetrahe_CSRformat(o_dat, U_N);
		}
		else
		{
			printf("ERROR: element Type for 3D must be 10 or 12.\n");
			exit(0);
		}*/
		assemblePDBEwork3D_CSRformat(o_dat, U_N);
	}
}

void pdsolve::assemblePDBEworkQuad_CSRformat(datModel& o_dat, double* U_N)
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
										if (o_dat.ci_solvFlag == 1)
										{
											//static solver
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
										else if (o_dat.ci_solvFlag == 2)
										{
											//quasi-static solver;
											eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
											//==internal force==may modified later;
											tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
											cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
											//===Ku and Ku*du;
											if (eqIndex_col != -1)
											{
												i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
												cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
											}
											//else
											//{
											//	//no need this part, as du is added into value of dof;
											//	//= may modified later;
											//	tempu = U_N[(NID_m - 1) * numDime + jj];// caution: here U_N is du;
											//	cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
											//}
										}
										else
										{
											//dynamic solver
											if (!o_dat.cb_Newmark)
											{
												tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
												cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
											}
											else
											{
												//Newmark method
												//tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
												// for residule force
												tempu = U_N[(NID_m - 1) * numDime + jj];
												cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
												//for KN
												eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
												if (eqIndex_col != -1)
												{
													i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
													cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
												}
											}
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

void pdsolve::assemblePDBEworkTetrahe_CSRformat(datModel& o_dat, double* U_N)
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
								if (o_dat.ci_solvFlag == 1)
								{
									//static solver
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
								else if (o_dat.ci_solvFlag == 2)
								{
									//quasi-static solver;
									eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
									//==internal force==may modified later;
									tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
									cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
									//===Ku and Ku*du;
									if (eqIndex_col != -1)
									{
										i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
									//else
									//{
									//	//no need this part, as du is added into value of dof;
									//	//= may modified later;
									//	tempu = U_N[(NID_m - 1) * numDime + jj];// caution: here U_N is du;
									//	cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
									//}
								}
								else
								{
									//dynamic solver
									if (!o_dat.cb_Newmark)
									{
										tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
										cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
									}
									else
									{
										//Newmark method;
										//tempu = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->d_getValue();
										//for residule force;
										tempu = U_N[(NID_m - 1) * numDime + jj];
										cdp_F[eqIndex_row] = cdp_F[eqIndex_row] - temp * tempu;
										//for KN
										eqIndex_col = o_dat.op_getNode(NID_m - 1)->op_getDof(jj)->i_getEqInde();
										if (eqIndex_col != -1)
										{
											i_temp = findCSRIndexOfMat(eqIndex_row, eqIndex_col);
											cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
										}
									}
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

void pdsolve::assembleSEDbyFEM_CSRformat(datModel& o_dat, double* U_N)
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
								if (o_dat.ci_solvFlag==1)
								{
									//cout << "== " << ci_solvFlag << endl;
									//static solver;
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
								else if (o_dat.ci_solvFlag==2)
								{
									//quasi-static solver;
									eqInde_j = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->i_getEqInde();
									//==internal force==may modified later;
									tempu = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->d_getValue();
									cdp_F[eqInde_i] = cdp_F[eqInde_i] - temp * tempu;
									//===Ku and Ku*du;
									if (eqInde_j != -1)
									{
										i_temp = findCSRIndexOfMat(eqInde_i, eqInde_j);
										cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;
									}
									//else
									//{
									//	//no need this part, as du is added into value of dof;
									//	//= may modified later
									//	tempu = U_N[(NID_j - 1) * numDime + jj];// caution: here U_N is du;
									//	cdp_F[eqInde_i] = cdp_F[eqInde_i] - temp * tempu;
									//}
									
								}
								else
								{
									//dynamic solver;
									if (!o_dat.cb_Newmark)
									{
										tempu = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->d_getValue();
										cdp_F[eqInde_i] = cdp_F[eqInde_i] - temp * tempu;
									}
									else
									{
										// if is newmark's method;
										//tempu = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->d_getValue();
										//===internal force
										tempu = U_N[(NID_j - 1) * numDime + jj];
										cdp_F[eqInde_i] = cdp_F[eqInde_i] - temp * tempu;
										//===KN=================
										eqInde_j = o_dat.op_getNode(NID_j - 1)->op_getDof(jj)->i_getEqInde();
										if (eqInde_j != -1)
										{
											i_temp = findCSRIndexOfMat(eqInde_i, eqInde_j);
											cdp_Ku[i_temp] = cdp_Ku[i_temp] + temp;// acceleration is 0;
										}
									}
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

void pdsolve::pdfemSolver(datModel& o_dat, fioFiles& o_files, char* argv[])
{
	if (o_dat.ci_solvFlag==1)
	{
		//static solver
		pdfemStaticSolver_CSRformat(o_dat, o_files,argv);
	}
	else if (o_dat.ci_solvFlag==0)
	{
		//dynamic solver;
		if (o_dat.cb_Newmark)
		{
			pdfemDynamicNewmarkSolver_CSRformat(o_dat, o_files,argv);
		}
		else
		{
			pdfemDynamicSolver_CSRformat(o_dat, o_files, argv);
		}
		
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

void pdsolve::timeIntegration_LM(datModel& o_dat,Vector* Vu_n, Vector* Vu_nm1, Vector* Vu_np1, int numEq)
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

void pdsolve::timeIntegration_NLM(datModel& o_dat, double* dA, Vector* Vu_n, Vector* Vu_nm1, Vector* Vu_np1, int numEq)
{
	int  startP, endP;
	startP = ci_rank * numEq / ci_numProce;
	endP = (ci_rank + 1) * numEq / ci_numProce;
	double tempA;
	double dt = o_dat.getTstep();
	for (int i = startP; i < endP; i++)
	{
		tempA = dt * dt * dA[i];
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
	double xc[3], tempu;
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

void pdsolve::calReacForc(datModel& o_dat, ofstream& fout, int stepN)
{
	/*for calculate node reaction force;
	only calculate total force of all node;
	can add part to obtain each node force;*/
	Matrix* Ke;
	Vector* Ue, * Fe;
	double tempu = 0;
	double NodeF = 0;
	double(*xN)[3];
	double LocTotReaF[3] = { 0 }, gloTotReaF[3] = { 0 };
	int ele, numNode, * conNid, numDime;
	numDime = o_dat.ci_Numdimen;
	int startP, endP;
	startP = ci_rank * o_dat.ci_numReacForceEle / ci_numProce;
	endP = (ci_rank + 1) * o_dat.ci_numReacForceEle / ci_numProce;
	for (int k = startP; k < endP; k++)
	{
		ele = o_dat.cop2_reacForceEle[k]->ci_eleIDX;
		numNode = o_dat.op_getEles(ele)->getNumNodes();
		conNid = new int[numNode];
		xN = new double[numNode][3];
		Ke = new Matrix(numDime * numNode, numDime * numNode);
		Fe = new Vector(numDime * numNode);
		Ue = new Vector(numDime * numNode);
		o_dat.op_getEles(ele)->getConNid(conNid);
		for (int nd = 0; nd < numNode; nd++)
		{
			o_dat.op_getNode(conNid[nd] - 1)->getcoor(xN[nd]);
			for (int j = 0; j < numDime; j++)
			{
				tempu = o_dat.op_getNode(conNid[nd] - 1)->op_getDof(j)->d_getValue();
				Ue->setCoeff(numDime * nd + j, tempu);
			}
		}
		o_dat.op_getEles(ele)->eleStiffMatFEM(Ke, cop_D, xN);
		matoperat.matMultiply(Ke, Ue, Fe);
		for (auto iter = o_dat.cop2_reacForceEle[k]->civ_reaForceEleNidx.begin(); iter!= o_dat.cop2_reacForceEle[k]->civ_reaForceEleNidx.end(); iter++)
		{
			for (int j = 0; j < numDime; j++)
			{
				LocTotReaF[j] = LocTotReaF[j] + Fe->d_getCoeff((*iter) * numDime + j);
			}
		}
		delete[] conNid, delete[]xN;
		delete Ke, delete Fe, delete Ue;
	}
	MPI_Reduce(LocTotReaF, gloTotReaF, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if (ci_rank==0)
	{
		double uc = 0;
		if (!o_dat.civ_reacForceOfessBCId.empty())
		{
			uc = o_dat.op_getEssenBC(o_dat.civ_reacForceOfessBCId.at(0))->cd_value;
		}
		fout << stepN << "\t" <<uc << "\t" << gloTotReaF[0] << "\t" << gloTotReaF[1] << "\t" << gloTotReaF[2] << endl;
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
	//cb_InteralForce = false;
	unique_ptr<double[]> U_N = make_unique<double[]>(1);
	setCSRIndexes_gloStiffMat(o_dat);
	assembleInterWorkPD_CSRformat(o_dat,U_N.get());
	assemblePDBEwork_CSRformat(o_dat, U_N.get());
	assembleSEDbyFEM_CSRformat(o_dat, U_N.get());
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

void pdsolve::pdfemStaticSolver_CSRformat(datModel& o_dat, fioFiles& o_files,char* argv[])
{
	//argv[2]: number of nodes,
	//argv[3]: test ID;
	//argv[4]: V_pd fraction;
	//Don't need block data any more for static solver;
	//o_dat.deleteBLOCK();
	double t1, t2;
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
		t1 = MPI_Wtime();
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
		//t2 = MPI_Wtime();
		//double memo = (((numEq + 1) + cip_ia[numEq]) * 8.0 + cip_ia[numEq] * 8.0) / 1024;
		//printf("DATA:time, %f, memo(kb), %f,numNon-zero, %d, testID, %d,nodes, %d, cores, %d, mRatio, %f, Vpdfrac, %d \n", t2 - t1, memo, cip_ia[numEq], atoi(argv[3]), atoi(argv[2]), ci_numProce,o_dat.op_getGeomP()->cd_factor, atoi(argv[4]));
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
	calLocalDamage(o_dat);
	bool addLoad;
	ofstream test;
	double maxv = failureProcess(o_dat, 8, addLoad, test);
	if (ci_rank == 0)
	{
		printf("MAX Value %e \n", maxv);
		//===========post processing=====================
		//ofstream cracPathOut;
		//cracPathOut.open("test.out");
		//b_cracPropag_qusiaStatic_BYKeq(o_dat, cracPathOut);//calculate KI,KII;
		printf("writing results........\n");
		o_files.writeResults(o_dat);
	}
	//====for calculate reaction force;
	/*string RFfilename = o_dat.cs_fileName + "RF.out";
	ofstream fout(RFfilename);
	fout << "Step\t" << "U\t" << "RFx\t" << "RFy\t" << "RFz\n";
	calReacForc(o_dat, fout, 0);
	fout.close();*/
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

void pdsolve::calinternalForce_CSRformat(datModel& o_dat, int numEq,double *U_N)
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
	//cb_InteralForce = true;
	assembleInterWorkPD_CSRformat(o_dat, U_N);
	assemblePDBEwork_CSRformat(o_dat, U_N);
	assembleSEDbyFEM_CSRformat(o_dat, U_N);
	//calExternalForce_CSRformat(o_dat);
	//MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void pdsolve::pdfemDynamicSolver_CSRformat(datModel& o_dat, fioFiles& o_files, char* argv[])
{
	int Tk = o_dat.ci_topk;
	bool addLoad;
	//====for calculate reaction force;
	string RFfilename = o_dat.cs_fileName + "RF.out";
	ofstream fout;
	if (ci_rank == 0)
	{
		fout.open(RFfilename);
		fout << "Step\t" << "U\t" << "RFx\t" << "RFy\t" << "RFz\n";
	}
	////===========static solver to set the start ebc of dynamic solver;
	//o_dat.ci_solvFlag = 1;
	//pdfemStaticSolver_CSRformat(o_dat, o_files, argv); // stress calculated;
	//calReacForc(o_dat, fout, -2);
	//double ratio=1, maxVal;
	//ofstream cracPath;
	//maxVal = failureProcess(o_dat, Tk, addLoad, cracPath);
	//if (o_dat.ci_failFlag == 0)
	//{
	//	double  KIC, G0, kappa, E, nu, delta, sc;
	//	E = o_dat.op_getmaterial()->getE();
	//	nu = o_dat.op_getmaterial()->getnu();
	//	KIC = o_dat.op_getmaterial()->getKIc();
	//	G0 = KIC * KIC / E;
	//	kappa = E / (3.0 * (1 - 2 * nu));
	//	delta = o_dat.op_getGeomP()->getmaxDelta();
	//	sc = sqrt(5 * G0 / (9 * kappa * delta));
	//	ratio = sc / maxVal * 1.0001;

	//}
	//else if (o_dat.ci_failFlag == 1)
	//{
	//	//===get ratio==
	//	ratio = o_dat.op_getmaterial()->getSigult() / maxVal * 1.0001;
	//}
	//else if (o_dat.ci_failFlag == 3)
	//{
	//	// maximum principal stress
	//	ratio = o_dat.op_getmaterial()->getSigult() / maxVal * 1.0001;
	//}
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
	//pdfemStaticSolver_CSRformat(o_dat, o_files, argv); // stress calculated;
	//calReacForc(o_dat, fout, -1);
	//// updat bond-state
	//o_dat.ci_solvFlag = 0;
	//failureProcess(o_dat, Tk, addLoad, cracPath);
	//================================================================
	//o_dat.ci_solvFlag = 0;// remove later;
	if (o_dat.ci_solvFlag != 0)
	{
		printf("ERROR: This is not dynamic solver\n");
		printf("Please reset the solver flag\n");
		exit(0);
	}
	ofstream cracPath;
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
	//============get 1th step displacement, u0;
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
	//==set varied NBC (NBC, not prescribed displacemets); 
	setVaryNaturalBC(o_dat);
	//====resultant force;
	unique_ptr<double[]> U_N = make_unique<double[]>(1);
	calinternalForce_CSRformat(o_dat, numEq, U_N.get());
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
	//===update displacement and store, u1;
	matoperat.matMultiply(Va_n, 0.5 * dt * dt, Vu_np1);
	matoperat.matAdd(Vu_np1, Vu_n, Vu_np1);
	storeDisplacementResult(o_dat, Vu_np1);
	//calculate stress==
	calGlobalNodeStresses(o_dat);
	//updat bond-state===
	failureProcess(o_dat, Tk, addLoad, cracPath);
	//===============
	calReacForc(o_dat, fout, 0);
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
		calinternalForce_CSRformat(o_dat, numEq, U_N.get());
		calExternalForce_CSRformat(o_dat);
		MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		//===time integration
		if (o_dat.cb_lumpedMass)
		{
			timeIntegration_LM(o_dat, Vu_n, Vu_nm1, Vu_np1, numEq);//for lumed mass;
		}
		else
		{
			//===solve for acceleration
			matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_MGlo,
				cdp_FGlo, dp_A, &comm);
			MPI_Bcast(dp_A, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			////===updat displacement and store the results;
			//matoperat.matMultiply(Va_n, dt* dt, Vu_np1); // a*t^2;
			//matoperat.matAdd(Vu_np1, 2.0, Vu_n, Vu_np1);// 2*U_n+a*t^2;
			//matoperat.matMinus(Vu_np1, Vu_nm1, Vu_np1);//2*U_n+a*t^2-U_(n-1);
			timeIntegration_NLM(o_dat,dp_A, Vu_n, Vu_nm1, Vu_np1, numEq);//for not lumed mass;
		}
		//==store displacement
		storeDisplacementResult(o_dat, Vu_np1);
		//calculate stress==
		calGlobalNodeStresses(o_dat);
		//updat bond-state===
		failureProcess(o_dat, Tk, addLoad, cracPath);
		////==calculta reaction force if needed;
		calReacForc(o_dat, fout, n);
		//write results===
		if (n%saveFren==0||n== numTstep)
		{
			calLocalDamage(o_dat);
			if (ci_rank == 0)
			{
				printf("Writing results at step of %d ......\n", n);
				o_files.writeReslutsTOTAL_vtk(o_dat, to_string(n)); //**set phi
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
	fout.close();
}

void pdsolve::pdfemDynamicNewmarkSolver_CSRformat(datModel& o_dat, fioFiles& o_files, char* argv[])
{
	int Tk = o_dat.ci_topk;
	bool addLoad;
	//====for calculate reaction force;
	string RFfilename = o_dat.cs_fileName + "RF.out";
	ofstream fout;
	if (ci_rank==0)
	{
		fout.open(RFfilename);
		fout << "Step\t" << "U\t" << "RFx\t" << "RFy\t" << "RFz\n";
	}
	
	////===========static solver to set the start ebc of dynamic solver;
	//o_dat.ci_solvFlag = 1;
	//pdfemStaticSolver_CSRformat(o_dat,o_files,argv); // stress calculated;
	//calReacForc(o_dat, fout, -1);
	//double ratio = 1, maxVal;
	//maxVal=failureProcess(o_dat, Tk, addLoad, cracPath);
	//if (o_dat.ci_failFlag == 0)
	//{
	//	double  KIC, G0, kappa, E, nu, delta, sc;
	//	E = o_dat.op_getmaterial()->getE();
	//	nu = o_dat.op_getmaterial()->getnu();
	//	KIC = o_dat.op_getmaterial()->getKIc();
	//	G0 = KIC * KIC / E;
	//	kappa = E / (3.0 * (1 - 2 * nu));
	//	delta = o_dat.op_getGeomP()->getmaxDelta();
	//	sc = sqrt(5 * G0 / (9 * kappa * delta));
	//	ratio = sc / maxVal * 1.0001;
	//	
	//}
	//else if (o_dat.ci_failFlag == 1)
	//{
	//	//===get ratio==
	//	ratio = o_dat.op_getmaterial()->getSigult() / maxVal * 1.0001;
	//}
	//else if (o_dat.ci_failFlag == 3)
	//{
	//	// maximum principal stress
	//	ratio = o_dat.op_getmaterial()->getSigult() / maxVal * 1.0001;
	//}
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
	//pdfemStaticSolver_CSRformat(o_dat, o_files,argv); // stress calculated;
	//calReacForc(o_dat, fout, 0);
	//// updat bond-state
	////////??????????need here???
	//o_dat.ci_solvFlag = 0;
	//failureProcess(o_dat, Tk, addLoad, cracPath);
	///*=======================================================================*/


	//Don't need block data any more for dynamic solver;
	//o_dat.deleteBLOCK();
	//================================================================

	ofstream cracPath;
	o_dat.ci_solvFlag = 0;// dynamic solver; may remove later;
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
	double* dp_A = new double[numEq];// accerations, Va_np1 point to this pointer;
	Vector* Vu_n, * Va_n, * Va_np1, * Vv_n;
	unique_ptr<double[]>U_N = make_unique<double[]>(o_dat.getTotnumNode() * numDime);
	Vu_n = new Vector(numEq);
	Va_n = new Vector(numEq);
	Va_np1 = new Vector();
	Va_np1->setNumRows(numEq);
	Va_np1->setCoeff(dp_A);
	Vv_n = new Vector(numEq);
	//=================initialize 0th step========
	//============get 0th step displacement;
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
	Vv_n->zero();
	Va_n->zero();
	//===============start dynamic solving==========
	double dt = o_dat.getTstep();
	if (ci_rank == 0)
	{
		printf("start to Newmark's Dynamic solving......\n");
		printf("Time step is %e\n", dt);
	}
	//===start loop=============;
	int comm = MPI_Comm_c2f(MPI_COMM_WORLD);
	
	for (int n = 1; n <= numTstep; n++) //must start from n=0;
	{
		if (ci_rank == 0)
		{
			cout << "time step n= " << n << endl;
		}
		
		//==============get a_np1==========
		//get K_N;
		//pdfemAssembleKN_CSRformat( o_dat, numEq, n);
		//get Residual force, and K_N;
		pdfemAssembleKNFR_CSRformat(o_dat, numEq, Va_n, Vv_n, U_N.get(),n-1); //n=0 is the first step
		//===solve for acceleration a_np1
		matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_KuGlo,
			cdp_FGlo, dp_A, &comm);
		MPI_Bcast(dp_A, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//===update displacement and store;
		updateDispVelo(Vu_n, Vv_n, Va_n, Va_np1, numEq, o_dat);
		storeDisplacementResult(o_dat, Vu_n);
		//==set varied prescribed displacemets;
		setVaryEssentialBC(o_dat);
		//==set varied NBC (NBC, not prescribed displacemets); 
		setVaryNaturalBC(o_dat);
		//update acceleration for next step;
		double* tempPointer = Va_n->cdp_vecCoeff;
		Va_n->cdp_vecCoeff= Va_np1->cdp_vecCoeff;
		Va_np1->cdp_vecCoeff = tempPointer;
		//calculate stress====
		calGlobalNodeStresses(o_dat);
		//updat bond-state===
		failureProcess(o_dat, Tk, addLoad, cracPath);
		//==calculta reaction force if needed;
		calReacForc(o_dat, fout, n);
		/***********************/
		//write results===
		//calLocalDamage(o_dat);
		if (n % saveFren == 0 || n == numTstep)
		{
			calLocalDamage(o_dat);
			if (ci_rank == 0)
			{
				printf("Writing results at step of %d ......\n", n);
				o_files.writeReslutsTOTAL_vtk(o_dat, to_string(n)); //**set phi
			}
		}
	}
	// release memory;
	delete Vu_n, delete Va_np1, delete Va_n, delete Vv_n; // no need to delete dp_A, if Va_n is deleted;
	Vu_n=NULL,  Va_np1=NULL,  Va_n=NULL,  Vv_n=NULL;
	delete[]cdp_M, delete[]cdp_Ku, delete[]cdp_KuGlo;
	delete[] cdp_F, delete[] cdp_FGlo;
	cdp_MGlo = NULL, cdp_Ku = NULL; cdp_KuGlo = NULL;
	cdp_F = NULL, cdp_FGlo = NULL; 
	if (cip_ia)
	{
		delete[] cip_ia, delete[] cip_ja;
		cip_ia = NULL, cip_ja = NULL;
	}
	fout.close();
}

void pdsolve::pdfemAssembleKN_CSRformat(datModel& o_dat, int numEq,int n)
{
	if (n==0)
	{
		setCSRIndexes_gloStiffMat(o_dat);
		//===initializa cdp_M;
		cdp_M = new double[cip_ia[numEq]];//*****don't foget to release as well as cdp_KU
		for (long long int i = 0; i < cip_ia[numEq]; i++)
		{
			cdp_M[i] = 0;
		}
		assembleElemassMatFEM_CSRformat(o_dat);
		assembleMassMatPD_CSRformat(o_dat);
	}
	//===assemble K====
	for (long long int i = 0; i < cip_ia[numEq]; i++)
	{
		cdp_Ku[i] = 0;
		cdp_KuGlo[i] = 0;
	}
	//cb_InteralForce = false;
	double dt =o_dat.cd_dt;
	unique_ptr<double[]> U_N = make_unique<double[]>(1);
	assembleInterWorkPD_CSRformat(o_dat, U_N.get()); 
	assemblePDBEwork_CSRformat(o_dat, U_N.get());
	assembleSEDbyFEM_CSRformat(o_dat, U_N.get());
	//==============add M====
	for (long long int i = 0; i < cip_ia[numEq]; i++)
	{
		cdp_Ku[i] = cdp_Ku[i]*cd_beta*dt*dt+cdp_M[i];
	}
	MPI_Allreduce(cdp_Ku, cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void pdsolve::pdfemAssembleKNFR_CSRformat(datModel& o_dat, int numEq, Vector* a_n, Vector* V_n,double *U_N,int n)
{
	//===get U_N====================================
	int numDime = o_dat.ci_Numdimen;
	int eqIdx;
	double dt = o_dat.cd_dt;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		for (int j = 0; j < numDime; j++)
		{
			eqIdx = o_dat.op_getNode(i)->op_getDof(j)->i_getEqInde();
			if (eqIdx!=-1)
			{
				U_N[i * numDime + j] = o_dat.op_getNode(i)->op_getDof(j)->d_getValue() +
					dt * V_n->d_getCoeff(eqIdx) + dt * dt * (0.5 - cd_beta) * a_n->d_getCoeff(eqIdx);
			}
			else
			{
				U_N[i * numDime + j] = o_dat.op_getNode(i)->op_getDof(j)->d_getValue();
			}
		}
	}
	//===varing ebc=============================================;
	int numEBCs, numNODE, i_dof, Nid;
	numEBCs = o_dat.getTotnumEssentialBCs();
	double dU;
	for (int i = 0; i < numEBCs; i++)
	{
		if (o_dat.op_getEssenBC(i)->cb_varing)
		{
			dU = dt * o_dat.op_getEssenBC(i)->cd_velocity;
			numNODE = o_dat.op_getEssenBC(i)->getNumNODE();
			i_dof = o_dat.op_getEssenBC(i)->getDOF();
			for (int j = 0; j < numNODE; j++)
			{
				Nid = o_dat.op_getEssenBC(i)->cip_NID[j];
				U_N[(Nid - 1) * numDime + i_dof] = U_N[(Nid - 1) * numDime + i_dof] + dU;
			}
		}
	}
	//=================get M, and allocate memory for cdp_Ku,cdp_Kuglo, only call it for the 1st step;
	if (n == 0)
	{
		setCSRIndexes_gloStiffMat(o_dat);
		//===initializa cdp_M;
		cdp_M = new double[cip_ia[numEq]];//*****don't foget to release as well as cdp_KU
		for (long long int i = 0; i < cip_ia[numEq]; i++)
		{
			cdp_M[i] = 0;
		}
		assembleElemassMatFEM_CSRformat(o_dat);
		assembleMassMatPD_CSRformat(o_dat);
	}
	for (long long int i = 0; i < cip_ia[numEq]; i++)
	{
		cdp_Ku[i] = 0;
		cdp_KuGlo[i] = 0;
	}
	//===============get internalforce and KN ===========================
	//cb_InteralForce = true;
	calinternalForce_CSRformat(o_dat, numEq, U_N);// This function also get part of KN in newmark's method;
	//may need set varing nbc pbc;
	calExternalForce_CSRformat(o_dat);
	MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//==============add M==========
	for (long long int i = 0; i < cip_ia[numEq]; i++)
	{
		cdp_Ku[i] = cdp_Ku[i] * cd_beta * dt * dt + cdp_M[i];
	}
	MPI_Allreduce(cdp_Ku, cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void pdsolve::updateDispVelo(Vector* Vu_n, Vector* Vv_n, Vector* Va_n, Vector* Va_np1, int numEq, datModel& o_dat)
{
	double temp1, dt, temp2;
	dt = o_dat.cd_dt;
	for (int i = 0; i < numEq; i++)
	{
		temp1 = Vu_n->d_getCoeff(i) + dt * Vv_n->d_getCoeff(i) +
			dt * dt * ((0.5 - cd_beta) * Va_n->d_getCoeff(i) + cd_beta * Va_np1->d_getCoeff(i));
		temp2 = Vv_n->d_getCoeff(i) + dt * ((1 - cd_gamma) * Va_n->d_getCoeff(i) + cd_gamma * Va_np1->d_getCoeff(i));
		Vu_n->setCoeff(i, temp1);//update displacement;
		Vv_n->setCoeff(i, temp2);//update velocity;
	}
}

void pdsolve::printArr(double* V, int numEq, ofstream& fout)
{
	fout << "numEq " << numEq << endl;
	for (int i = 0; i < numEq; i++)
	{
		fout << V[i] << endl;
	}
}


void pdsolve::pdfemQuasiStaticAssembleEquasSys_CSRformat(datModel& o_dat, int numEq, bool AddLoad)
{
	//here U_N is the dU;
	o_dat.ci_solvFlag = 2;
	//initial K matrix
	for (long long int i = 0; i < cip_ia[numEq]; i++)
	{
		cdp_Ku[i] = 0;
		cdp_KuGlo[i] = 0;
	}
	//===add load=======================================================
	//==================================================================
	if (AddLoad)
	{
		//==set varied NBC (NBC, not prescribed displacemets); 
		setVaryNaturalBC(o_dat);
		//==set varied prescribed displacemets;
		setVaryEssentialBC(o_dat);
	}
	//=================================================================================
	//====internal force;
	unique_ptr<double[]> U_N = make_unique<double[]>(1);
	calinternalForce_CSRformat(o_dat, numEq, U_N.get()); //This function get the internal force and the K simultaneously
	//==external force;
	calExternalForce_CSRformat(o_dat);
	//==all reduce sum;
	MPI_Allreduce(cdp_Ku, cdp_KuGlo, cip_ia[numEq], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(cdp_F, cdp_FGlo, numEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

void pdsolve::pdfemQuasiStaticSolver_CSRformat(datModel& o_dat, fioFiles& o_files)
{
	//====for calculate reaction force;
	string RFfilename = o_dat.cs_fileName + "RF.out";
	ofstream fout;
	if (ci_rank == 0)
	{
		fout.open(RFfilename);
		fout << "Step\t" << "U\t" << "RFx\t" << "RFy\t" << "RFz\n";
	}
	//===============================
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
	ofstream cracPathout;
	//===============start quasi-static solving==========
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
	bool addLoad = true;
	int Tk = o_dat.ci_topk;
	double Keq;
	//===start loop=============;
	for (int n = 1; n < numTstep + 1; n++)
	{
		if (ci_rank == 0)
		{
			printf("load step n=%d. ", n);
			if (addLoad)
			{
				printf("Add load");
			}
			printf("\n");
		}
		//assenble Euqation system;
		pdfemQuasiStaticAssembleEquasSys_CSRformat(o_dat, numEq,addLoad);

		// solver du;
		matoperat.cluster_PARDISO_64Solver(numEq_long, cip_ia, cip_ja, cdp_KuGlo,
			cdp_FGlo, Vdu->cdp_vecCoeff, &comm);

		/*store results========;*/
		MPI_Bcast(Vdu->cdp_vecCoeff, numEq, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		matoperat.matAdd(Vu_n, Vdu, Vu_n);
		storeDisplacementResult(o_dat, Vu_n);
		/*==stresses==========*/
		calGlobalNodeStresses(o_dat);
		////===updat bond-state by falure criterion=======
		failureProcess(o_dat, Tk, addLoad, cracPathout);
		//======calculta reaction force if needed;
		calReacForc(o_dat, fout, n);
		/*********/
		//write results===
		if (n % saveFren == 0 || n == numTstep)
		{
			calLocalDamage(o_dat);
			if (ci_rank == 0)
			{
				printf("Writing results at step of %d ......\n", n);
				o_files.writeReslutsTOTAL_vtk(o_dat, to_string(n)); //**set phi
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
	fout.close();
}

double pdsolve::failureProcess(datModel& o_dat, int Tk, bool& addLoad, ofstream& cracPath)
{
	if (o_dat.ci_Numdimen==3)
	{
		return 1;
	}
	double maxVal = 0;
	if (o_dat.ci_failFlag == 1)
	{
		//1--maximum circumferential tensile stress;
		addLoad = !(b_cracPropag_qusiaStatic_BYKeq(o_dat, cracPath, maxVal));
	}
	else if (o_dat.ci_failFlag == 2)
	{
		//2-- maximum principal stress: return max stress;
		maxVal = failureCriterion_maxPriSig(o_dat, Tk, addLoad);
	}
	//else if (o_dat.ci_failFlag == 3)
	//{
	//	//stetch criterion: return max stretch
	//	maxVal = failureCriterion_stretch(o_dat, Tk, addLoad);
	//}
	//else if (o_dat.ci_failFlag == 4)
	//{
	//	//equivalent stress: return max stress
	//	maxVal = failureCriterion_stress(o_dat, Tk, addLoad);
	//}
	
	return maxVal;
}

double pdsolve::failureCriterion_stretch(datModel& o_dat, int Tk, bool& addLoad)
{
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	pdFamily* temP_fami;
	/*vector<int>mp_fami(numFami);
	vector<int>mp_fami_glo(numFami);*/
	int numNodeOfFami, NID_k, NID_m, mu_km;
	double xk[3] = { 0 }, xm[3] = { 0 }, uk[3] = { 0 }, um[3] = { 0 }, xi[3] = { 0 }, eta[3] = { 0 };
	double mag_ref,mag, s,s_max=0, sc;//get latter
	double  KIC, G0, kappa, E, nu, delta_k, delta_m, delta;
	double fac = o_dat.op_getGeomP()->getFactor();
	E = o_dat.op_getmaterial()->getE();
	nu = o_dat.op_getmaterial()->getnu();
	KIC = o_dat.op_getmaterial()->getKIc();
	G0 = KIC * KIC / E;
	kappa = E / (3.0 * (1 - 2 * nu));
	priority_queue<critStru> TopK;
	double topVal;
	for (int famkk = startP; famkk < endP; famkk++)
	{
		temP_fami = o_dat.op_getFami(famkk);
		if (temP_fami->cb_allowFail)
		{
			numNodeOfFami = temP_fami->getNumNode();
			NID_k = temP_fami->getNodeID(0);
			if (o_dat.op_getNode(NID_k-1)->getLocalDamage()>0.42)
			{
				continue;
			}
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
						//delta = 0.5 * (delta_k + delta_m);
						delta = o_dat.op_getGeomP()->getmaxDelta();
						sc = sqrt(5 * G0 / (9 * kappa * delta));
						if (s > sc)
						{
							//if (Tk>0) //top K value;
							{
								critStru temStr(s, famkk, m);
								if (Tk==0|| int(TopK.size()) < Tk)
								{
									TopK.push(temStr);
								}
								else
								{
									topVal = TopK.top().sd_value;
									if (s > topVal)
									{
										TopK.pop();
										TopK.push(temStr);
									}
								}
							}
							//else
							//{
							//	temP_fami->setbondstate(m, 0);//only non-static solver set mu=0;
							//	mp_fami[famkk] = 1;
							//}
							
						}

					}
					else//static solver
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
		if (Tk>0)
		{
			addLoad = Finally_TopK_andUpdate_bondStaus(o_dat, Tk, TopK, 1);//flag==1
		}
		/*else
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
						MPI_Bcast(&(o_dat.op_getFami(famk)->cip_bondState[0]), numNodeOfFami, MPI_INT, rank, MPI_COMM_WORLD);
					}
				}
			}
		}*/
		
	}
	else if (o_dat.ci_solvFlag==1)//static solver
	{
		MPI_Allreduce(&s_max, &maxS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	return maxS;
}

double pdsolve::failureCriterion_stress(datModel& o_dat, int Tk, bool& addLoad)
{
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	//startP = 0;
	//endP = numFami;

	pdFamily* temP_fami;
	/*vector<int>mp_fami(numFami);
	vector<int>mp_fami_glo(numFami);*/
	int numNodeOfFami, NID_k, NID_m, mu_km;
	double sig[6], sigma_k[6], sigma_m[6], sigma1, sigult, sigma1_max = 0;
	sigult = o_dat.op_getmaterial()->getSigult();
	priority_queue<critStru> TopK;
	double topVal;
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
							critStru temStr(sigma1, famkk, m);
							if (Tk == 0 || int(TopK.size()) < Tk)
							{
								TopK.push(temStr);
							}
							else
							{
								topVal = TopK.top().sd_value;
								if (sigma1 > topVal)
								{
									TopK.pop();
									TopK.push(temStr);
								}
							}
							//temP_fami->setbondstate(m, 0);//only dynamic solver set mu=0;
							//mp_fami[famkk] = 1;
						}
					}
					else//static solver
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
		addLoad = Finally_TopK_andUpdate_bondStaus(o_dat, Tk, TopK, 1);//flag==1
		/*MPI_Allreduce(&(mp_fami[0]), &(mp_fami_glo[0]), numFami, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
		}*/
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

double pdsolve::failureCriterion_TopK_Stress(datModel& o_dat,int Tk, bool& addLoad)
{
	// equivalent stress failure criterion;
	int numDime = o_dat.ci_Numdimen;
	int numFami, startP, endP;
	numFami = o_dat.getTotnumFami();
	startP = ci_rank * numFami / ci_numProce;
	endP = (ci_rank + 1) * numFami / ci_numProce;
	//========
	int numNodeOfFami, NID_k, NID_m, mu_km;
	double  sigma_k[6], sigma_m[6], sig[6], sigma1, sigult, sigma1_max = 0, topVal;
	sigult = o_dat.op_getmaterial()->getSigult();
	pdFamily* temP_fami;
	priority_queue<critStru> TopK;
	for (int famkk = startP; famkk < endP; famkk++)
	{
		temP_fami = o_dat.op_getFami(famkk);
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
				for (int ii = 0; ii < 6; ii++)
				{
					sig[ii] = 0.5 * (sigma_k[ii] + sigma_m[ii]);
				}
				sigma1 = sqrt(0.5 * ((sig[0] - sig[1]) * (sig[0] - sig[1]) +
					(sig[0] - sig[2]) * (sig[0] - sig[2]) + (sig[1] - sig[2]) * (sig[1] - sig[2])) +
					3 * (sig[3] * sig[3] + sig[4] * sig[4] + sig[5] * sig[5]));
				if (o_dat.ci_solvFlag != 1)//non-static solver
				{
					if (sigma1 > sigult)
					{
						critStru temStr(sigma1, famkk, m);
						if (int(TopK.size())<Tk) 
						{
							TopK.push(temStr);
						}
						else
						{
							topVal = TopK.top().sd_value;
							if (sigma1>topVal)
							{
								TopK.pop();
								TopK.push(temStr);
							}
						}
					}
				}
				else//static solver
				{
					if (sigma1 > sigma1_max)
					{
						sigma1_max = sigma1;
					}
				}
			}
		}
	}
	//==== sent and receive data; and finally find the top k bonds;
	double maxSig = 0;
	addLoad = true;
	if (o_dat.ci_solvFlag != 1 && numFami > 0)//non-static solver
	{
		addLoad=Finally_TopK_andUpdate_bondStaus(o_dat, Tk, TopK,1);//flag==1
	}
	else if (o_dat.ci_solvFlag == 1)//static solver;
	{
		MPI_Allreduce(&sigma1_max, &maxSig, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	return maxSig;
}

bool pdsolve::Finally_TopK_andUpdate_bondStaus(datModel& o_dat, int Tk, priority_queue<critStru>& TopK, int flag)
{
	/*flag == 1: directly break the bond;
	flag==2: for  maximum principal stress*/
	int finaSize = 0;
	priority_queue<critStru> finaTopK;
	if (ci_rank != 0)
	{
		//send data to rank 0;
		int siz = TopK.size();
		MPI_Send(&siz, 1, MPI_INT, 0, ci_rank + 100, MPI_COMM_WORLD);
		double* dp_topSig = new double[siz];
		int* ip_famk = new int[siz];
		int* ip_m = new int[siz];
		int count = 0;
		while (!TopK.empty())
		{
			dp_topSig[count] = TopK.top().sd_value;
			ip_famk[count] = TopK.top().si_famk;
			ip_m[count] = TopK.top().si_m;
			TopK.pop();
			count++;
		}
		if (siz > 0)
		{
			MPI_Send(dp_topSig, siz, MPI_DOUBLE, 0, ci_rank + 200, MPI_COMM_WORLD);
			MPI_Send(ip_famk, siz, MPI_INT, 0, ci_rank + 201, MPI_COMM_WORLD);
			MPI_Send(ip_m, siz, MPI_INT, 0, ci_rank + 202, MPI_COMM_WORLD);
		}
		delete[] dp_topSig, delete[] ip_famk, delete[]ip_m;
	}
	else
	{
		//rank==0, recieve data and finally find the top k bonds;
		int* ip_size = new int[ci_numProce];
		ip_size[0] = TopK.size();
		int totSize = ip_size[0];
		for (int rank = 1; rank < ci_numProce; rank++)
		{
			MPI_Recv(&(ip_size[rank]), 1, MPI_INT, rank, rank + 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			totSize = totSize + ip_size[rank];
		}
		double* totSig = new double[totSize];
		int* totFamk = new int[totSize];
		int* tot_m = new int[totSize];
		//==get the data of rank 0;
		int count = 0;
		while (!TopK.empty())
		{
			totSig[count] = TopK.top().sd_value;
			totFamk[count] = TopK.top().si_famk;
			tot_m[count] = TopK.top().si_m;
			TopK.pop();
			count++;
		}
		//====get data from other ranks.
		for (int rank = 1; rank < ci_numProce; rank++)
		{
			if (ip_size[rank] > 0)
			{
				MPI_Recv(&(totSig[count]), ip_size[rank], MPI_DOUBLE, rank, rank + 200, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				MPI_Recv(&(totFamk[count]), ip_size[rank], MPI_INT, rank, rank + 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&(tot_m[count]), ip_size[rank], MPI_INT, rank, rank + 202, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				count = count + ip_size[rank];
			}

		}
		//========finally find the top k bonds;
		double topVal;
		for (int i = 0; i < totSize; i++)
		{
			critStru temStr(totSig[i], totFamk[i], tot_m[i]);
			if (Tk==0||int(finaTopK.size()) < Tk)
			{
				finaTopK.push(temStr);
			}
			else
			{
				topVal = finaTopK.top().sd_value;
				if (totSig[i] > topVal)
				{
					finaTopK.pop();
					finaTopK.push(temStr);
				}
			}
		}
		finaSize = finaTopK.size();
		delete[] ip_size;
		delete[]totSig, delete[]totFamk, delete[]tot_m;
	}
	//=============Becast the results to all rankers====;
	MPI_Bcast(&finaSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int* finaFamk = new int[finaSize];
	int* fina_m = new int[finaSize];
	if (ci_rank == 0)
	{
		int con = 0;
		while (!finaTopK.empty())
		{
			finaFamk[con] = finaTopK.top().si_famk;
			fina_m[con] = finaTopK.top().si_m;
			finaTopK.pop();
			con++;
		}
	}
	MPI_Bcast(finaFamk, finaSize, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(fina_m, finaSize, MPI_INT, 0, MPI_COMM_WORLD);
	//====updat bond status====================;
	if (flag==1) // set 
	{
		for (int i = 0; i < finaSize; i++)
		{
			if (o_dat.op_getFami(finaFamk[i])->cb_allowFail)
			{
				o_dat.op_getFami(finaFamk[i])->setbondstate(fina_m[i], 0);
			}
		}
	}
	else if (flag==2)
	{
		// this part may be parallelize in future;
		//===block and family
		double XXk[3], XXm[3], xc[3], crack[3][3] = { 0 }, maxDelta;
		maxDelta =0.5* o_dat.op_getGeomP()->getmaxDelta();
		int Nidk, Nidm;
		int  maxIDX;
		double sigma_k[6], sigma_m[6], eigV[3] = { 0 }, temp;
		int Matrc = 0, numDim = o_dat.ci_Numdimen;
		Matrc = numDim;
		Matrix* sigma_kg = new Matrix(Matrc, Matrc);
		Matrix* simga_priaxis = new Matrix(Matrc, Matrc);
		Vector* sigma_pri = new Vector(Matrc);
		for (int k = 0; k < finaSize; k++)
		{
			Nidk=o_dat.op_getFami(finaFamk[k])->getNodeID(0);
			Nidm = o_dat.op_getFami(finaFamk[k])->getNodeID(fina_m[k]);
			o_dat.op_getNode(Nidk - 1)->getcoor(XXk);
			o_dat.op_getNode(Nidm - 1)->getcoor(XXm);
			for (int i = 0; i < 3; i++)
			{
				xc[i] = 0.5 * (XXk[i] + XXm[i]);
			}
			//get stress===
			o_dat.op_getNode(Nidk - 1)->getStress(sigma_k);
			o_dat.op_getNode(Nidm - 1)->getStress(sigma_m);
			if (numDim==3)
			{
				for (int ii = 0; ii < 3; ii++)
				{
					sigma_kg->setCoeff(ii, ii, 0.5 * (sigma_k[ii] + sigma_m[ii]));
				}
				sigma_kg->setCoeff(0, 1, 0.5 * (sigma_k[3] + sigma_m[3]));
				sigma_kg->setCoeff(0, 2, 0.5 * (sigma_k[5] + sigma_m[5]));
				sigma_kg->setCoeff(1, 2, 0.5 * (sigma_k[4] + sigma_m[4]));
			}
			else if (numDim == 2)
			{
				for (int ii = 0; ii < 2; ii++)
				{
					sigma_kg->setCoeff(ii, ii, 0.5 * (sigma_k[ii] + sigma_m[ii]));
				}
				sigma_kg->setCoeff(0, 1, 0.5 * (sigma_k[3] + sigma_m[3]));
			}
			matoperat.dSymeEigenV('V', sigma_kg, sigma_pri, simga_priaxis);
			maxIDX = max_element(sigma_pri->cdp_vecCoeff, sigma_pri->cdp_vecCoeff + Matrc) - sigma_pri->cdp_vecCoeff;
			for (int ii = 0; ii < numDim; ii++)
			{
				eigV[ii] = simga_priaxis->d_getCoeff(ii, maxIDX);
			}
			temp =sqrt(eigV[0] * eigV[0] + eigV[1] * eigV[1] + eigV[2] * eigV[2]);
			for (int ii = 0; ii < numDim; ii++)
			{
				eigV[ii] = eigV[ii]/temp;
			}
			//=====judge bond===
			if (numDim==3)
			{
				updateBondState_CirclePlane(xc, eigV,maxDelta, o_dat);
			}
			else if (numDim == 2)
			{
				crack[0][0] = xc[0] - maxDelta * (eigV[1]);
				crack[0][1] = xc[1] + maxDelta * eigV[0];
				crack[1][0] = xc[0] + maxDelta * (eigV[1]);
				crack[1][1] = xc[1] - maxDelta * eigV[0];
				updateBondstate(crack, o_dat);
			}
		}
		delete sigma_kg, delete sigma_pri, delete simga_priaxis;
	}
	delete[] finaFamk, delete[]fina_m;
	if (finaSize>0)
	{
		return false;
	}
	else
	{
		return true;
	}
}



double pdsolve::failureCriterion_maxPriSig(datModel& o_dat, int Tk, bool& addLoad)
{
	//maximum principal stress failure criterion;
	int numDim = o_dat.ci_Numdimen;
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
	int Matrc = numDim;
	Matrix* sigma_kg = new Matrix(Matrc, Matrc);
	Matrix* simga_priaxis = new Matrix(Matrc, Matrc);// dummy here.
	Vector* sigma_pri = new Vector(Matrc);
	priority_queue<critStru> TopK;
	double topVal;
	for (int famkk = startP; famkk < endP; famkk++)
	{
		temP_fami = o_dat.op_getFami(famkk);
		numNodeOfFami = temP_fami->getNumNode();
		NID_k = temP_fami->getNodeID(0);
		o_dat.op_getNode(NID_k - 1)->getStress(sigma_k);
		for (int m = 1; m < numNodeOfFami; m++)
		{
			mu_km = temP_fami->getbondstate(m);
			if (mu_km==1)
			{
				NID_m = temP_fami->civ_NID[m];
				o_dat.op_getNode(NID_m - 1)->getStress(sigma_m);
				// symmetry matrix, store by upper triangle
				if (numDim==3)
				{
					for (int ii = 0; ii < 3; ii++)
					{
						sigma_kg->setCoeff(ii, ii, 0.5 * (sigma_k[ii] + sigma_m[ii]));
					}
					sigma_kg->setCoeff(0, 1, 0.5 * (sigma_k[3] + sigma_m[3]));
					sigma_kg->setCoeff(0, 2, 0.5 * (sigma_k[5] + sigma_m[5]));
					sigma_kg->setCoeff(1, 2, 0.5 * (sigma_k[4] + sigma_m[4]));
				}
				else if (numDim == 2)
				{
					for (int ii = 0; ii < 2; ii++)
					{
						sigma_kg->setCoeff(ii, ii, 0.5 * (sigma_k[ii] + sigma_m[ii]));
					}
					sigma_kg->setCoeff(0, 1, 0.5 * (sigma_k[3] + sigma_m[3]));
				}
				
				matoperat.dSymeEigenV('N', sigma_kg, sigma_pri, simga_priaxis);
				sigma1 = *max_element(sigma_pri->cdp_vecCoeff, sigma_pri->cdp_vecCoeff + Matrc);
				/*maxIDX = max_element(sigma_pri->cdp_vecCoeff, sigma_pri->cdp_vecCoeff + 3)- sigma_pri->cdp_vecCoeff;
				for (int ii = 0; ii < 3; ii++)
				{
					eigV[ii] = simga_priaxis->d_getCoeff(ii, maxIDX);
				}*/
				if (o_dat.ci_solvFlag != 1)//non-static solver
				{
					if (sigma1 > sigult)
					{
						critStru temStr(sigma1, famkk, m);
						if (Tk==0||int(TopK.size()) < Tk)
						{
							TopK.push(temStr);
						}
						else
						{
							topVal = TopK.top().sd_value;
							if (sigma1 > topVal)
							{
								TopK.pop();
								TopK.push(temStr);
							}
						}
					}
				}
				else//static solver
				{
					if (sigma1 > sigma1_max)
					{
						sigma1_max = sigma1;
					}
				}
			}
			
		}
	}
	delete sigma_kg, delete sigma_pri, delete simga_priaxis;

	double maxSig = 0;
	addLoad = true;
	if (o_dat.ci_solvFlag != 1 && numFami > 0)//non-static solver
	{
		addLoad = Finally_TopK_andUpdate_bondStaus(o_dat, Tk, TopK,2);//flag==2
	}
	else if (o_dat.ci_solvFlag == 1)//static solver;
	{
		MPI_Allreduce(&sigma1_max, &maxSig, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
	return maxSig;
}

bool pdsolve::segCirclePlanIntersect(double xp1[], double xp2[], double xc[], double norm[], double Del)
{
	double v1[3], v2[3];
	double temp1, temp2;
	for (int i = 0; i < 3; i++)
	{
		v1[i] = xp1[i] - xc[i];
		v2[i] = xp2[i] - xc[i];
	}
	temp1 = v1[0] * norm[0] + v1[1] * norm[1] + v1[2] * norm[2];
	temp2 = v2[0] * norm[0] + v2[1] * norm[1] + v2[2] * norm[2];
	if (temp1 * temp2 > 0)// if p1 and p2 are at the same side of the crack plane
	{
		return false;
	}
	else if (temp1 * temp2 == 0)
	{
		//printf("Warning: node on crack plane\n");
		return true;
	}
	else
	{
		double tempV1[3], tempV2[3];
		for (int i = 0; i < 3; i++)
		{
			tempV1[i] = xc[i] - xp1[i];
			tempV2[i] = xp2[i] - xp1[i];
		}
		double d = (tempV1[0] * norm[0] + tempV1[1] * norm[1] + tempV1[2] * norm[2]) /
			(tempV2[0] * norm[0] + tempV2[1] * norm[1] + tempV2[2] * norm[2]);
		double Xisc[3], R2;
		for (int i = 0; i < 3; i++)
		{
			Xisc[i] = xp1[i] + d * tempV2[i];// the intersecting point of p1p2 with the crack plane
		}
		R2 = (Xisc[0] - xc[0]) * (Xisc[0] - xc[0]) + (Xisc[1] - xc[1]) * (Xisc[1] - xc[1]) + (Xisc[2] - xc[2]) * (Xisc[2] - xc[2]);
		if (R2<0.25*Del*Del)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
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
	//if (ci_rank==0)
	{
		int NID_k;
		for (int famk = 0; famk < numFami; famk++)
		{
			NID_k = o_dat.op_getFami(famk)->getNodeID(0);
			o_dat.op_getNode(NID_k - 1)->setLocalDamage(glo_locaDama[famk]);
		}
	}
}

bool pdsolve::b_cracPropag_qusiaStatic_BYKeq(datModel& o_dat, ofstream& cracPath, double &KEQ)
{
	bool b_startPropage = false;
	int* ip_startPropage = new int[o_dat.ci_numCrack];
	for (int i = 0; i < o_dat.ci_numCrack; i++)
	{
		ip_startPropage[i] = 0;
	}
	//======
	double blocSize = o_dat.op_getGeomP()->getBlockSize();
	double lbc[3], rtc[3];
	o_dat.op_getGeomP()->getlbc(lbc);
	o_dat.op_getGeomP()->getrtc(rtc);
	//=====
	double KIc = o_dat.op_getmaterial()->getKIc();
	//vector<double*> ::iterator iter; //*iter, what store in the vector
	double min1stDisSqua, Xtip[3], XNode[3], disSqua, delt;
	double KI, KII, theta_c, R, Keq;
	int Idex[3], numNodeOfBloc, NID, NID_near1st = -1, famk, blockIndex, numBlocks[3];
	o_dat.getNumOfBLock(numBlocks);
	//===============set R ========================
	double m_R = o_dat.cd_mr;
	//===============================================
	//for (iter = cv_cracTips.begin(); iter != cv_cracTips.end(); iter++)
	for (int ck = 0; ck < o_dat.ci_numCrack; ck++)
	{

		//==find out the nearest PD node of crack tip;
		min1stDisSqua = 1.0e200;
		/*Xtip[0] = (*iter)[2];
		Xtip[1] = (*iter)[3];*/
		for (int i = 0; i < 3; i++)
		{
			Xtip[i] = o_dat.cdp_crack[ck][1][i];
			Idex[i] = (Xtip[i] - lbc[i]) / blocSize;
		}
		blockIndex = Idex[0] + Idex[1] * numBlocks[0] +
			Idex[2] * numBlocks[0] * numBlocks[1];
		numNodeOfBloc = o_dat.op_getBlock(blockIndex)->getNumNodeoB();
		for (int n = 0; n < numNodeOfBloc; n++)
		{
			NID = o_dat.op_getBlock(blockIndex)->getNodeoB(n);
			o_dat.op_getNode(NID - 1)->getcoor(XNode);
			disSqua = (XNode[0] - Xtip[0]) * (XNode[0] - Xtip[0]) + (XNode[1] - Xtip[1]) * (XNode[1] - Xtip[1]);
			if (disSqua <= min1stDisSqua)
			{
				min1stDisSqua = disSqua;
				NID_near1st = NID;
			}
		}
		famk = o_dat.op_getNode(NID_near1st - 1)->getFamID() - 1;
		//cout << famk << "================================" << endl;
		delt = o_dat.op_getFami(famk)->gethorizon() / o_dat.op_getGeomP()->getFactor();
		R = m_R * delt;

		// get SIFs and theta_c;
		SIFsAndPropaDire(KI, KII, theta_c, R, m_R, o_dat, ck, cracPath);
		Keq = cos(0.5 * theta_c) * (KI * cos(0.5 * theta_c) * cos(0.5 * theta_c) - 1.5 * KII * sin(theta_c));
		KEQ = Keq;
		// you may uncomment here to see value;
		if (ci_rank==0&&o_dat.ci_solvFlag==2)
		{
			printf("KI= %e, KII=%e, Keq=%e\n", KI, KII, Keq);
		}
		////==update new crack tip;
		if (Keq > KIc)
		{
			b_startPropage = true;
			ip_startPropage[ck] = 1;
			//double vec_tip[2] = { (*iter)[2] - (*iter)[0],(*iter)[3] - (*iter)[1] };
			double vec_tip[2] = { o_dat.cdp_crack[ck][1][0] - o_dat.cdp_crack[ck][0][0], o_dat.cdp_crack[ck][1][1] - o_dat.cdp_crack[ck][0][1] };
			double alpha = atan2(vec_tip[1], vec_tip[0]);
			double cracDire[2];
			cracDire[0] = cos(alpha) * cos(theta_c) - sin(alpha) * sin(theta_c);
			cracDire[1] = sin(alpha) * cos(theta_c) + cos(alpha) * sin(theta_c);
			double DeltaMin =o_dat.cd_dcf* o_dat.op_getGeomP()->getmaxDelta();
			//double DeltaMin = 1.0e-3;//d_c
			double xNewTip[2] = { Xtip[0] + cracDire[0] * DeltaMin, Xtip[1] + cracDire[1] * DeltaMin };
			o_dat.cdp_crack[ck][0][0] = Xtip[0];
			o_dat.cdp_crack[ck][0][1] = Xtip[1];
			o_dat.cdp_crack[ck][1][0] = xNewTip[0];
			o_dat.cdp_crack[ck][1][1] = xNewTip[1];
		}
	}

	//update the bond state;
	if (o_dat.ci_solvFlag==2)
	{
		for (int ck = 0; ck < o_dat.ci_numCrack; ck++)
		{
			if (ip_startPropage[ck] == 1)
			{
				// you may uncomment here to see value;
				if (ci_rank == 0 && o_dat.ci_solvFlag == 2)
				{
					printf("Xtip =( %e, %e),  XnewTip = (%e, %e) \n", o_dat.cdp_crack[ck][0][0], o_dat.cdp_crack[ck][0][1], o_dat.cdp_crack[ck][1][0], o_dat.cdp_crack[ck][1][1]);
				}
				updateBondstate(o_dat.cdp_crack[ck], o_dat);
			}
		}
	}

	// output file path to file
	for (int ck = 0; ck < o_dat.ci_numCrack; ck++)
	{
		cracPath << o_dat.cdp_crack[ck][1][0] << "\t" << o_dat.cdp_crack[ck][1][1] << "\t";
	}
	cracPath << endl;

	delete[] ip_startPropage;
	ip_startPropage = NULL;
	return b_startPropage;
}

void pdsolve::SIFsAndPropaDire(double& KI, double& KII, double& theta_c, double R, double m_R, datModel& o_dat, int ck, ofstream& fout)
{
	vector<int> v_ele;
	vector<int*>v_NoType_Jinte;
	double Xtip[3] = {o_dat.cdp_crack[ck][1][0], o_dat.cdp_crack[ck][1][1],o_dat.cdp_crack[ck][1][2] };
	//printf("tip=%f,%f\n", Xtip[0], Xtip[1]);
	int IDX[3], numEleoB, EleID, numNIn, NID[4], * nT_J;
	double blockSize = o_dat.op_getGeomP()->getBlockSize();
	double lbc[3], rtc[3], xN[4][3], dis;
	o_dat.op_getGeomP()->getlbc(lbc);
	o_dat.op_getGeomP()->getrtc(rtc);
	//=============find out the path for integral=========.
	// using block to search, but if R is too large this way doesn't work;
	if (m_R < o_dat.op_getGeomP()->getFactor() && m_R>0)
	{
		//R is small using block to search element;
		int numBlocks[3], blockIndex;
		o_dat.getNumOfBLock(numBlocks);
		for (int i = 0; i < 3; i++)
		{
			IDX[i] = (Xtip[i] - lbc[i]) / blockSize;
		}
		for (int idx = IDX[0] - 1; idx < IDX[0] + 2; idx++)
		{
			for (int idy = IDX[1] - 1; idy < IDX[1] + 2; idy++)
			{
				if (idx >= 0 && idx < (numBlocks[0]) && idy>=0 && idy < (numBlocks[1]))
				{
					blockIndex = idx + idy * numBlocks[0];
					numEleoB = o_dat.op_getBlock(blockIndex)->getNumEleoB();
					for (int k = 0; k < numEleoB; k++)
					{
						numNIn = 0;
						EleID = o_dat.op_getBlock(blockIndex)->getEleoB(k);

						o_dat.op_getEles(EleID - 1)->getConNid(NID);//NID [4], be caution;
						if (o_dat.op_getEles(EleID - 1)->ci_numNodes==3)
						{
							NID[3] = NID[2];
						}
						nT_J = new int[4];
						for (int n = 0; n < 4; n++)
						{
							nT_J[n] = 1;
						}
						for (int n = 0; n < 4; n++)
						{
							o_dat.op_getNode(NID[n] - 1)->getcoor(xN[n]);
							dis = (xN[n][0] - Xtip[0]) * (xN[n][0] - Xtip[0]) + (xN[n][1] - Xtip[1]) * (xN[n][1] - Xtip[1]);
							if (dis < R * R)
							{
								nT_J[n] = -1;
								numNIn++;
							}
						}
						if (numNIn == 0 || numNIn == 4 || numNIn == 3)
						{
							delete[] nT_J;
							nT_J = NULL;
						}
						else
						{
							v_NoType_Jinte.push_back(nT_J);
							v_ele.push_back(EleID);
						}
					}
				}
			}
		}
	}
	else
	{
		// if R is too large, search all PD element;
		for (int ele = 0; ele < o_dat.getTotnumEle(); ele++)
		{
			if (o_dat.op_getEles(ele)->getAlgoType() == 1)
			{
				//PD element;
				numNIn = 0;
				EleID = ele + 1;
				o_dat.op_getEles(ele)->getConNid(NID);//NID [4], be caution;
				if (o_dat.op_getEles(ele)->ci_numNodes == 3)
				{
					NID[3] = NID[2];
				}
				nT_J = new int[4];
				for (int n = 0; n < 4; n++)
				{
					nT_J[n] = 1;
				}
				for (int n = 0; n < 4; n++)
				{
					o_dat.op_getNode(NID[n] - 1)->getcoor(xN[n]);
					dis = (xN[n][0] - Xtip[0]) * (xN[n][0] - Xtip[0]) + (xN[n][1] - Xtip[1]) * (xN[n][1] - Xtip[1]);
					if (dis < R * R)
					{
						nT_J[n] = -1;
						numNIn++;
					}
				}
				if (numNIn == 0 || numNIn == 4 || numNIn == 3)
				{
					delete[] nT_J;
					nT_J = NULL;
				}
				else
				{
					v_NoType_Jinte.push_back(nT_J);
					v_ele.push_back(EleID);
				}
			}
		}
	}


	//======J_integral to get K_I and K_II====
	double I_m1 = 0, I_m2 = 0, temp_I1, temp_I2, d_J = 0, temp_J;
	int numPath = v_ele.size();
	int nNext;
	double vp[2], vt[2], normV[2], temp;
	for (int i = 0; i < numPath; i++)
	{
		EleID = v_ele[i];
		nT_J = v_NoType_Jinte[i];
		o_dat.op_getEles(EleID - 1)->getConNid(NID);//NID [4], be caution;
		if (o_dat.op_getEles(EleID - 1)->ci_numNodes == 3)
		{
			NID[3] = NID[2];
		}
		for (int n = 0; n < 4; n++)
		{
			o_dat.op_getNode(NID[n] - 1)->getcoor(xN[n]);
		}
		for (int n = 0; n < 4; n++)
		{
			if (nT_J[n] == 1)
			{
				if (n == 3)
				{
					nNext = 0;
				}
				else
				{
					nNext = n + 1;
				}
				if (nT_J[nNext] == 1)
				{
					//==get normal vector;
					//printf("eleID=%d, NID1=%d, NID2=%d ", EleID, NID[n], NID[nNext]);
					if (NID[n] != NID[nNext]) // not triangle PD element
					{
						for (int ii = 0; ii < 2; ii++)
						{
							vp[ii] = xN[nNext][ii] - xN[n][ii];
							vt[ii] = xN[nNext][ii] - Xtip[ii];
						}
						temp = vp[0] * vt[1] - vp[1] * vt[0];
						if (temp > 0)
						{
							normV[0] = -vp[1];  // this norm vector is not unit norm vector;
							normV[1] = vp[0];
						}
						else if (temp < 0)
						{
							normV[0] = vp[1];
							normV[1] = -vp[0];
						}
						else
						{
							cout << "Error: crack tip located on J integral path" << endl;
							exit(0);
						}
						//printf(" normal vector=%f,%f\n", normV[0], normV[1]);
						// calculate KI and KII;
						//Jintegrand(temp_I1, temp_I2, temp_J, NID[n], NID[nNext], normV, cracTip, o_dat, test);
						Jintegrand(temp_I1, temp_I2, temp_J, NID[n], NID[nNext], normV, ck, o_dat, fout);
						I_m1 = I_m1 + temp_I1;
						I_m2 = I_m2 + temp_I2;
						d_J = d_J + temp_J;
					}
				}
			}
		}
	}
	//==get theta_c===
	int proType = o_dat.getProType();
	double E = o_dat.op_getmaterial()->getE();
	double nu = o_dat.op_getmaterial()->getnu();
	double E_prim;
	if (proType == 1)
	{
		E_prim = E;  // plane stress;
	}
	else if (proType == 2)
	{
		E_prim = E / (1 - nu * nu);//plane strain;
	}
	KI = 0.5 * E_prim * I_m1;
	KII = 0.5 * E_prim * I_m2;
	//printf("KI= %e, KII=%e, J=%e, JJ=%e\n", KI, KII, d_J,(KI*KI+KII*KII)/E_prim);
	//=======get theta_c===============
	//		sometimes KI<0, this way doesn't work
	/*if (abs(KI)<1.0e-15)
	{
		if (KII>0)
		{
			theta_c = 2 * atan(-0.5 * sqrt(2));
		}
		else
		{
			theta_c = 2 * atan(0.5 * sqrt(2));
		}
	}
	else
	{
		temp = -2.0 * KII / KI / (1 + sqrt(1.0 + 8.0 * KII * KII / KI / KI));
		theta_c = 2.0 * atan(temp);
	}*/
	//====maybe this way is better.
	if (KII > 0)
	{
		theta_c = 2 * atan(0.25 * (KI / KII - sqrt(KI * KI / KII / KII + 8.0)));
	}
	else if (KII < 0)
	{
		theta_c = 2 * atan(0.25 * (KI / KII + sqrt(KI * KI / KII / KII + 8.0)));
	}
	else
	{
		theta_c = 0;
	}


	// release memory
	vector<int*>::iterator iter_nTJ;
	for (iter_nTJ = v_NoType_Jinte.begin(); iter_nTJ != v_NoType_Jinte.end(); iter_nTJ++)
	{
		delete[] * iter_nTJ;
		*iter_nTJ = NULL;
	}
	v_NoType_Jinte.clear();
	v_ele.clear();

}

void pdsolve::Jintegrand(double& I_m1, double& I_m2, double& J_s1, int NID1, int NID2, double norV[], int ck, datModel& o_dat, ofstream& fout)
{
	int proType = o_dat.getProType();
	double E = o_dat.op_getmaterial()->getE();
	double nu = o_dat.op_getmaterial()->getnu();
	double mu = 0.5 * E / (1 + nu);
	double kappa;
	if (proType == 1)
	{
		//plane stress;
		kappa = (3 - nu) / (1 + nu);// plane stress;
	}
	else if (proType == 2)
	{
		//plane strain;
		kappa = (3 - 4 * nu);
	}
	else
	{
		cout << "problem type is not 1 or 2" << endl;
		exit(0);
	}

	//sigma_glo is state 1, need to be gransfor to local coordinate;
	double w_m1[2], w_m2[2], sigma_glo[2][3], eps_m1[2][3], eps_m2[2][3];// be caution: this eps is not engineering strain;
	double sig_t[2][6];
	double sigma_loc[2][3];
	double sig_m1[2][3], sig_m2[2][3], sig_z;
	double r[2], theta[2];
	o_dat.op_getNode(NID1 - 1)->getStress(sig_t[0]);
	o_dat.op_getNode(NID2 - 1)->getStress(sig_t[1]);
	for (int i = 0; i < 2; i++)
	{
		sigma_glo[i][0] = sig_t[i][0];
		sigma_glo[i][1] = sig_t[i][1];
		sigma_glo[i][2] = sig_t[i][3];
	}
	//==find r and theta====
	double XoldTip[3], XnewTip[3], vt[3];
	/*XoldTip[0] = cracTip[0]; XoldTip[1] = cracTip[1];
	XnewTip[0] = cracTip[2]; XnewTip[1] = cracTip[3];*/
	for (int i = 0; i < 3; i++)
	{
		XoldTip[i] = o_dat.cdp_crack[ck][0][i];
		XnewTip[i] = o_dat.cdp_crack[ck][1][i];
		vt[i] = XnewTip[i] - XoldTip[i];
	}
	double xN[2][3], vp[3], alpha, vpT[3];
	alpha = atan2(vt[1], vt[0]);
	o_dat.op_getNode(NID1 - 1)->getcoor(xN[0]);
	o_dat.op_getNode(NID2 - 1)->getcoor(xN[1]);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			vp[j] = xN[i][j] - XnewTip[j];
		}
		r[i] = sqrt(vp[0] * vp[0] + vp[1] * vp[1]);
		vpT[0] = cos(alpha) * vp[0] + sin(alpha) * vp[1];
		vpT[1] = -sin(alpha) * vp[0] + cos(alpha) * vp[1];
		theta[i] = atan2(vpT[1], vpT[0]);
	}
	//====Trans  from global to local coordinae;
	// trans normal vector
	double n_loc[2];
	n_loc[0] = cos(alpha) * norV[0] + sin(alpha) * norV[1];
	n_loc[1] = -sin(alpha) * norV[0] + cos(alpha) * norV[1];
	//trans sigma of state 1;
	for (int i = 0; i < 2; i++)
	{
		sigma_loc[i][0] = (cos(alpha) * sigma_glo[i][0] + sin(alpha) * sigma_glo[i][2]) * cos(alpha) +
			(cos(alpha) * sigma_glo[i][2] + sin(alpha) * sigma_glo[i][1]) * sin(alpha);
		sigma_loc[i][1] = -(-sin(alpha) * sigma_glo[i][0] + cos(alpha) * sigma_glo[i][2]) * sin(alpha) +
			(-sin(alpha) * sigma_glo[i][2] + cos(alpha) * sigma_glo[i][1]) * cos(alpha);
		sigma_loc[i][2] = -(cos(alpha) * sigma_glo[i][0] + sin(alpha) * sigma_glo[i][2]) * sin(alpha) +
			(cos(alpha) * sigma_glo[i][2] + sin(alpha) * sigma_glo[i][1]) * cos(alpha);
	}
	//trans displacment gradient;
	double dg_glo[2][4], PUxPx_loc[2], PUyPx_loc[2];
	int NID[2] = { NID1,NID2 };
	for (int i = 0; i < 2; i++)
	{
		DispGrad(dg_glo[i], NID[i], o_dat);
		PUxPx_loc[i] = (cos(alpha) * dg_glo[i][0] + sin(alpha) * dg_glo[i][3]) * cos(alpha) + (cos(alpha) * dg_glo[i][2] + sin(alpha) * dg_glo[i][1]) * sin(alpha);
		PUyPx_loc[i] = (-sin(alpha) * dg_glo[i][0] + cos(alpha) * dg_glo[i][3]) * cos(alpha) + (-sin(alpha) * dg_glo[i][2] + cos(alpha) * dg_glo[i][1]) * sin(alpha);
	}
	//===J integral===
	double w_J[2], secItem_J[2], eps_loc[2][3];
	for (int i = 0; i < 2; i++)
	{
		if (proType == 1)
		{
			sig_z = 0;  // plane stress;
		}
		else if (proType == 2)
		{
			sig_z = nu * (sigma_loc[i][0] + sigma_loc[i][1]);  // plane strain;
		}
		eps_loc[i][0] = 1.0 / E * (sigma_loc[i][0] - nu * (sigma_loc[i][1] + sig_z));
		eps_loc[i][1] = 1.0 / E * (sigma_loc[i][1] - nu * (sigma_loc[i][0] + sig_z));
		eps_loc[i][2] = sigma_loc[i][2] * (1 + nu) / E;
	}
	for (int i = 0; i < 2; i++)
	{
		w_J[i] = 0.5 * (sigma_loc[i][0] * eps_loc[i][0] + sigma_loc[i][1] * eps_loc[i][1] + 2 * sigma_loc[i][2] * eps_loc[i][2]) * n_loc[0];
		secItem_J[i] = sigma_loc[i][0] * PUxPx_loc[i] * n_loc[0] + sigma_loc[i][2] * (PUxPx_loc[i] * n_loc[1] + PUyPx_loc[i] * n_loc[0]) + sigma_loc[i][1] * PUyPx_loc[i] * n_loc[1];
	}
	J_s1 = (w_J[0] + w_J[1]) * 0.5 - (secItem_J[0] + secItem_J[1]) * 0.5;
	//======================MODE I ====================================
	//== mode 1; interaction strain energy
	for (int i = 0; i < 2; i++)
	{
		sig_m1[i][0] = 1.0 / sqrt(2 * PI * r[i]) * cos(0.5 * theta[i]) * (1 - sin(0.5 * theta[i]) * sin(1.5 * theta[i]));
		sig_m1[i][1] = 1.0 / sqrt(2 * PI * r[i]) * cos(0.5 * theta[i]) * (1 + sin(0.5 * theta[i]) * sin(1.5 * theta[i]));
		sig_m1[i][2] = 1.0 / sqrt(2 * PI * r[i]) * sin(0.5 * theta[i]) * cos(0.5 * theta[i]) * cos(1.5 * theta[i]);
		if (proType == 1)
		{
			sig_z = 0;  // plane stress;
		}
		else if (proType == 2)
		{
			sig_z = nu * (sig_m1[i][0] + sig_m1[i][1]);  // plane strain;
		}
		eps_m1[i][0] = 1.0 / E * (sig_m1[i][0] - nu * (sig_m1[i][1] + sig_z));
		eps_m1[i][1] = 1.0 / E * (sig_m1[i][1] - nu * (sig_m1[i][0] + sig_z));
		eps_m1[i][2] = sig_m1[i][2] * (1 + nu) / E;
	}
	for (int i = 0; i < 2; i++)
	{
		w_m1[i] = (sigma_loc[i][0] * eps_m1[i][0] + sigma_loc[i][1] * eps_m1[i][1] + 2 * sigma_loc[i][2] * eps_m1[i][2]) * n_loc[0];
	}
	//===Mode 1, the second Item;
	double PuxPx_m1[2], PuyPx_m1[2], secItem[2];
	for (int i = 0; i < 2; i++)
	{
		PuxPx_m1[i] = (1 + nu) / E / sqrt(2 * PI * r[i]) * cos(0.5 * theta[i]) * ((kappa - 1) * 0.5 - sin(0.5 * theta[i]) * sin(1.5 * theta[i]));
		PuyPx_m1[i] = (1 + nu) / E / sqrt(2 * PI * r[i]) * sin(0.5 * theta[i]) * (-0.5 * (kappa + 1) + cos(0.5 * theta[i]) * cos(0.5 * theta[i]) - sin(theta[i]) * sin(theta[i]));
	}
	for (int i = 0; i < 2; i++)
	{
		secItem[i] = sigma_loc[i][0] * PuxPx_m1[i] * n_loc[0] + sigma_loc[i][2] * (PuxPx_m1[i] * n_loc[1] + PuyPx_m1[i] * n_loc[0]) + sigma_loc[i][1] * PuyPx_m1[i] * n_loc[1];
	}
	//===Mode 1, the Third Item;
	double ThirItem[2];
	for (int i = 0; i < 2; i++)
	{
		ThirItem[i] = sig_m1[i][0] * PUxPx_loc[i] * n_loc[0] + sig_m1[i][2] * (PUxPx_loc[i] * n_loc[1] + PUyPx_loc[i] * n_loc[0]) + sig_m1[i][1] * PUyPx_loc[i] * n_loc[1];
	}
	//== get I_m1==;
	I_m1 = (w_m1[0] + w_m1[1]) * 0.5 - (secItem[0] + secItem[1]) * 0.5 - (ThirItem[0] + ThirItem[1]) * 0.5;
	//================MODE II=======================================
	//== mode 2; interaction strain energy
	for (int i = 0; i < 2; i++)
	{
		sig_m2[i][0] = -1.0 / sqrt(2 * PI * r[i]) * sin(0.5 * theta[i]) * (2 + cos(0.5 * theta[i]) * cos(1.5 * theta[i]));
		sig_m2[i][1] = 1.0 / sqrt(2 * PI * r[i]) * sin(0.5 * theta[i]) * cos(0.5 * theta[i]) * cos(1.5 * theta[i]);
		sig_m2[i][2] = 1.0 / sqrt(2 * PI * r[i]) * cos(0.5 * theta[i]) * (1 - sin(0.5 * theta[i]) * sin(1.5 * theta[i]));
		if (proType == 1)
		{
			sig_z = 0;  // plane stress;
		}
		else if (proType == 2)
		{
			sig_z = nu * (sig_m1[i][0] + sig_m1[i][1]);  // plane strain;
		}
		eps_m2[i][0] = 1.0 / E * (sig_m2[i][0] - nu * (sig_m2[i][1] + sig_z));
		eps_m2[i][1] = 1.0 / E * (sig_m2[i][1] - nu * (sig_m2[i][0] + sig_z));
		eps_m2[i][2] = sig_m2[i][2] * (1 + nu) / E;
	}
	for (int i = 0; i < 2; i++)
	{
		w_m2[i] = (sigma_loc[i][0] * eps_m2[i][0] + sigma_loc[i][1] * eps_m2[i][1] + 2 * sigma_loc[i][2] * eps_m2[i][2]) * n_loc[0];
	}
	//===Mode 2, the second Item;
	double PuxPx_m2[2], PuyPx_m2[2];
	for (int i = 0; i < 2; i++)
	{
		PuxPx_m2[i] = (1 + nu) / E / sqrt(2 * PI * r[i]) * sin(0.5 * theta[i]) * (sin(theta[i]) * sin(theta[i]) - cos(0.5 * theta[i]) * cos(0.5 * theta[i]) - 0.5 * (kappa + 1));
		PuyPx_m2[i] = (1 + nu) / E / sqrt(2 * PI * r[i]) * cos(0.5 * theta[i]) * (-sin(0.5 * theta[i]) * sin(1.5 * theta[i]) - 0.5 * (kappa - 1));
	}
	for (int i = 0; i < 2; i++)
	{
		secItem[i] = sigma_loc[i][0] * PuxPx_m2[i] * n_loc[0] + sigma_loc[i][2] * (PuxPx_m2[i] * n_loc[1] + PuyPx_m2[i] * n_loc[0]) + sigma_loc[i][1] * PuyPx_m2[i] * n_loc[1];
	}
	//===Mode 2, the Third Item;
	for (int i = 0; i < 2; i++)
	{
		ThirItem[i] = sig_m2[i][0] * PUxPx_loc[i] * n_loc[0] + sig_m2[i][2] * (PUxPx_loc[i] * n_loc[1] + PUyPx_loc[i] * n_loc[0]) + sig_m2[i][1] * PUyPx_loc[i] * n_loc[1];
	}
	//== get I_m2==;
	I_m2 = (w_m2[0] + w_m2[1]) * 0.5 - (secItem[0] + secItem[1]) * 0.5 - (ThirItem[0] + ThirItem[1]) * 0.5;
}

void pdsolve::updateBondstate(double crack[][3], datModel & o_dat)
{
	//update bond state;
	
	//================================
	//==crack 
	int numDime = o_dat.ci_Numdimen;
	double Xmin[3], Xmax[3], xk[3], xm[3];
	//===block and family
	double blockSize, lbc[3], temp[3];
	blockSize = o_dat.op_getGeomP()->getBlockSize();
	o_dat.op_getGeomP()->getlbc(lbc);
	int i_bIndexMin[3], i_bIndexMax[3], numBlocks[3], blockIndex, numNodeoB, NodeIDoB;
	int famID, numNodeFami, NID_k, NID_m;
	int fn = 0;
	if (numDime == 3)
	{
		fn = 3;
	}
	else if (numDime == 2)
	{
		fn = 2;
	}
	bool crosCrack;
	o_dat.getNumOfBLock(numBlocks);
	
	for (int kk = 0; kk < 3; kk++)
	{
		temp[0] = crack[0][kk];
		temp[1] = crack[1][kk];
		temp[2] = crack[2][kk];
		Xmin[kk] = *min_element(temp, temp + fn);
		Xmax[kk] = *max_element(temp, temp + fn);
	}
	for (int i = 0; i < 3; i++)
	{
		i_bIndexMin[i] = (Xmin[i] - lbc[i]) / blockSize;
		i_bIndexMax[i] = (Xmax[i] - lbc[i]) / blockSize;
	}
	if (numDime==2)
	{
		i_bIndexMin[2] = 0;
		i_bIndexMax[2] = 0;
	}
	if (ci_rank==0)
	{

	}
	for (int xIdex = i_bIndexMin[0] - 1; xIdex < i_bIndexMax[0] + 2; xIdex++)
	{
		for (int yIdex = i_bIndexMin[1] - 1; yIdex < i_bIndexMax[1] + 2; yIdex++)
		{
			for (int zIdex = i_bIndexMin[2] - 1; zIdex < i_bIndexMax[2] + 2; zIdex++)
			{
				if (xIdex >= 0 && xIdex < (numBlocks[0]) &&yIdex >= 0 && yIdex < (numBlocks[1]) &&zIdex >= 0 && zIdex < (numBlocks[2]))	
				{
					blockIndex = xIdex + yIdex * numBlocks[0] + zIdex * numBlocks[0] * numBlocks[1];
					numNodeoB = o_dat.op_getBlock(blockIndex)->getNumNodeoB();
					for (int ii = 0; ii < numNodeoB; ii++)
					{
						NodeIDoB = o_dat.op_getBlock(blockIndex)->getNodeoB(ii);
						famID = o_dat.op_getNode(NodeIDoB - 1)->getFamID();
						if (famID==-1)
						{
							continue;
						}
						if (o_dat.op_getFami(famID - 1)->cb_allowFail)
						{
							numNodeFami = o_dat.op_getFami(famID - 1)->getNumNode();
							NID_k = o_dat.op_getFami(famID - 1)->getNodeID(0);
							o_dat.op_getNode(NID_k - 1)->getcoor(xk);
							for (int m = 1; m < numNodeFami; m++)
							{
								NID_m = o_dat.op_getFami(famID - 1)->getNodeID(m);
								o_dat.op_getNode(NID_m - 1)->getcoor(xm);
								if (o_dat.op_getFami(famID - 1)->cip_bondState[m] == 1)
								{
									if (numDime == 3)
									{
										crosCrack = segmentPlaneIntersection(xk, xm, crack);
									}
									else if (numDime == 2)
									{
										crosCrack = intersection(xk, xm, crack[0], crack[1]);
									}
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

void pdsolve::updateBondState_CirclePlane(double xc[], double norm[],double maxDelta, datModel& o_dat)
{
	bool crosCrack;
	int i_bIndex[3], numBlocks[3], blockIndex, numNodeoB, NodeIDoB, famID, numNodeFami, NID_k, NID_m;
	double blockSize, lbc[3], xk[3], xm[3];
	blockSize = o_dat.op_getGeomP()->getBlockSize();
	o_dat.op_getGeomP()->getlbc(lbc);
	o_dat.getNumOfBLock(numBlocks);
	for (int i = 0; i < 3; i++)
	{
		i_bIndex[i] = (xc[i] - lbc[i]) / blockSize;
	}
	for (int xIdex = i_bIndex[0] - 1; xIdex < i_bIndex[0] + 2; xIdex++)
	{
		for (int yIdex = i_bIndex[1] - 1; yIdex < i_bIndex[1] + 2; yIdex++)
		{
			for (int zIdex = i_bIndex[2] - 1; zIdex < i_bIndex[2] + 2; zIdex++)
			{
				if (xIdex >= 0 && xIdex < (numBlocks[0]) && yIdex >= 0 && yIdex < (numBlocks[1]) && zIdex >= 0 && zIdex < (numBlocks[2]))
				{
					blockIndex = xIdex + yIdex * numBlocks[0] + zIdex * numBlocks[0] * numBlocks[1];
					numNodeoB = o_dat.op_getBlock(blockIndex)->getNumNodeoB();
					for (int ii = 0; ii < numNodeoB; ii++)
					{
						NodeIDoB = o_dat.op_getBlock(blockIndex)->getNodeoB(ii);
						famID = o_dat.op_getNode(NodeIDoB - 1)->getFamID();
						if (o_dat.op_getFami(famID - 1)->cb_allowFail)
						{
							numNodeFami = o_dat.op_getFami(famID - 1)->getNumNode();
							NID_k = o_dat.op_getFami(famID - 1)->getNodeID(0);
							o_dat.op_getNode(NID_k - 1)->getcoor(xk);
							for (int m = 1; m < numNodeFami; m++)
							{
								if (o_dat.op_getFami(famID - 1)->cip_bondState[m] == 1)
								{
									NID_m = o_dat.op_getFami(famID - 1)->getNodeID(m);
									o_dat.op_getNode(NID_m - 1)->getcoor(xm);
									crosCrack = segCirclePlanIntersect(xk, xm, xc, norm, maxDelta);
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

bool pdsolve::intersection(double L1X1[], double L1X2[], double L2X1[], double L2X2[])
{
	// Cheack two segments intersect with each other or not?
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
		int numNodeEle, * conNID, ele;
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
		int numNodeEle, * conNID, ele;
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
		int numNodeEle, * conNID, ele;
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
		int numNodeEle, * conNID, ele;
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
		if (numDime==2&&o_dat.ci_proType==2)
		{
			double nu = o_dat.op_getmaterial()->getnu();
			for (int i = 0; i < totNumNODE; i++)
			{
				o_dat.op_getNode(i)->calSigzz(nu);
			}
		}
	}
	delete[] count, delete[] gloCount;
	count = NULL, gloCount = NULL;
}

void pdsolve::postProcessing(datModel & o_dat, ofstream & test)
{
	//calGlobalNodeStresses(o_dat, test);
}
