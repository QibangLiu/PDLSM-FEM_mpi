#include "pdNode.h"


pdNode::pdNode(int id, double x[],dataLev2*p_datLev2)
{
	cop_datLev2 = p_datLev2;
	ci_NodeId = id;
	cop_datLev2->setNodeCOOR(id - 1, x);
	for (int i = 0; i < 3; i++)
	{
		cop_dof[i] = new pdDof();
	}
	cd_dv = 0;
	ci_nodeType = 0;
	ci_famID = -1;
	//cdp_sigma = &(cop_datLev2->cdp_sigma[6 * (id - 1)]);
	for (int i = 0; i < 6; i++)
	{
		(cop_datLev2->cdp_sigma[6 * (id - 1) + i]) = 0;
	}
	cd_localDamage = 0;
}

pdNode::~pdNode()
{
	for (int i = 0; i < 3; i++)
	{
		delete cop_dof[i];
		cop_dof[i] = NULL;
	}
}

int pdNode::getId() const
{
	return ci_NodeId;
}

void pdNode::getcoor(double x[]) const
{
	/*for (int i = 0; i < 3; i++)
	{
		x[i] = cd_x[i];
	}*/
	cop_datLev2->getNodeCOOR(ci_NodeId - 1, x);
}


double pdNode::getvolume() const
{
	return cd_dv;
}

void pdNode::setVolume(double dv)
{
	cd_dv = dv;
}

void pdNode::addvolume(double dv)
{
	cd_dv = cd_dv + dv;
}


pdDof * pdNode::op_getDof(int i) const
{
	return cop_dof[i];
}

void pdNode::print(ofstream & fout)
{
	double xN[3];
	cop_datLev2->getNodeCOOR(ci_NodeId - 1, xN);
	fout << std::left << setiosflags(ios::scientific)
		<< setprecision(5);
	fout<< ci_NodeId << "\t" << xN[0]
	<< "\t" << xN[1] << "\t" << xN[2] << endl;
}

void pdNode::printStress(ofstream & fout)
{
	double xN[3];
	cop_datLev2->getNodeCOOR(ci_NodeId - 1, xN);
	fout << ci_NodeId << "\t";
	fout <<  xN[0] << "\t" << xN[1] << "\t" << xN[2];
	for (int i = 0; i < 6; i++)
	{
		fout << "\t" << (cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]);
	}
	fout << endl;
}

void pdNode::printStressTensor_vtk(ofstream& fout)
{
	float cdp_sigma[6];
	for (int i = 0; i < 6; i++)
	{
		cdp_sigma[i] = (cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]);
	}
	fout << cdp_sigma[0] << ' ' << cdp_sigma[3] << ' ' << cdp_sigma[5] << endl;
	fout << cdp_sigma[3] << ' ' << cdp_sigma[1] << ' ' << cdp_sigma[4] << endl;
	fout << cdp_sigma[5] << ' ' << cdp_sigma[4] << ' ' << cdp_sigma[2] << endl;
	//0-xx, 1-yy, 2-zz, 3-xy, 4-yz, 5-zx;
}



void pdNode::printFinalU(ofstream & PDout)
{
	double xN[3];
	cop_datLev2->getNodeCOOR(ci_NodeId - 1, xN);
	PDout << std::left << setiosflags(ios::scientific)<< setprecision(4);
	PDout<< ci_NodeId << "\t" << xN[0] << "\t" << xN[1]<<
		"\t" << xN[2] << "\t";
	//PDout << setprecision(8);
	PDout << cop_dof[0]->d_getValue() << "\t"
		<< cop_dof[1]->d_getValue() << "\t" << cop_dof[2]->d_getValue() << endl;
}

void pdNode::printDamage(ofstream & fout)
{
	double cd_x[3];
	cop_datLev2->getNodeCOOR(ci_NodeId - 1, cd_x);
	fout << ci_NodeId;
	fout << "\t" << cd_x[0] << "\t" << cd_x[1] << "\t" << cd_localDamage << endl;
}

void pdNode::setLocalDamage(double phi)
{
	cd_localDamage = phi;
}

double pdNode::getLocalDamage()
{
	return cd_localDamage;
}

void pdNode::setNodeType(int val)
{
	ci_nodeType = val;
}

int pdNode::getNodeType()
{
	return ci_nodeType;
}


void pdNode::setFamID(int famID)
{
	ci_famID = famID;
}

int pdNode::getFamID() const
{
	return ci_famID;
}



void pdNode::calSigzz(const double& nu)
{
	//calculate sigma_zz for plane strain problem;
	cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + 2] = nu * (cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + 0] +
		cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + 1]);

}

void pdNode::setStress(int i, double val)
{
	(cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]) = val;
}

void pdNode::getStress(double sigma[])
{
	for (int i = 0; i < 6; i++)
	{
		sigma[i] = (cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]);
	}
}

void pdNode::addStress(int i, double val)
{
	(cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]) =
		(cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]) + val;
}

void pdNode::calAverageStress(int count)
{
	for (int i = 0; i < 6; i++)
	{
		(cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]) =
			(cop_datLev2->cdp_sigma[6 * (ci_NodeId - 1) + i]) / count;
	}
}
