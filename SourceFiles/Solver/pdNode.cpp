#include "pdNode.h"


pdNode::pdNode(int id, dataLev2 *op_datLev2)
{
	ci_NodeId = id;
	cop_datLev2 = op_datLev2;
	for (int i = 0; i < 3; i++)
	{
		cop_dof[i] = new pdDof();
	}
	cd_dv = 0;
	ci_nodeType = 0;
	ci_famID = -1;
	for (int i = 0; i < 6; i++)
	{
		cd_sigma[i] = 0;
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
	double *cd_x = &(cop_datLev2->cd_X[3 * (ci_NodeId - 1)]);
	for (int i = 0; i < 3; i++)
	{
		x[i] = cd_x[i];
	}
}


double pdNode::getvolume() const
{
	return cd_dv;
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
	double* cd_x = &(cop_datLev2->cd_X[3 * (ci_NodeId - 1)]);
	fout << std::left << setiosflags(ios::scientific)
		<< setprecision(5);
	fout<< ci_NodeId << "\t" << cd_x[0]
	<< "\t" << cd_x[1] << "\t" << cd_x[2] << endl;
}

void pdNode::printStress(ofstream & fout)
{
	double* cd_x = &(cop_datLev2->cd_X[3 * (ci_NodeId - 1)]);
	fout << ci_NodeId << "\t";
	fout <<  cd_x[0] << "\t" << cd_x[1] << "\t" << cd_x[2];
	for (int i = 0; i < 6; i++)
	{
		fout << "\t" << cd_sigma[i];
	}
	fout << endl;
}

void pdNode::printStressTensor_vtk(ofstream& fout)
{
	fout << cd_sigma[0] << ' ' << cd_sigma[3] << ' ' << cd_sigma[5] << endl;
	fout << cd_sigma[3] << ' ' << cd_sigma[1] << ' ' << cd_sigma[4] << endl;
	fout << cd_sigma[5] << ' ' << cd_sigma[4] << ' ' << cd_sigma[2] << endl;
}



void pdNode::printFinalU(ofstream & PDout)
{
	double* cd_x = &(cop_datLev2->cd_X[3 * (ci_NodeId - 1)]);
	PDout << std::left << setiosflags(ios::scientific)<< setprecision(4);
	PDout<< ci_NodeId << "\t" << cd_x[0] << "\t" << cd_x[1]<<
		"\t" << cd_x[2] << "\t";
	//PDout << setprecision(8);
	PDout << cop_dof[0]->d_getValue() << "\t"
		<< cop_dof[1]->d_getValue() << "\t" << cop_dof[2]->d_getValue() << endl;
}

void pdNode::printDamage(ofstream & fout)
{
	double* cd_x = &(cop_datLev2->cd_X[3 * (ci_NodeId - 1)]);
	fout << ci_NodeId;
	fout << "\t" << cd_x[0] << "\t" << cd_x[1] << "\t" << cd_localDamage << endl;
}

void pdNode::setLocalDamage(double phi)
{
	cd_localDamage = phi;
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



void pdNode::setStress(int i, double val)
{
	cd_sigma[i] = val;
}

void pdNode::getStress(double sigma[])
{
	for (int i = 0; i < 6; i++)
	{
		sigma[i] = cd_sigma[i];
	}
}

void pdNode::addStress(int i, double val)
{
	cd_sigma[i] = cd_sigma[i] + val;
}

void pdNode::calAverageStress(int count)
{
	for (int i = 0; i < 6; i++)
	{
		cd_sigma[i] = cd_sigma[i] / count;
	}
}
