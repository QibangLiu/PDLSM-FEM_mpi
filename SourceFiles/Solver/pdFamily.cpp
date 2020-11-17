#include "pdFamily.h"

pdFamily::pdFamily()
{
	cd_delta = 0;
	ci_famID = 0;
	cb_allowFail = true;
}

pdFamily::~pdFamily()
{
	delete[] cip_bondState;
	cip_bondState = NULL;
}

void pdFamily::setID(int ID)
{
	ci_famID = ID;
}

void pdFamily::putNodeIntoFami(int NID)
{
	civ_NID.push_back(NID);
}

void pdFamily::printNODE(ofstream & fout)
{
	fout <<ci_famID<<"** "<< civ_NID.size() << "**\t";
	for (int i = 0; i < civ_NID.size(); i++)
	{
		fout <<civ_NID[i]<<' ';
	}
	fout << endl;
}


int pdFamily::getNodeID(int i)
{
	return civ_NID[i];
}

int pdFamily::getNumNode() const
{
	return civ_NID.size();
}

double pdFamily::gethorizon() const
{
	return cd_delta;
}

void pdFamily::sethorizon(double val)
{
	cd_delta = val;
}

void pdFamily::initialBondState()
{
	int numNode = civ_NID.size();
	cip_bondState = new int[numNode];
	for (int i = 0; i < numNode; i++)
	{
		cip_bondState[i] = 1;
	}
}

int pdFamily::getbondstate(int m)
{
	return cip_bondState[m];
}

void pdFamily::setbondstate(int m, int val)
{
	cip_bondState[m] = val;
}

vector<int> pdFamily::GetvecNid()
{
	return vector<int>(civ_NID);
}





