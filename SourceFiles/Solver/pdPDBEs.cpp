#include "pdPDBEs.h"

pdPDBEs::pdPDBEs(int nodeID[], int numNode)
{
	ci_numNode = numNode;
	ci_nodeID = new int[ci_numNode];
	for (int i = 0; i < ci_numNode; i++)
	{
		ci_nodeID[i] = nodeID[i];
	}
}

pdPDBEs::~pdPDBEs()
{
	delete[] ci_nodeID;
	ci_nodeID = nullptr;
}

void pdPDBEs::getNodeID(int nID[])
{
	for (int i = 0; i < ci_numNode; i++)
	{
		nID[i] = ci_nodeID[i];
	}

}

int pdPDBEs::i_getNumnode()
{
	return ci_numNode;
}
