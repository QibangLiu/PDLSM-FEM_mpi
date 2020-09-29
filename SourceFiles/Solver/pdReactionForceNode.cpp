#include "pdReactionForceNode.h"

pdReactionForceNode::pdReactionForceNode(int Nid)
{
	ci_NID = Nid;
	cd_Fp = 0;
	ci_numEle = 0;
}

void pdReactionForceNode::putEleInside(int eleID,int Index)
{
	civ_eleID.push_back(eleID);
	civ_NIndexOfEle.push_back(Index);
	ci_numEle = ci_numEle + 1;
}

int pdReactionForceNode::getNid() const
{
	return ci_NID;
}

int pdReactionForceNode::getnumELE() const
{
	return ci_numEle;
}

int pdReactionForceNode::getEleID(int i)
{
	return	civ_eleID[i];
}

int pdReactionForceNode::getNIndexOfEle(int i)
{
	return civ_NIndexOfEle[i];
}
