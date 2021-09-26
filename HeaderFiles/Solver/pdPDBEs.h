#pragma once
class pdPDBEs
{
public:
	pdPDBEs(int nodeID[], int numNode);
	~pdPDBEs();
	void getNodeID(int nID[]);
	int i_getNumnode();
private:
	pdPDBEs();
	int *ci_nodeID;
	int /*ci_normDire,*/ ci_numNode;
};

