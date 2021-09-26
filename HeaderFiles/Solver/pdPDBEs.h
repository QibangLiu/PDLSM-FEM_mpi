#pragma once
class pdPDBEs
{
public:
	pdPDBEs(int nodeID[], int numNode);
	~pdPDBEs();
	void getNodeID(int nID[]);
	int i_getNumnode();
	int* ci_nodeID;
private:
	pdPDBEs();
	
	int /*ci_normDire,*/ ci_numNode;
};

