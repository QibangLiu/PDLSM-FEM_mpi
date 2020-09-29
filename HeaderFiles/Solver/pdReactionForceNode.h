#pragma once
#include<vector>
using namespace std;
class pdReactionForceNode
{
public:
	pdReactionForceNode(int Nid);
	void putEleInside(int eleID, int Index);
	int getNid()const;
	int getnumELE()const;
	int getEleID(int i);
	int getNIndexOfEle(int i);
private:
	pdReactionForceNode();
	vector<int> civ_eleID;
	vector<int> civ_NIndexOfEle;
	int ci_numEle;
	int ci_NID;
	double cd_Fp; //point force;
};

