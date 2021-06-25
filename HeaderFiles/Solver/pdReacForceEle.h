#pragma once
using namespace std;
#include<vector>
class pdReacForceEle
{
public:
	pdReacForceEle(int eleIDX);
	int ci_eleIDX;
	vector<int> civ_reaForceEleNidx;
	vector<int> civ_reaForceCPNDnodeIDX;// idx of civ_reaForceNID, corresponding to the idx of node of element;
};

