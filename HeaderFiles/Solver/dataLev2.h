#pragma once
/*This class is named second level data, is for storing some data in contigous array (1D array)*/
using namespace std;
class dataLev2
{
public:
	dataLev2();
	~dataLev2();
	void getNodeCOOR(int NIDX, double x[]);
	void setNodeCOOR(int NIDX, double x[]);
	//====node data====
	double* cdp_sigma;
	double* cdp_X;
};

