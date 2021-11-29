#pragma once

#include<iomanip>
#include<vector>
#include<fstream>
using namespace std;
class pdFamily
{
public:
	pdFamily();
	~pdFamily();
	void setID(int ID);
	void putNodeIntoFami(int NID);
	void printNODE(ofstream&fout);
	int getNodeID(int i);
	int getNumNode()const;
	double gethorizon()const;
	void sethorizon(double val);
	void initialBondState();
	int getbondstate(int m);
	void setbondstate(int m, int val);
	vector<int> GetvecNid();

public:
	vector<int>civ_NID;
	int* cip_bondState;//bond state ;
	bool cb_allowFail;
	double cd_v;
private: 
	double cd_delta;
	int ci_famID;//no need any more;
	
	
};
