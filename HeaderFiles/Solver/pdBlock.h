#pragma once
#include<vector>
#include<fstream>
#include<iomanip>
using namespace std;

class pdBlock
{
public:
	pdBlock(int id);
	~pdBlock();
	void putNodeInBlock(int id);//put the id of element into block
	void putEleInBlock(int eleID);
	int getNodeoB(int i);//get the id of i_th element
	int getNumNodeoB()const;
	int getEleoB(int i);
	int getNumEleoB()const;
	void print(ofstream &fout);
	vector<int>cv_Nodeid;//for store the MP id in this block;
	vector<int>cv_eleID;
private:
	pdBlock();
	int ci_ID; 	
};
