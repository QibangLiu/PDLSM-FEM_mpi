#pragma once
//Geometry parameters
#include<fstream>
using namespace std;

class pdGeomP
{
public:
	
	pdGeomP(double factor);
	~pdGeomP();
	void setRTC(double rtc[]);
	void setLBC(double lbc[]);
	void print(ofstream &fout);
	void getrtc(double rtc[]);
	void getlbc(double lbc[]);
	double getFactor()const;
	double getmaxDelta()const;
	double getminDelta()const;
	double getBlockSize()const;
	void setBlockSize(double val);	
	void setMaxDelta(double del);
	void setMinDelta(double del);
	//===
	double cd_factor;//horizon size factor
private:
	pdGeomP();
	double cd_rtc[3];//right top cornor coordinate;
	double cd_lbc[3];//left bottom cornor coordinate
	double cd_maxDelta, cd_minDelta;//the max,min delta of PD node;
	double cd_blockSize;
};

