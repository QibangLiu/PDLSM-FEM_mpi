#include "pdGeomP.h"

pdGeomP::pdGeomP( double factor)
{
	for (int i = 0; i < 3; i++)
	{
		cd_lbc[i] = 0;
		cd_rtc[i] = 0;
	}
	cd_factor = factor;
	cd_maxDelta = 0;
	cd_minDelta = 0;
	
}

pdGeomP::~pdGeomP()
{
}

void pdGeomP::setRTC(double rtc[])
{
	for (int i = 0; i < 3; i++)
	{
		cd_rtc[i] = rtc[i];
	}
}

void pdGeomP::setLBC(double lbc[])
{
	for (int i = 0; i < 3; i++)
	{
		cd_lbc[i] = lbc[i];
	}
}

void pdGeomP::print(ofstream & fout)
{
	fout << "Diagonal point: (" << cd_lbc[0]<<", "<< cd_lbc[1]<< ")" << endl;
	fout << "         (" << cd_rtc[0] << ", " << cd_rtc[1] << ")" << endl;
	
	//fout << "volume of material point=" << cd_Volume << endl;
	
}

void pdGeomP::getrtc(double rtc[])
{
	rtc[0] = cd_rtc[0];
	rtc[1] = cd_rtc[1];
	rtc[2] = cd_rtc[2];
}

void pdGeomP::getlbc(double lbc[])
{
	lbc[0] = cd_lbc[0];
	lbc[1] = cd_lbc[1];
	lbc[2] = cd_lbc[2];
}

double pdGeomP::getFactor() const
{
	return cd_factor;
}

double pdGeomP::getmaxDelta() const
{
	return cd_maxDelta;
}

double pdGeomP::getminDelta() const
{
	return cd_minDelta;
}

double pdGeomP::getBlockSize() const
{
	return cd_blockSize;
}

void pdGeomP::setBlockSize(double val)
{
	cd_blockSize = val;
}

void pdGeomP::setMaxDelta(double del)
{
	cd_maxDelta = del;
}

void pdGeomP::setMinDelta(double del)
{
	cd_minDelta = del;
}
