#include "dataLev2.h"

dataLev2::dataLev2()
{
	cdp_sigma = nullptr;
	cdp_X = nullptr;
}

dataLev2::~dataLev2()
{
	if (cdp_sigma)
	{
		delete[]cdp_sigma;
		cdp_sigma = nullptr;
	}
	if (cdp_X)
	{
		delete[]cdp_X;
		cdp_X = nullptr;
	}
}

void dataLev2::getNodeCOOR(int NIDX, double x[])
{
	//NIDX-- node index;
	for (int i = 0; i < 3; i++)
	{
		x[i] = cdp_X[3 * NIDX + i];
	}
}

void dataLev2::setNodeCOOR(int NIDX, double x[])
{
	//NIDX-- node index;
	for (int i = 0; i < 3; i++)
	{
		 cdp_X[3 * NIDX + i]=x[i];
	}
}
