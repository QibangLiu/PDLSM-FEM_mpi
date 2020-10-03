#include "dataLev2.h"

dataLev2::dataLev2()
{
	cdp_sigma = nullptr;
}

dataLev2::~dataLev2()
{
	if (cdp_sigma)
	{
		delete[]cdp_sigma;
		cdp_sigma = nullptr;
	}
}
