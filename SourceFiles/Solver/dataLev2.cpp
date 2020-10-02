#include "dataLev2.h"

dataLev2::dataLev2()
{
	cd_X = nullptr;
}

dataLev2::~dataLev2()
{
	if (cd_X)
	{
		delete[] cd_X;
		cd_X = nullptr;
	}
}
