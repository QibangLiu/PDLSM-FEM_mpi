#include "pdGaussPt.h"

pdGaussPt::pdGaussPt(int numGP)
{
	/*ci_numGp = 1;
	cd_point[0] = 0.0;
	cd_weight[0] = 2.0;*/

	ci_numGp = numGP;
	cdp_point = new double[numGP];
	cdp_weight = new double[numGP];
	switch (numGP)
	{
	case 1:
	{
		cdp_point[0] = 0;
		cdp_weight[0] = 2.0;
	}
	case 2:
	{
		cdp_point[0] = -1 / sqrt(3.0);
		cdp_point[1] = 1 / sqrt(3.0);
		cdp_weight[0] = 1.0;
		cdp_weight[1] = 1.0;
		break;
	}
	case 3:
	{
		cdp_point[0] = -sqrt(0.6);
		cdp_point[1] = 0;
		cdp_point[2] = sqrt(0.6);
		cdp_weight[0] = 5.0 / 9;
		cdp_weight[1] = 8.0 / 9;
		cdp_weight[2] = 5.0 / 9;
		break;
	}
	case 4:
	{
		cdp_point[0] = -sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5));
		cdp_point[1] = -sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5));
		cdp_point[2] = sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5));
		cdp_point[3] = sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5));

		cdp_weight[0] = 0.5 - sqrt(30.0) / 36.0;
		cdp_weight[1] = 0.5 + sqrt(30.0) / 36.0;
		cdp_weight[2] = 0.5 + sqrt(30.0) / 36.0;
		cdp_weight[3] = 0.5 - sqrt(30.0) / 36.0;
		break;
	}
	case 5:
	{
		cdp_point[0] = -1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10 / 7));
		cdp_point[1] = -1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10 / 7));
		cdp_point[2] = 0.0;
		cdp_point[3] = 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10 / 7));
		cdp_point[4] = 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10 / 7));

		cdp_weight[0] = (322.0 - 13 * sqrt(70.0)) / 900.0;
		cdp_weight[1] = (322.0 + 13 * sqrt(70.0)) / 900.0;
		cdp_weight[2] = 128.0 / 225.0;
		cdp_weight[3] = (322.0 + 13 * sqrt(70.0)) / 900.0;
		cdp_weight[4] = (322.0 - 13 * sqrt(70.0)) / 900.0;
		break;
	}
	default:
		break;
	}
	//if (numGP==2)
	//{
	//	cdp_point[0] = -1 / sqrt(3.0);
	//	cdp_point[1] = 1 / sqrt(3.0);
	//	cdp_weight[0] = 1.0;
	//	cdp_weight[1] = 1.0;
	//}
	//else if (numGP==3)
	//{
	//	cdp_point[0] = -sqrt(0.6);
	//	cdp_point[1] = 0;
	//	cdp_point[2] = sqrt(0.6);
	//	cdp_weight[0] = 5.0 / 9;
	//	cdp_weight[1] = 8.0 / 9;
	//	cdp_weight[2] = 5.0 / 9;
	//}
	//else if (numGP==4)
	//{
	//	cdp_point[0] = -sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5));
	//	cdp_point[1] = -sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5));
	//	cdp_point[2] = sqrt(3.0 / 7 - 2.0 / 7 * sqrt(6.0 / 5));
	//	cdp_point[3] = sqrt(3.0 / 7 + 2.0 / 7 * sqrt(6.0 / 5));
	//	cdp_weight[0] = 0.5 - sqrt(30.0) / 36.0;
	//	cdp_weight[1] = 0.5 + sqrt(30.0) / 36.0;
	//	cdp_weight[2] = 0.5 + sqrt(30.0) / 36.0;
	//	cdp_weight[3] = 0.5 - sqrt(30.0) / 36.0;
	//}
	//else if (numGP==5)
	//{
	//	cdp_point[0] = -1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10 / 7));
	//	cdp_point[1] = -1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10 / 7));
	//	cdp_point[2] = 0.0;
	//	cdp_point[3] = 1.0 / 3.0 * sqrt(5.0 - 2.0 * sqrt(10 / 7));
	//	cdp_point[4] = 1.0 / 3.0 * sqrt(5.0 + 2.0 * sqrt(10 / 7));
	//	cdp_weight[0] = (322.0 - 13 * sqrt(70.0)) / 900.0;
	//	cdp_weight[1] = (322.0 + 13 * sqrt(70.0)) / 900.0;
	//	cdp_weight[2] = 128.0 / 225.0;
	//	cdp_weight[3] = (322.0 + 13 * sqrt(70.0)) / 900.0;
	//	cdp_weight[4] = (322.0 - 13 * sqrt(70.0)) / 900.0;
	//}
	
}

pdGaussPt::~pdGaussPt()
{
	delete[] cdp_point, cdp_weight;
	cdp_point = NULL, cdp_weight = NULL;
}

   int pdGaussPt::i_getNumPts() const
{
	return ci_numGp;
}

   double pdGaussPt::d_getGaussPt(int gp) const
{
	return cdp_point[gp];
}

   double pdGaussPt::d_getWeight(int gp) const
{
	return cdp_weight[gp];
}
