/* femGauss3Pt
This class stores the gauss point data for three point integration.

*/

#pragma once

//#include "stdafx.h"
#include<iomanip>
#include<cmath>
class pdGaussPt
{
public:
	//Constructor
	pdGaussPt(int numGP);
	//Destructor
	~pdGaussPt();

	// Functions
	int i_getNumPts() const;				//	return ci_numGp;
	double d_getGaussPt(int gp) const;	//	return cd_point[gp];
	double d_getWeight(int gp) const;		//	return cd_weight[gp];
private:
	//Constructor
	pdGaussPt();
	int ci_numGp;		//number of Gauss integration points, ci_numGp = 3;
	double *cdp_point;	//location of Gauss point, i.e. cd_point[0] = -0.774596669241483 etc.
	double *cdp_weight;	//weight, i e. cd_weight[0] = 0.555555555555556 etc.
};


