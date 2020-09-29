#pragma once
using namespace std;

class pdPointBC
{
private:
	int ci_NodeID;
	int ci_fDof; //0 is FX, 1 is FY
	double cd_value;			//value of the point force
	pdPointBC();	//never used constructors
public:
	pdPointBC(int NID, int fdof, double value);
	int  i_getfDOF() const;		// return ce_pbctype;
	double d_getValue() const;			// return cd_value;
	int i_getNID()const;
	void setForceV(double f);
};

