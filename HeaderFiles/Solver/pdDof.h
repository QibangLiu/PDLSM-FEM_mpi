#pragma once
using namespace std;
class pdDof
{
public:
	pdDof();
	~pdDof();
	void setNotActive();
	bool b_isActive()const;
	void setValue(double u);
	double d_getValue() const;

	
	void setEqInde(int eqInde);
	int i_getEqInde()const;
private:
	bool cb_active;
	double cd_u;
	int ci_eqInde;// the index of equation, 0,1,2,3...if this dof is not active, ci_eqInde=-1;
	//int ci_fixInde;// index of this dof which is fixed;
	
	
};
