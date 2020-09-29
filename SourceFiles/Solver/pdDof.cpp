#include "pdDof.h"

pdDof::pdDof()
{
	cb_active = true;
	cd_u = 0;
	ci_eqInde = 0;
	
}

pdDof::~pdDof()
{
}

void pdDof::setNotActive()
{
	cb_active = false;
}

bool pdDof::b_isActive() const
{
	return cb_active;
}

void pdDof::setValue(double u)
{
	cd_u = u;
}

double pdDof::d_getValue() const
{
	return cd_u;
}

void pdDof::setEqInde(int eqInde)
{
	ci_eqInde = eqInde;
}

int pdDof::i_getEqInde() const
{
	return ci_eqInde;
}


