#include "pdVaryEssentialBCs.h"



pdVaryEssentialBCs::pdVaryEssentialBCs(int Nid, int direc)
{
	ci_Nid = Nid;
	ci_direction = direc;
}


pdVaryEssentialBCs::~pdVaryEssentialBCs()
{
}

int pdVaryEssentialBCs::getNid() const
{
	return ci_Nid;
}

int pdVaryEssentialBCs::getDirec() const
{
	return ci_direction;
}
