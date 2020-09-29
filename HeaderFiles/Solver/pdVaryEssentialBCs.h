#pragma once
class pdVaryEssentialBCs
{
public:
	pdVaryEssentialBCs(int Nid,int direc);
	~pdVaryEssentialBCs();
	int getNid()const;
	int getDirec()const;
private:
	int ci_Nid;
	int ci_direction;
};

