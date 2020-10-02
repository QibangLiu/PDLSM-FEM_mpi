#pragma once
/*This class is named data level 2, for storing the data be in contigous adress (1D array),
/so that they can be easy pass throungh by MPIand stroing;*/
using namespace std;
class dataLev2
{
public:
	dataLev2();
	~dataLev2();
	//====node data;
	double* cd_X;
	//===elemement data;
	int* ci_eleAlgoType;// algotype;
	int* ci_eleNodeID;//connective node id
	int* ci_ANN;//accumulated number of node; size =numELE+1;
private:

};

