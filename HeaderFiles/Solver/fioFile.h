/* fioFiles.h
fils operations.  files for input and output

*/

#pragma once

//#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<string>
#include"datModel.h"
#ifdef _WIN32
#include <direct.h>
#elif __linux
#include <unistd.h>
#include <sys/stat.h>
#define _MAX_PATH	777
#endif


using namespace std;

class fioFiles
{
public:
	void writeResults(datModel &o_dat, int wflag);
private:
	void writeUofNode(datModel &o_dat);
	void writeSigofNode(datModel &o_dat);
	void writeReslutsTOTAL_vtk(datModel &o_dat);
	void writeResultsPD_vtk(datModel &o_dat);
	void writeResultsFEM_vtk(datModel& o_dat);
	void getTitle(datModel &o_dat,char Titl[]);

};
