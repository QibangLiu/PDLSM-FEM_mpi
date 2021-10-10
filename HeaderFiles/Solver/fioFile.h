/* fioFiles.h
fils operations.  files for input and output

*/

#pragma once

//#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<string>
#include <sstream>
#include <iterator>
#include"datModel.h"
#ifdef _WIN32
#include <direct.h>
#elif __linux
#include <unistd.h>
#include <sys/stat.h>
#define _MAX_PATH	777
#elif __APPLE__
#include <unistd.h>
#include <sys/stat.h>
#define _MAX_PATH	777
#endif


using namespace std;

class fioFiles
{
public:
	fioFiles(int rank);
	void CMDfile(datModel& o_dat,ifstream &fin);
	void excuteCMD(datModel& o_dat, vector<string>& tokens);
	void writeResults(datModel &o_dat);
	void writeReslutsTOTAL_vtk(datModel& o_dat,  string nTstep);
	void writeReslutsTOTAL_vtk_ASCII(datModel& o_dat, string nTstep);
	void writeReslutsTOTAL_vtk_Binary(datModel& o_dat, string nTstep);
private:
	void writeUofNode(datModel &o_dat);
	void writeSigofNode(datModel &o_dat);
	void writeResultsPD_vtk(datModel &o_dat);
	void writeResultsFEM_vtk(datModel& o_dat);
	void setTitle(datModel &o_dat);
	char cc_fileNAME[_MAX_PATH];
	bool is_little_endian();
	bool cb_littleEndian;
	int ci_rank;
	int ci_wflag;//write flag;
	template <typename T>
	void SwapEnd(T& var);

};

template<typename T>
inline void fioFiles::SwapEnd(T& var)
{
	char* varArray = reinterpret_cast<char*>(&var);
	for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
		std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}
