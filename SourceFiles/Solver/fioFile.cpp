#include "fioFile.h"
#include<iostream>
using namespace std;


fioFiles::fioFiles(int rank)
{
	ci_rank = rank;
	ci_wflag = 3;
}

void fioFiles::CMDfile(datModel& o_dat, ifstream& fin)
{
	string sentence;
	vector<string> tokens;
	while (getline(fin, sentence))
	{
		istringstream iss(sentence);
		copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter(tokens));
		excuteCMD(o_dat, tokens);
		tokens.clear();
	}


}

void fioFiles::excuteCMD(datModel& o_dat, vector<string>& tokens)
{
	
	if (tokens.size() != 0)
	{
		if (tokens[0][0] != '#')// not comment cmd;
		{
			if (tokens[0]=="MSHFILE")  //input mesh file
			{
				ifstream fin;		//fin for input mesh file
				fin.open(tokens[1]);
				if (ci_rank == 0)
				{
					if (!fin.is_open())
					{
						printf("File %s is not exist\n", tokens[1]);
						exit(0);
					}
				}
				printf("Reading data on %d-th core....\n", ci_rank);
				o_dat.readdata(fin);
			}
			else if (tokens[0] == "SOLVER") //solver cmd
			{
				if (tokens.size()<2)
				{
					if (ci_rank==0)
					{
						printf("Miss solver. \n");
						exit(0);
					}
				}
				else
				{
					if (tokens[1]=="DYNAMIC")
					{
						o_dat.ci_solvFlag = 0;
					}
					else if (tokens[1] == "STATIC")
					{
						o_dat.ci_solvFlag = 1;
					}
					else if (tokens[1] == "QUASI-STATIC")
					{
						o_dat.ci_solvFlag = 2;
					}
					else
					{
						if (ci_rank == 0)
						{
							printf("Don't exist command of %s.\n", tokens[1]);
							exit(0);
						}
					}
				}
			}
			else if (tokens[0] == "SETSOLVING") //setsolver cmd
			{
				if (tokens.size() < 5)
				{
					if (ci_rank == 0)
					{
						printf("Miss SETSOLVING CMD. \n");
						exit(0);
					}
				}
				else
				{
					o_dat.cd_dt = atof(tokens[1].c_str());
					o_dat.ci_numTstep= atoi(tokens[2].c_str());
					o_dat.ci_savefrequence= atoi(tokens[3].c_str());
					o_dat.op_getGeomP()->cd_factor= atof(tokens[4].c_str());
				}
			}
			else if (tokens[0] == "VEBC")// essential bc
			{
				if (tokens.size() < 3)
				{
					if (ci_rank == 0)
					{
						printf("Miss VEBC CMD. \n");
						exit(0);
					}
				}
				else
				{
					int ID;
					double velo;
					ID= atoi(tokens[1].c_str());
					velo= atof(tokens[2].c_str());
					o_dat.op_getEssenBC(ID)->cb_varing = true;
					o_dat.op_getEssenBC(ID)->cd_velocity = velo;
				}
			}
			else if (tokens[0] == "VNBC")// natural bc
			{
				if (tokens.size() < 3)
				{
					if (ci_rank == 0)
					{
						printf("Miss VNBC CMD. \n");
						exit(0);
					}
				}
				else
				{
					int ID;
					double velo;
					ID = atoi(tokens[1].c_str());
					velo = atof(tokens[2].c_str());
					o_dat.op_getNaturalBC(ID)->cb_varing = true;
					o_dat.op_getNaturalBC(ID)->cd_velocity = velo;
				}
			}
			else
			{
				if (ci_rank == 0)
				{
					printf("WARNING: ignore unknown command of %s. \n", tokens[0].c_str());
				}
			}
			
		}


	}
}

void fioFiles::writeResults(datModel &o_dat)
{
	struct stat info;
	char s_path[_MAX_PATH];
	getcwd(s_path, _MAX_PATH);
	sprintf(s_path, "%s/Results", s_path);
	int i_stat=stat(s_path, &info);
	if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows 
	{
		//printf("%s is a directory\n", s_path);
		if (i_stat != 0)
		{
			printf("%s is a directory but cannot access.\n", s_path);
		}
	}
	else
	{
		printf("%s is no exist, and made it automaticly.\n", s_path);
		#ifdef __linux__
			mkdir(s_path, 777); 
		#else
			cout<<"=="<<s_path<<endl;
			mkdir(s_path);
		#endif
	}
	if (ci_wflag==1)
	{
		writeReslutsTOTAL_vtk(o_dat,0);
	}
	else if (ci_wflag ==2)
	{
		writeResultsPD_vtk(o_dat);
		writeResultsFEM_vtk(o_dat);
	}
	else if (ci_wflag ==3)
	{
		writeReslutsTOTAL_vtk(o_dat,0);
		writeResultsPD_vtk(o_dat);
		writeResultsFEM_vtk(o_dat);
		writeUofNode(o_dat);
		writeSigofNode(o_dat);
	}
	

}


void fioFiles::writeUofNode(datModel &o_dat)
{
	char fileNAME[_MAX_PATH];
	getTitle(o_dat, fileNAME);
	sprintf(fileNAME, "%s.disp", fileNAME);
	ofstream fout(fileNAME);
	fout << "*****The final displacement of  point*****" << endl;
	fout << std::left;
	fout<< "ID" << "\t"
		<< "x" << "\t" << "y" << "\t" << "z" <<
		"\t" << "Ux" << "\t" << "Uy" << "\t" << "Uz" << endl;


	for (int k = 0; k < o_dat.getTotnumNode(); k++)
	{

		o_dat.op_getNode(k)->printFinalU(fout);

	}
	fout.close();
}

void fioFiles::getTitle(datModel &o_dat, char Titl[])
{
	char s_path[_MAX_PATH];
	getcwd(s_path, _MAX_PATH);
	sprintf(Titl, "%s/Results/%s", s_path, o_dat.cs_title.c_str());

	
}

void fioFiles::writeSigofNode(datModel &o_dat)
{
	char fileNAME[_MAX_PATH];
	getTitle(o_dat, fileNAME);
	sprintf(fileNAME,"%s.sig", fileNAME);
	ofstream fout(fileNAME);
	fout << "*****The stress of  point*****" << endl;
	fout << std::left;
	fout << "NID" << "\t"
		<< "x" << "\t" << "y" << "\t" << "z" <<
		"\t" << "Sxx" << "\t" << "Syy" <<
		"\t" << "Szz" << "\t" << "Sxy" << "\t" << "Syz" <<
		"\t" << "Szx" << endl;
	fout << setiosflags(ios::scientific) << setprecision(5);
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		o_dat.op_getNode(i)->printStress(fout);
	}
	fout.close();
}

void fioFiles::writeReslutsTOTAL_vtk(datModel &o_dat,int nTstep)
{
	char fileNAME[_MAX_PATH];
	getTitle(o_dat, fileNAME);
	sprintf(fileNAME, "%s_%d.vtk", fileNAME, nTstep);
	ofstream fout(fileNAME);
	//=============header, title, data type( ASCII or BINARY)============
	fout << setiosflags(ios::scientific) << setprecision(4);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << o_dat.cs_title << endl;
	fout << "ASCII" << endl;
	//=================Geometry/topology=================
	// ==point data==
	double xN[3];
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << o_dat.getTotnumNode() << " float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		o_dat.op_getNode(i)->getcoor(xN);
		fout << xN[0] << ' ' << xN[1] << ' ' << xN[2] << endl;
	}
	fout << endl;
	//===element data==;
	// connectivity;
	int sizeList = 0;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		sizeList = sizeList + o_dat.op_getEles(i)->getNumNodes_vtk();
	}
	sizeList = sizeList + o_dat.getTotnumEle();
	fout << "CELLS " << o_dat.getTotnumEle() << ' ' << sizeList << endl;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		o_dat.op_getEles(i)->print_vtk(fout, o_dat.op_getEles(i)->cip_EleNodeId);
	}
	fout << endl;
	//element type;
	fout << "CELL_TYPES " << o_dat.getTotnumEle() << endl;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		fout<<o_dat.op_getEles(i)->ci_eleType << endl;
	}
	//================Dataset attributes============
	//displacement;
	fout << "POINT_DATA " << o_dat.getTotnumNode() << endl;
	fout << "VECTORS " << "U " << "float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fout << o_dat.op_getNode(i)->op_getDof(j)->d_getValue() << ' ';
		}
		fout << endl;
	}
	//stress;
	fout << "TENSORS " << "sigma " << "float\n";
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		o_dat.op_getNode(i)->printStressTensor_vtk(fout);
		fout << endl;
	}
	//phi
	fout << "SCALARS " << "phi " << "float " << 1 << endl;
	fout << "LOOKUP_TABLE  " << "default  "<< endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		fout << o_dat.op_getNode(i)->getLocalDamage() << ' ';
		if (i%20==0)
		{
			fout << endl;
		}
	}
	fout.close();
}

void fioFiles::writeResultsPD_vtk(datModel& o_dat)
{
	
	int* mapNodeID = new int[o_dat.getTotnumNode()];
	char fileNAME[_MAX_PATH];
	getTitle(o_dat, fileNAME);
	sprintf(fileNAME, "%s_PD.vtk", fileNAME);
	ofstream fout(fileNAME);
	//=============header, title, data type( ASCII or BINARY)============
	fout << setiosflags(ios::scientific) << setprecision(4);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << o_dat.cs_title << endl;
	fout << "ASCII" << endl;
	//=================Geometry/topology=================
	// ==PD point data==
	int count = 1;//caution, one- based;
	double xN[3];
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << o_dat.getTotnumFami() << " float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType())
		{
			o_dat.op_getNode(i)->getcoor(xN);
			fout << xN[0] << ' ' << xN[1] << ' ' << xN[2] << endl;
			mapNodeID[i] = count;
			count++;
		}
		
	}
	fout << endl;
	//===element data==;
	// connectivity;
	int sizeList = 0, algoType;
	count = 0;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		algoType = o_dat.op_getEles(i)->ci_algorithmType;
		if (algoType==1)
		{
			sizeList = sizeList + o_dat.op_getEles(i)->getNumNodes_vtk();
			count++;
		}
		
	}
	sizeList = sizeList + count;
	fout << "CELLS " << count << ' ' << sizeList << endl;
	int* mapEleNodeID, numNodeEle;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		algoType = o_dat.op_getEles(i)->ci_algorithmType;
		if (algoType==1)
		{
			numNodeEle = o_dat.op_getEles(i)->ci_numNodes;
			mapEleNodeID = new int[numNodeEle];
			for (int ele = 0; ele < numNodeEle; ele++)
			{
				mapEleNodeID[ele] = mapNodeID[(o_dat.op_getEles(i)->cip_EleNodeId)[ele] - 1];
			}
			o_dat.op_getEles(i)->print_vtk(fout, mapEleNodeID);
			delete[] mapEleNodeID;
			mapEleNodeID = NULL;
		}
	}
	fout << endl;
	//element type;
	fout << "CELL_TYPES " << count << endl;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		algoType = o_dat.op_getEles(i)->ci_algorithmType;
		if (algoType == 1)
		{
			fout << o_dat.op_getEles(i)->ci_eleType << endl;
		}
	}
	//================Dataset attributes============
	//displacement;
	fout << "POINT_DATA " << o_dat.getTotnumFami() << endl;
	fout << "VECTORS " << "U " << "float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType())// pd node;
		{
			for (int j = 0; j < 3; j++)
			{
				fout << o_dat.op_getNode(i)->op_getDof(j)->d_getValue() << ' ';
			}
			fout << endl;
		}
		
	}
	//stress;
	fout << "TENSORS " << "sigma " << "float\n";
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType())// pd node;
		{
			o_dat.op_getNode(i)->printStressTensor_vtk(fout);
			fout << endl;
		}
	}
	fout.close();
	delete[] mapNodeID;
	mapNodeID = NULL;
}

void fioFiles::writeResultsFEM_vtk(datModel& o_dat)
{
	//=======count number of fem node;
	int numFEMnode = 0;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType() != 2)
		{
			numFEMnode++;
		}
	}
	int* mapNodeID = new int[o_dat.getTotnumNode()];
	//===file name;
	char fileNAME[_MAX_PATH];
	getTitle(o_dat, fileNAME);
	sprintf(fileNAME, "%s_FEM.vtk", fileNAME);
	ofstream fout(fileNAME);
	//=============header, title, data type( ASCII or BINARY)============
	fout << setiosflags(ios::scientific) << setprecision(4);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << o_dat.cs_title << endl;
	fout << "ASCII" << endl;
	//=================Geometry/topology=================
	// ==PD point data==
	int count = 1;//caution, one- based;
	double xN[3];
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << numFEMnode << " float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType() != 2)
		{
			o_dat.op_getNode(i)->getcoor(xN);
			fout << xN[0] << ' ' << xN[1] << ' ' << xN[2] << endl;
			mapNodeID[i] = count;
			count++;
		}
	}
	fout << endl;
	//===element data==;
	// connectivity;
	int sizeList = 0, algoType;
	count = 0;//count for number of fem element;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		algoType = o_dat.op_getEles(i)->ci_algorithmType;
		if (algoType == 2)
		{
			sizeList = sizeList + o_dat.op_getEles(i)->getNumNodes_vtk();
			count++;
		}
	}
	sizeList = sizeList + count;
	fout << "CELLS " << count << ' ' << sizeList << endl;
	int* mapEleNodeID, numNodeEle;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		algoType = o_dat.op_getEles(i)->ci_algorithmType;
		if (algoType == 2)
		{
			numNodeEle = o_dat.op_getEles(i)->ci_numNodes;
			mapEleNodeID = new int[numNodeEle];
			for (int ele = 0; ele < numNodeEle; ele++)
			{
				mapEleNodeID[ele] = mapNodeID[(o_dat.op_getEles(i)->cip_EleNodeId)[ele] - 1];
			}
			o_dat.op_getEles(i)->print_vtk(fout, mapEleNodeID);
			delete[] mapEleNodeID;
			mapEleNodeID = NULL;
		}
	}
	fout << endl;
	//element type;
	fout << "CELL_TYPES " << count << endl;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		algoType = o_dat.op_getEles(i)->ci_algorithmType;
		if (algoType == 2)
		{
			fout << o_dat.op_getEles(i)->ci_eleType << endl;
		}
	}
	//================Dataset attributes============
	//displacement;
	fout << "POINT_DATA " << numFEMnode << endl;
	fout << "VECTORS " << "U " << "float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType()!=2)// pd node;
		{
			for (int j = 0; j < 3; j++)
			{
				fout << o_dat.op_getNode(i)->op_getDof(j)->d_getValue() << ' ';
			}
			fout << endl;
		}
	}
	//stress;
	fout << "TENSORS " << "sigma " << "float\n";
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType()!=2)// pd node;
		{
			o_dat.op_getNode(i)->printStressTensor_vtk(fout);
			fout << endl;
		}
	}
	fout.close();
	delete[]mapNodeID;
	mapNodeID = NULL;
}