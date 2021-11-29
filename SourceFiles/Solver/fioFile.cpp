/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */


#include "fioFile.h"
#include<iostream>
using namespace std;


fioFiles::fioFiles(int rank)
{
	ci_rank = rank;
	ci_wflag = 3;
	cb_littleEndian = is_little_endian();
	// make dir ~/Results to store results
	struct stat info;
	char s_path[_MAX_PATH];
	getcwd(s_path, _MAX_PATH);
	sprintf(s_path, "%s/Results", s_path);
	int i_stat = stat(s_path, &info);
	if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows 
	{
		//printf("%s is a directory\n", s_path);
		if (i_stat != 0 && ci_rank == 0)
		{
			printf("%s is a directory but cannot access.\n", s_path);
		}
	}
	else
	{
		if (ci_rank == 0)
		{
			printf("%s is no exist, and made it automaticly.\n", s_path);
		}
#ifdef __linux__
		mkdir(s_path, 777);
#elif __APPLE__
		mkdir(s_path, 777);
#else
		mkdir(s_path);
#endif
	}
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
				//if (ci_rank == 0)
				{
					if (!fin.is_open())
					{
						printf("ERROR: File \"%s\" is not exist.\n", tokens[1].c_str());
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
						printf("ERROR: Miss solver. \n");
					}
					exit(0);
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
							printf("ERROR: Don't exist command of \"%s\".\n", tokens[1].c_str());
						}
						exit(0);
					}
				}
			}
			else if (tokens[0] == "SETSOLVING") //setsolver cmd
			{
				if (tokens.size() < 6)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss SETSOLVING CMD. \n");
					}
				}
				else
				{
					o_dat.cd_dt = atof(tokens[1].c_str());
					if (o_dat.ci_solvFlag==2)
					{
						o_dat.cd_dt = 1;
					}
					o_dat.ci_numTstep= atoi(tokens[2].c_str());
					o_dat.ci_savefrequence= atoi(tokens[3].c_str());
					o_dat.op_getGeomP()->cd_factor= atof(tokens[4].c_str());
					o_dat.cd_NLF = atof(tokens[5].c_str());
					if (tokens.size()>6)
					{
						o_dat.cd_me = atof(tokens[6].c_str());
					}
				}
			}
			else if (tokens[0] == "NEWMARK") //setsolver cmd
			{
				if (tokens.size() < 2)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss NEWMARK CMD. Command is ignord. \n");
						//exit(0);
					}
				}
				else if (tokens[1] == "ON")
				{
					o_dat.cb_Newmark = true;
				}
				else if (tokens[1] == "OFF")
				{
					o_dat.cb_Newmark = false;
				}
				else
				{
					if (ci_rank == 0)
					{
						printf("WARNING: ignore unknown command: \"%s %s\"! \n", tokens[0].c_str(), tokens[1].c_str());
					}
				}
				if (tokens.size()>3)
				{
					o_dat.cd_gamma= atof(tokens[2].c_str());
					o_dat.cd_beta= atof(tokens[3].c_str());
				}
			}
			else if (tokens[0] == "FENSF") //setsolver cmd
			{
				if (tokens.size() < 2)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss FENSF CMD. Command is ignord. \n");
						//exit(0);
					}
				}
				else if (tokens[1] == "ON")
				{
					o_dat.cb_FENSF = true;
				}
				else if (tokens[1] == "OFF")
				{
					o_dat.cb_FENSF = false;
				}
				else
				{
					if (ci_rank == 0)
					{
						printf("WARNING: ignore unknown command: \"%s %s\"! \n", tokens[0].c_str(), tokens[1].c_str());
					}
				}
			}
			else if (tokens[0] == "LUMPMASS") //setsolver cmd
			{
				if (tokens.size() < 2)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss NEWMARK CMD. Command is ignord. \n");
					}
				}
				else if (tokens[1] == "ON")
				{
					o_dat.cb_lumpedMass = true;
				}
				else if (tokens[1] == "OFF")
				{
					o_dat.cb_lumpedMass = false;
				}
				else
				{
					if (ci_rank == 0)
					{
						printf("WARNING: ignore unknown command: \"%s %s\"! \n", tokens[0].c_str(), tokens[1].c_str());
					}
				}
			}
			else if (tokens[0] == "VEBC")// essential bc
			{
				if (tokens.size() < 3)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss VEBC CMD. Command is ignord. \n");
					}
				}
				else
				{
					int ID;
					double velo;
					ID= atoi(tokens[1].c_str());
					if (ID>-1&&ID<o_dat.getTotnumEssentialBCs())
					{
						velo = atof(tokens[2].c_str());
						o_dat.op_getEssenBC(ID)->cb_varing = true;
						o_dat.op_getEssenBC(ID)->cd_velocity = velo;
						o_dat.op_getEssenBC(ID)->cd_maxV = velo;
					}
					else if (ci_rank==0)
					{
						printf("Warning: the setting of essencital BC %d is ignored!\n", ID);
					}
				
				}
			}
			else if (tokens[0] == "EBC")// essential bc
			{
				if (tokens.size() < 3)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss EBC CMD. Command is ignord. \n");
					}
				}
				else
				{
					int ID;
					double val;
					ID = atoi(tokens[1].c_str());
					if (ID > -1 && ID < o_dat.getTotnumEssentialBCs())
					{
						val = atof(tokens[2].c_str());
						o_dat.setEssentialBC(ID, val);
					}
					else if (ci_rank == 0)
					{
						printf("Warning: the setting of essencital BC %d is ignored!\n", ID);
					}

				}
			}
			else if (tokens[0] == "VNBC")// natural bc
			{
				if (tokens.size() < 3)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss VNBC CMD. Command is ignord. \n");
					}
				}
				else
				{
					int ID;
					double velo;
					ID = atoi(tokens[1].c_str());
					if (ID>-1&&ID<o_dat.getTotnumNaturalBCs())
					{
						velo = atof(tokens[2].c_str());
						o_dat.op_getNaturalBC(ID)->cb_varing = true;
						o_dat.op_getNaturalBC(ID)->cd_velocity = velo;
					}
					else if (ci_rank == 0)
					{
						printf("Warning: the setting of natural BC %d is ignored!\n", ID);
					}
					
				}
			}
			else if (tokens[0] == "NBC")// natural bc
			{
			if (tokens.size() < 3)
			{
				if (ci_rank == 0)
				{
					printf("Warning: Miss NBC CMD, command is ignord.\n");
					//exit(0);
				}
			}
			else
			{
				int ID;
				double val;
				ID = atoi(tokens[1].c_str());
				if (ID > -1 && ID < o_dat.getTotnumNaturalBCs())
				{
					val = atof(tokens[2].c_str());
					o_dat.op_getNaturalBC(ID)->cd_value = val;
				}
				else if (ci_rank == 0)
				{
					printf("Warning: the setting of natural BC %d is ignored!\n", ID);
				}

			}
			}
			else if (tokens[0] == "FC")
			{
				if (tokens.size() < 2)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss failure flag, command is ignord. \n");
						//exit(0);
					}
				}
				else
				{
					int failFlag= atoi(tokens[1].c_str());
					o_dat.ci_failFlag = failFlag;
					if (failFlag < 0 || failFlag>3)
					{
						if (ci_rank == 0)
						{
							printf("Warning: the failure criterion %d is ignored!\n", failFlag);
						}
					}
					else if (failFlag==1)
					{
						//Flags: 1--maximum circumferential tensile stress;
						// may set m_r, d_c;
						if (ci_rank == 0 && o_dat.ci_Numdimen != 2)
						{
							printf("Warning: the maximum circumferential tensile stress failure criterion is for 2D problems!\n", failFlag);
						}
						else if (tokens.size()>3)
						{
							o_dat.cd_mr = atof(tokens[2].c_str());
							o_dat.cd_dcf= atof(tokens[3].c_str());
						}
						
					}
					else if (failFlag == 2)
					{
						//#Flags: 2--maximum principal stress.
						if (tokens.size() > 2)
						{
							o_dat.ci_topk = atoi(tokens[2].c_str());
						}
						
					}
				}
			}
			else if (tokens[0] == "TOPK")
			{
				if (tokens.size() < 2)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss k value for TOPK. \n");
					}
				}
				else
				{
					o_dat.ci_topk = atoi(tokens[1].c_str());
				}
			}
			else if (tokens[0] == "RF")
			{
				if (tokens.size() < 2)
				{
					if (ci_rank == 0)
					{
						printf("Warning: Miss ID value for RF. \n");
					}
				}
				else
				{
					int id;
					for (int i = 1; i < tokens.size(); i++)
					{
						id = atoi(tokens[i].c_str());
						o_dat.civ_reacForceOfessBCId.push_back(id);
					}
					o_dat.setReacForcNode();
				}
			}
			else if (tokens[0] == "VTKFORMAT")
			{
			if (tokens.size() < 2)
			{
				if (ci_rank == 0)
				{
					printf("Warning: Miss formate value for VTKFORMAT. \n");
				}
			}
			else
			{
				if (tokens[1] == "ASCII")
				{
					o_dat.cb_vtkBinary = false;
				}
				else if (tokens[1] == "BINARY")
				{
					o_dat.cb_vtkBinary = true;;
				}
				else if (ci_rank == 0)
				{
					printf("WARNING: ignore unknown command \"%s %s\". \n", tokens[0].c_str(), tokens[1].c_str());
				}
			
			}
			}
			else
			{
				if (ci_rank == 0)
				{
					printf("WARNING: ignore unknown command \"%s\". \n", tokens[0].c_str());
				}
			}
			
		}


	}
	//======set the filename here;
	setTitle(o_dat);
}

void fioFiles::writeResults(datModel &o_dat)
{
	if (ci_wflag==1)
	{
		writeReslutsTOTAL_vtk(o_dat,"static");
	}
	else if (ci_wflag ==2)
	{
		writeResultsPD_vtk(o_dat);
		writeResultsFEM_vtk(o_dat);
	}
	else if (ci_wflag ==3)
	{
		writeReslutsTOTAL_vtk(o_dat, "static");
		writeResultsPD_vtk(o_dat);
		writeResultsFEM_vtk(o_dat);
	}
	

}


void fioFiles::writeUofNode(datModel &o_dat, string nTstep)
{
	char fileNAME[_MAX_PATH];
	sprintf(fileNAME, "%s_%s.disp", cc_fileNAME, nTstep.c_str());
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

void fioFiles::setTitle(datModel &o_dat)
{
	char s_path[_MAX_PATH];
	getcwd(s_path, _MAX_PATH);
	sprintf(cc_fileNAME, "%s/Results/%s", s_path, o_dat.cs_title.c_str());
	o_dat.cs_fileName = cc_fileNAME;

	
}

bool fioFiles::is_little_endian()
{
	union w
	{
		int a;
		char b;
	}c;
	c.a = 1;
	return (c.b == 1);
}

void fioFiles::writeSigofNode(datModel &o_dat, string nTstep)
{
	char fileNAME[_MAX_PATH];
	sprintf(fileNAME,"%s_%s.sig", cc_fileNAME, nTstep.c_str());
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

void fioFiles::writeReslutsTOTAL_vtk(datModel &o_dat, string nTstep)
{
	if (o_dat.cb_vtkBinary)
	{
		writeReslutsTOTAL_vtk_Binary(o_dat, nTstep);
	}
	else
	{
		writeReslutsTOTAL_vtk_ASCII(o_dat, nTstep);
		writeSigofNode(o_dat, nTstep);
		writeUofNode(o_dat, nTstep);
	}
}

void fioFiles::writeReslutsTOTAL_vtk_ASCII(datModel& o_dat, string nTstep)
{
	char fileNAME[_MAX_PATH];
	sprintf(fileNAME, "%s_%s.vtk", cc_fileNAME, nTstep.c_str());
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
		fout << o_dat.op_getEles(i)->ci_eleType << endl;
	}
	//================Dataset attributes============
	//displacement;
	fout << "POINT_DATA " << o_dat.getTotnumNode() << endl;
	fout << "VECTORS " << "U " << "float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fout << float(o_dat.op_getNode(i)->op_getDof(j)->d_getValue()) << ' ';
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
	fout << "LOOKUP_TABLE  " << "default  " << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		fout << float(o_dat.op_getNode(i)->getLocalDamage()) << ' ';
		if (i % 20 == 0)
		{
			fout << endl;
		}
	}
	fout.close();
}

void fioFiles::writeReslutsTOTAL_vtk_Binary(datModel& o_dat, string nTstep)
{
	/*write results in vtk file with Binary format;
	vtk using big endian to write binary file;
	if machine is little endial, swap byte is needed;*/
	char fileNAME[_MAX_PATH];
	sprintf(fileNAME, "%s_%s.vtk", cc_fileNAME, nTstep.c_str());
	ofstream fout(fileNAME, std::ios::binary);
	//=============header, title, data type( ASCII or BINARY)============
	//fout << setiosflags(ios::scientific) << setprecision(4);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << o_dat.cs_title << endl;
	fout << "BINARY" << endl;
	//=================Geometry/topology=================
	// ==point data==
	double xN[3];
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << o_dat.getTotnumNode() << " float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		o_dat.op_getNode(i)->getcoor(xN);
		float X[3] = { float(xN[0]), float(xN[1]),float(xN[2]) };
		if (cb_littleEndian)
		{
			SwapEnd(X[0]);
			SwapEnd(X[1]);
			SwapEnd(X[2]);
		}
		fout.write(reinterpret_cast<char*>(X), sizeof(float)*3);
		/*fout.write(reinterpret_cast<char*>(&X[[0]), sizeof(float));
		fout.write(reinterpret_cast<char*>(&X[[1]), sizeof(float));
		fout.write(reinterpret_cast<char*>(&X[[2]), sizeof(float));*/
	}
	fout << endl;
	//===element data==;
	// connectivity;
	int sizeList = 0, Nen, * conN;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		sizeList = sizeList + o_dat.op_getEles(i)->getNumNodes_vtk();
	}
	sizeList = sizeList + o_dat.getTotnumEle();
	fout << "CELLS " << o_dat.getTotnumEle() << ' ' << sizeList << endl;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		Nen = o_dat.op_getEles(i)->getNumNodes();
		int Nen1 = Nen;
		conN = new int[Nen];
		o_dat.op_getEles(i)->getConNid(conN);
		for (int j = 0; j < Nen; j++)
		{
			conN[j]--;
		}
		if (cb_littleEndian)
		{
			SwapEnd(Nen1);
			for (int j = 0; j < Nen; j++)
			{
				SwapEnd(conN[j]);
			}
		}
		fout.write(reinterpret_cast<char*>(&Nen1), sizeof(int));
		fout.write(reinterpret_cast<char*>(conN), sizeof(int)* Nen);
		delete[] conN;
		conN = NULL;
	}
	fout << endl;
	//element type;
	fout << "CELL_TYPES " << o_dat.getTotnumEle() << endl;
	for (int i = 0; i < o_dat.getTotnumEle(); i++)
	{
		int eleT = o_dat.op_getEles(i)->ci_eleType;
		if (cb_littleEndian)
		{
			SwapEnd(eleT);
		}
		fout.write(reinterpret_cast<char*>(&eleT), sizeof(int));
	}
	fout << endl;
	////================Dataset attributes============
	////displacement;
	fout << "POINT_DATA " << o_dat.getTotnumNode() << endl;
	fout << "VECTORS " << "U " << "float" << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		float U[3] = {float(o_dat.op_getNode(i)->op_getDof(0)->d_getValue()), 
		float(o_dat.op_getNode(i)->op_getDof(1)->d_getValue()),
		float(o_dat.op_getNode(i)->op_getDof(2)->d_getValue()) };
		if (cb_littleEndian)
		{
			SwapEnd(U[0]);
			SwapEnd(U[1]);
			SwapEnd(U[2]);
		}
		fout.write(reinterpret_cast<char*>(U), sizeof(float) * 3);
	}
	fout << endl;
	////stress;
	fout << "TENSORS " << "sigma " << "float\n";
	double sig[6];
	int idx[9] = { 0,3,5,3,1,4,5,4,2 };
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		o_dat.op_getNode(i)->getStress(sig);
		float sigma[6] = { float(sig[0]),float(sig[1]),float(sig[2]),float(sig[3]),float(sig[4]),float(sig[5]) };
		if (cb_littleEndian)
		{
			for (int j = 0; j < 6; j++)
			{
				SwapEnd(sigma[j]);
			}
		}
		for (int j = 0; j < 9; j++)
		{
			fout.write(reinterpret_cast<char*>(&sigma[idx[j]]), sizeof(float));
		}
		
	}
	fout << endl;
	//phi
	fout << "SCALARS " << "phi " << "float " << 1 << endl;
	fout << "LOOKUP_TABLE  " << "default  " << endl;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		float phi = o_dat.op_getNode(i)->getLocalDamage();
		if (cb_littleEndian)
		{
			SwapEnd(phi);
		}
		fout.write(reinterpret_cast<char*>(&phi), sizeof(float));
	}
	fout << endl;
	fout.close();
}

void fioFiles::writeResultsPD_vtk(datModel& o_dat)
{
	
	int* mapNodeID = new int[o_dat.getTotnumNode()];
	char fileNAME[_MAX_PATH];
	sprintf(fileNAME, "%s_PD.vtk", cc_fileNAME);
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
				fout << float(o_dat.op_getNode(i)->op_getDof(j)->d_getValue()) << ' ';
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
	//phi
	fout << "SCALARS " << "phi " << "float " << 1 << endl;
	fout << "LOOKUP_TABLE  " << "default  " << endl;
	count = 0;
	for (int i = 0; i < o_dat.getTotnumNode(); i++)
	{
		if (o_dat.op_getNode(i)->getNodeType())// pd node;
		{
			fout << float(o_dat.op_getNode(i)->getLocalDamage()) << ' ';
			if (count % 20 == 0)
			{
				fout << endl;
			}
			count++;
		}
	}
	fout.close();
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
	sprintf(fileNAME, "%s_FEM.vtk", cc_fileNAME);
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
				fout << float(o_dat.op_getNode(i)->op_getDof(j)->d_getValue()) << ' ';
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