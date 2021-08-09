/* ----------------------------------------------------------------------------
   PDLSM-FEM: Peridynamics least squares minimization-finite element method

   Author: Qibang Liu, qibangliu@ksu.edu
   Department of mechanical and nuclear engineering, Kansas State University

   Copyright (2021). This software is distributed under
   the GNU General Public License.

   See the PDLSM-FEM_manual.pdf file in the top-level directory for instructions.
-------------------------------------------------------------------------------- */

#include "pdfemEles.h"

pdfemEles::pdfemEles(int id, int numNodes, int *nId,int algoType, dataLev2* p_datLev2)
{
	cop_datLev2 = p_datLev2;
	ci_EleId = id;
	ci_numNodes = numNodes;
	cip_EleNodeId = new int[numNodes];
	for (int i = 0; i < numNodes; i++)
	{
		cip_EleNodeId[i] = nId[i];
	}
	ci_algorithmType = algoType;
}

pdfemEles::~pdfemEles()
{
	delete[] cip_EleNodeId;
	cip_EleNodeId = NULL;
}

void pdfemEles::getConNid(int Nid[])
{
	for (int i = 0; i < ci_numNodes; i++)
	{
		Nid[i] = cip_EleNodeId[i];
	}
}


void pdfemEles::print(ofstream & fout)
{
	fout << setw(6)<<ci_EleId;
	for (int i = 0; i < ci_numNodes; i++)
	{
		fout << setw(6) << cip_EleNodeId[i];
	}
	fout << endl;
}

//void pdfemEles::print_vtk(ofstream& fout)
//{
//	//vtk element is zero-based;
//	switch (ci_eleType)
//	{
//	default:
//		cout << "ERROR: element type is wrong\n";
//		exit(0);
//	case 5:
//	{
//		fout << 3 << ' ';
//		for (int i = 0; i < 3; i++)
//		{
//			fout << cip_EleNodeId[i] - 1 << ' ';
//		}
//		fout << endl;
//		break; // 3;   triangle
//	}
//	case 22:
//
//		break;// 6;  6n_triangle
//	case 9:
//	{
//		fout << 4 << ' ';
//		for (int i = 0; i < 4; i++)
//		{
//			fout << cip_EleNodeId[i] - 1 << ' ';
//		}
//		fout << endl;
//		break;// 4;  quad
//	}
//
//	case 23:
//		break;// 8;  8n_quad;
//	case 10:
//		break;// 4;  4n_tetra
//	case 24:
//		break;// 10; 4n_tetra
//	case 12:
//		break;// 8; hexahedron
//	case 25:
//		break;// 20;  20n_HEXAHEDRON
//	case 13:
//		break;// 6;  wedge;
//	case 26:
//		break;// 15;  15n_wedge
//	case 14:
//		break;// 5; pyramid;
//	case 27:
//		break;// 13; 13n_PYRAMID;
//	}
//
//}




int pdfemEles::getNumNodes_vtk() const
{
	//switch (ci_eleType)
	//{
	//default:
	//	cout << "ERROR: element type is wrong\n";
	//	exit(0);
	//case 5:
	//	return 3; // triangle
	//case 22:
	//	return 6; //6n_triangle
	//case 9:
	//	return 4; //quad
	//case 23:
	//	return 8; //8n_quad;
	//case 10:
	//	return 4; //4n_tetra
	//case 24:
	//	return 10;//10n_tetra
	//case 12:
	//	return 8;//hexahedron
	//case 25:
	//	return 20;// 20n_HEXAHEDRON
	//case 13:
	//	return 6; //wedge;
	//case 26:
	//	return 15; //15n_wedge
	//case 14:
	//	return 5;//pyramid;
	//case 27:
	//	return 13;//13n_PYRAMID;
	//}
	return ci_numNodes;
}

int pdfemEles::getAlgoType() const
{
	return ci_algorithmType;
}

int pdfemEles::getNumNodes() const
{
	return ci_numNodes;
}

double pdfemEles::detMat33(double mat[][3])
{
	//==det(J)===
	double det_J = 0.0;
	int rr, cc;
	double ans1, ans2;
	for (int i = 0; i < 3; i++) {
		rr = 0, cc = i;
		ans1 = 1., ans2 = 1.;
		for (int j = 0; j < 3; j++) {
			ans1 *= mat[rr + j][(cc + j) % 3];
			ans2 *= mat[rr + j][(cc - j + 3) % 3];
		}
		det_J = det_J + ans1 - ans2;
	}
	return (det_J);
}
