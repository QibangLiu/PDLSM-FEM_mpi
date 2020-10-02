#pragma once
#include<fstream>
#include<iomanip>
#include"Matrix.h"
#include"Vector.h"
#include"calMatrixOperations.h"
#include"pdGaussPt.h"
using namespace std;

/*This class is the parent of elements
The element type ci_eleType is compatiable with vtk library (see vtkCell.h).
No.	element			vtk_name				ci_eleType
1	3n_triangle		VTK_TRIANGLE			(=5) 
2	6n_triangle		VTK_QUADRATIC_TRIANGLE	(=22)
3	4n_quad			VTK_QUAD				(=9)
4	8n_quad			VTK_QUADRATIC_QUAD		(=23)
5	4n_tetra		VTK_TETRA				(=10)
6	10n_tetra		VTK_QUADRATIC_TETRA		(=24)
7	8n_HEXAHEDRON	VTK_HEXAHEDRON			(=12)
8	20n_HEXAHEDRON	VTK_QUADRATIC_HEXAHEDRON(=25)
9	6n_WEDGE		VTK_WEDGE				(=13)
10	15n_wedge		VTK_QUADRATIC_WEDGE		(=26)
11	5n_PYRAMID		VTK_PYRAMID				(=14)
12	13n_PYRAMID		VTK_QUADRATIC_PYRAMID	(=27) */



class pdfemEles
{
public:
	pdfemEles(int id, int numNodes, int *nId, int algoType);
	virtual ~pdfemEles();
	void getConNid(int Nid[]);
	void print(ofstream &fout);
	int getNumNodes_vtk()const;
	int getAlgoType()const;
	int getNumNodes()const;
	virtual double detJacobi(double xN[][3], double p, double q, double r) = 0;
	virtual void shapeFunction(double N[], double p, double q, double r) = 0;
	virtual void eleStiffMatFEM(Matrix* Ke, Matrix *D, double xN[][3]) = 0;
	virtual void eleMassMat(Matrix* Me,double rho, double xN[][3]) = 0;
	virtual void eleEquivNodalForce(Vector* Fe, double t, double xN[][3]) = 0;
	virtual void eleFitStresses(int flag, Vector* Nsigma[], Matrix* D, Matrix* L, Vector* Ue, double xN[][3]) = 0;
	virtual void print_vtk(ofstream& fout,int *eleNodeID) = 0;

public:
	int ci_numNodes;
	int* cip_EleNodeId;
	int ci_EleId;
	int ci_algorithmType;
	// algorithm type, 1 is for PD, 2 is traditional FEM;
	//-1 is NBCs elements
	int ci_eleType;

private:
	pdfemEles();
	
};


