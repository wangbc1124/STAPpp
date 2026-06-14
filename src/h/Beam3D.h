/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//!	Beam3D element class (3D Euler-Bernoulli beam: axial + torsion + biaxial bending)
class CBeam3D : public CElement
{
public:

//!	Constructor
	CBeam3D();

//!	Desconstructor
	~CBeam3D();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Generate location matrix (all 6 DOFs per node)
	virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix (12x12 upper triangular, col-by-col)
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);
};
