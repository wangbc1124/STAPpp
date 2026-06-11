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

//! Beam3DTimoshenko element class (3D Timoshenko beam: axial + torsion + biaxial bending + shear)
class CBeam3DTimoshenko : public CElement
{
public:

//!	Constructor
	CBeam3DTimoshenko();

//!	Desconstructor
	~CBeam3DTimoshenko();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//!	Write element data to stream
	virtual void Write(COutputter& output);

//!	Generate location matrix (all 6 DOFs per node)
	virtual void GenerateLocationMatrix();

//!	Calculate element stiffness matrix (12x12 upper triangular, col-by-col)
	virtual void ElementStiffness(double* Matrix);

//!	Calculate element end actions / stress resultants
	virtual void ElementStress(double* stress, double* Displacement);
};
