#pragma once

#include "Element.h"

//! H8 element class: 8-node 3D hexahedral solid element
class CH8 : public CElement
{
public:
	CH8();
	~CH8();

	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
	virtual void Write(COutputter& output);
	virtual void GenerateLocationMatrix();
	virtual void ElementStiffness(double* Matrix);
	virtual void ElementStress(double* stress, double* Displacement);
};
