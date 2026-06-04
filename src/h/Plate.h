#pragma once

#include "Element.h"

//! Plate element class: 4-node Mindlin-Reissner plate bending element
class CPlate : public CElement
{
public:
	CPlate();
	~CPlate();

	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
	virtual void Write(COutputter& output);
	virtual void GenerateLocationMatrix();
	virtual void ElementStiffness(double* Matrix);
	virtual void ElementStress(double* stress, double* Displacement);
};
