#pragma once

#include "Element.h"

//! H8 element class: 8-node 3D hexahedral solid element
class CH8 : public CElement
{
protected:
	virtual double HourglassAlphaBase() const;
	virtual double HourglassAlphaMin() const;
	virtual double GlobalYStiffnessScale() const;
	virtual bool UseSelectiveReducedIntegration() const;
	virtual double SRIHourglassBlend() const;
	virtual double FBHourglassScale() const;
	virtual bool UseOrthogonalHourglass() const;
	virtual bool UseFastReducedIntegration() const;
	virtual double OrthogonalHourglassScale() const;
	virtual double OrthogonalHourglassDirectionScale(int direction) const;

public:
	CH8();
	~CH8();

	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
	virtual void Write(COutputter& output);
	virtual void GenerateLocationMatrix();
	virtual void ElementStiffness(double* Matrix);
	virtual void ElementStress(double* stress, double* Displacement);
};

//! H8R element class: reduced-integration 8-node brick.
//! The implementation shares CH8's C3D8R-like one-point integration with
//! stabilization, while exposing a distinct input type for Abaqus C3D8R.
class CH8R : public CH8
{
protected:
	virtual bool UseFastReducedIntegration() const;
};

//! H8RPier element class: local C3D8R-like trial for pier solids.
class CH8RPier : public CH8R
{
protected:
	virtual double HourglassAlphaBase() const;
	virtual double HourglassAlphaMin() const;
	virtual double GlobalYStiffnessScale() const;
	virtual bool UseSelectiveReducedIntegration() const;
	virtual double SRIHourglassBlend() const;
	virtual double FBHourglassScale() const;
	virtual bool UseOrthogonalHourglass() const;
	virtual double OrthogonalHourglassScale() const;
	virtual double OrthogonalHourglassDirectionScale(int direction) const;
};
