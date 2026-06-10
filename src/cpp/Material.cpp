/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;
	Input >> E >> Nu >> Thickness;

	return true;
}

//	Read beam material data from stream Input
bool CBeamMaterial::Read(ifstream& Input)
{
	Input >> nset;
	Input >> E >> Area >> I;

	return true;
}

//	Write beam material data to Stream
void CBeamMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << setw(16) << I << endl;
}

//	Read H8 material data from stream Input
bool CH8Material::Read(ifstream& Input)
{
	Input >> nset;
	Input >> E >> Nu;

	return true;
}

//	Write H8 material data to Stream
void CH8Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << endl;
}

//	Read Beam3D material data from stream Input
bool CBeam3DMaterial::Read(ifstream& Input)
{
	Input >> nset;
	Input >> E >> Nu >> Area >> Iy >> Iz >> J;
	Input >> n1[0] >> n1[1] >> n1[2];

	return true;
}

//	Write Beam3D material data to Stream
void CBeam3DMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << Area
		   << setw(16) << Iy << setw(16) << Iz << setw(16) << J
		   << setw(16) << n1[0] << setw(16) << n1[1] << setw(16) << n1[2] << endl;
}

//	Write material data to Stream
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << Thickness << endl;
}
