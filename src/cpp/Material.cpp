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

//	Write material data to Stream
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Nu << setw(16) << Thickness << endl;
}
