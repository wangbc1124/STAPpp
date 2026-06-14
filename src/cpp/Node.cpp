/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;

    for (unsigned int i = 0; i < NDF; i++)
        bcode[i] = 0;	// Boundary codes
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
	string line;
	vector<double> values;
	while (std::getline(Input, line))
	{
		if (line.find_first_not_of(" \t\r\n") == string::npos)
			continue;

		istringstream iss(line);
		double value = 0.0;
		while (iss >> value)
			values.push_back(value);
		break;
	}

	if (values.size() == 7)
	{
		NodeNumber = static_cast<unsigned int>(values[0]);
		bcode[0] = static_cast<unsigned int>(values[1]);
		bcode[1] = static_cast<unsigned int>(values[2]);
		bcode[2] = static_cast<unsigned int>(values[3]);
		bcode[3] = 1;
		bcode[4] = 1;
		bcode[5] = 1;
		XYZ[0] = values[4];
		XYZ[1] = values[5];
		XYZ[2] = values[6];
		return true;
	}

	if (values.size() == 10)
	{
		NodeNumber = static_cast<unsigned int>(values[0]);
		for (unsigned int i = 0; i < NDF; i++)
			bcode[i] = static_cast<unsigned int>(values[i + 1]);
		XYZ[0] = values[7];
		XYZ[1] = values[8];
		XYZ[2] = values[9];
		return true;
	}

	return false;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(5) << bcode[3] << setw(5) << bcode[4] << setw(5) << bcode[5]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
