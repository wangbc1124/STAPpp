/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Material.h"
#include <cmath>
#include <functional>
#include <sstream>
#include <string>

using namespace std;

CDomain* CDomain::_instance = nullptr;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = nullptr;
	DofAliases.clear();
	ConstraintEquations.clear();
	ConstraintPenalty = 0.0;
	
	NUMEG = 0;
	EleGrpList = nullptr;
	
	NLCASE = 0;
	NLOAD = nullptr;
	LoadCases = nullptr;
	
	NEQ = 0;

	Force = nullptr;
	StiffnessMatrix = nullptr;
}

//	Desconstructor
CDomain::~CDomain()
{
	delete [] NodeList;

	delete [] EleGrpList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] Force;
	delete StiffnessMatrix;
}

//	Return pointer to the instance of the Domain class
CDomain* CDomain::GetInstance()
{
	if (!_instance) 
		_instance = new CDomain();
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	COutputter* Output = COutputter::GetInstance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	if (ReadNodalPoints())
        Output->OutputNodeInfo();
    else
        return false;

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
	if (ReadLoadCases())
        Output->OutputLoadInfo();
    else
        return false;

//	Read optional MPC equations
	if (!ReadConstraintEquations())
		return false;

//	Read element data
	if (ReadElements())
        Output->OutputElementInfo();
    else
        return false;

	return true;
}

bool CDomain::ReadDofAliases()
{
	DofAliases.assign(NUMNP * CNode::NDF, CDofAlias());
	return true;
}

bool CDomain::ReadConstraintEquations()
{
	Input >> ws;
	streampos pos = Input.tellg();
	string tag;
	if (!(Input >> tag))
	{
		Input.clear();
		return true;
	}
	if (tag != "MPC")
	{
		Input.clear();
		Input.seekg(pos);
		return true;
	}

	unsigned int nEquation = 0;
	Input >> nEquation;
	if (!Input)
	{
		cerr << "*** Error *** Invalid MPC section header." << endl;
		return false;
	}

	ConstraintEquations.clear();
	ConstraintEquations.reserve(nEquation);
	for (unsigned int i = 0; i < nEquation; ++i)
	{
		unsigned int nTerm = 0;
		double rhs = 0.0;
		Input >> nTerm >> rhs;
		if (!Input)
		{
			cerr << "*** Error *** Invalid MPC header record in input file." << endl;
			return false;
		}

		CMPCEquation equation;
		equation.rhs = rhs;
		equation.terms.reserve(nTerm);
		for (unsigned int t = 0; t < nTerm; ++t)
		{
			CMPCEquationTerm term;
			Input >> term.node >> term.dof >> term.coefficient;
			if (!Input)
			{
				cerr << "*** Error *** Invalid MPC term record in input file." << endl;
				return false;
			}
			if (term.node < 1 || term.node > NUMNP || term.dof < 1 || term.dof > CNode::NDF)
			{
				cerr << "*** Error *** MPC term out of range: "
					 << term.node << ' ' << term.dof << endl;
				return false;
			}
			equation.terms.push_back(term);
		}

		ConstraintEquations.push_back(equation);
	}

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new CNode[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
    {
		if (!NodeList[np].Read(Input))
			return false;
    
        if (NodeList[np].NodeNumber != np + 1)
        {
            cerr << "*** Error *** Nodes must be inputted in order !" << endl
            << "   Expected node number : " << np + 1 << endl
            << "   Provided node number : " << NodeList[np].NodeNumber << endl;
        
            return false;
        }
    }

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	if (DofAliases.size() != NUMNP * CNode::NDF)
		DofAliases.assign(NUMNP * CNode::NDF, CDofAlias());

	vector<unsigned int> equation(NUMNP * CNode::NDF, 0);
	vector<unsigned char> state(NUMNP * CNode::NDF, 0);

	std::function<unsigned int(unsigned int, unsigned int)> assign_eq =
		[&](unsigned int node, unsigned int dof) -> unsigned int
	{
		const unsigned int idx = node * CNode::NDF + dof;
		if (NodeList[node].bcode[dof])
			return 0;
		if (equation[idx])
			return equation[idx];
		if (state[idx] == 1)
		{
			cerr << "*** Error *** Cyclic ALIAS definition detected at node "
				 << (node + 1) << ", dof " << (dof + 1) << endl;
			exit(4);
		}

		state[idx] = 1;
		const CDofAlias& alias = DofAliases[idx];
		if (alias.master_node)
		{
			equation[idx] = assign_eq(alias.master_node - 1, alias.master_dof - 1);
		}
		else
		{
			equation[idx] = ++NEQ;
		}
		state[idx] = 2;
		return equation[idx];
	};

	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
		{
			NodeList[np].bcode[dof] = assign_eq(np, dof);
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new CLoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
    {
        unsigned int LL;
        Input >> LL;
        
        if (LL != lcase + 1)
        {
            cerr << "*** Error *** Load case must be inputted in order !" << endl
            << "   Expected load case : " << lcase + 1 << endl
            << "   Provided load case : " << LL << endl;
            
            return false;
        }

        LoadCases[lcase].Read(Input);
    }

	return true;
}

// Read element data
bool CDomain::ReadElements()
{
    EleGrpList = new CElementGroup[NUMEG];

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
        if (!EleGrpList[EleGrp].Read(Input))
            return false;
    
    return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
#ifdef _DEBUG_
    COutputter* Output = COutputter::GetInstance();
    *Output << setw(9) << "Ele = " << setw(22) << "Location Matrix" << endl;
#endif

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
    {
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();
        
		for (unsigned int Ele = 0; Ele < NUME; Ele++)	//	Loop over for all elements in group EleGrp
        {
            CElement& Element = ElementGrp[Ele];

            // Generate location matrix
            Element.GenerateLocationMatrix();
            
#ifdef _DEBUG_
            unsigned int* LocationMatrix = Element.GetLocationMatrix();
            
            *Output << setw(9) << Ele+1;
            for (int i=0; i<Element.GetND(); i++)
                *Output << setw(5) << LocationMatrix[i];
            *Output << endl;
#endif

            StiffnessMatrix->CalculateColumnHeight(Element.GetLocationMatrix(), Element.GetND());
        }
    }

	for (size_t eq = 0; eq < ConstraintEquations.size(); ++eq)
	{
		vector<unsigned int> location;
		location.reserve(ConstraintEquations[eq].terms.size());
		for (size_t term = 0; term < ConstraintEquations[eq].terms.size(); ++term)
		{
			const CMPCEquationTerm& item = ConstraintEquations[eq].terms[term];
			unsigned int equationNo = NodeList[item.node - 1].bcode[item.dof - 1];
			if (equationNo)
				location.push_back(equationNo);
		}
		if (!location.empty())
			StiffnessMatrix->CalculateColumnHeight(location.data(), location.size());
	}
    
    StiffnessMatrix->CalculateMaximumHalfBandwidth();
    
#ifdef _DEBUG_
    *Output << endl;
	Output->PrintColumnHeights();
#endif

}

//    Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//    and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
    //    Allocate for global force/displacement vector
    Force = new double[NEQ];
    
    //  Create the banded stiffness matrix
    StiffnessMatrix = new CSkylineMatrix<double>(NEQ);
    
    //    Calculate column heights
    CalculateColumnHeights();
    
    //    Calculate address of diagonal elements in banded matrix
    StiffnessMatrix->CalculateDiagnoalAddress();
    
    //    Allocate for banded global stiffness matrix
    StiffnessMatrix->Allocate();
    
    COutputter* Output = COutputter::GetInstance();
    Output->OutputTotalSystemData();
}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
        CElementGroup& ElementGrp = EleGrpList[EleGrp];
        unsigned int NUME = ElementGrp.GetNUME();

		unsigned int size = ElementGrp[0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME; Ele++)
        {
            CElement& Element = ElementGrp[Ele];
            Element.ElementStiffness(Matrix);
            StiffnessMatrix->Assembly(Matrix, Element.GetLocationMatrix(), Element.GetND());
        }

		delete[] Matrix;
		Matrix = nullptr;
	}

	ApplyConstraintEquations();

#ifdef _DEBUG_
	COutputter* Output = COutputter::GetInstance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];

    clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		unsigned int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
        
        if(dof) // The DOF is activated
            Force[dof - 1] += LoadData->load[lnum];
	}

	if (ConstraintPenalty > 0.0)
	{
		for (size_t eq = 0; eq < ConstraintEquations.size(); ++eq)
		{
			const CMPCEquation& equation = ConstraintEquations[eq];
			if (fabs(equation.rhs) <= 0.0)
				continue;
			for (size_t term = 0; term < equation.terms.size(); ++term)
			{
				const CMPCEquationTerm& item = equation.terms[term];
				unsigned int equationNo = NodeList[item.node - 1].bcode[item.dof - 1];
				if (equationNo)
					Force[equationNo - 1] += ConstraintPenalty * item.coefficient * equation.rhs;
			}
		}
	}

	return true;
}

void CDomain::ApplyConstraintEquations()
{
	if (ConstraintEquations.empty())
		return;

	double max_diag = 0.0;
	for (unsigned int i = 1; i <= NEQ; ++i)
	{
		double diag = fabs((*StiffnessMatrix)(i, i));
		if (diag > max_diag)
			max_diag = diag;
	}
	if (max_diag <= 0.0)
		max_diag = 1.0;

	ConstraintPenalty = max_diag * 1.0e8;

	for (size_t eq = 0; eq < ConstraintEquations.size(); ++eq)
	{
		vector<unsigned int> equations;
		vector<double> coefficients;

		for (size_t term = 0; term < ConstraintEquations[eq].terms.size(); ++term)
		{
			const CMPCEquationTerm& item = ConstraintEquations[eq].terms[term];
			unsigned int equationNo = NodeList[item.node - 1].bcode[item.dof - 1];
			if (!equationNo)
				continue;

			bool merged = false;
			for (size_t i = 0; i < equations.size(); ++i)
			{
				if (equations[i] == equationNo)
				{
					coefficients[i] += item.coefficient;
					merged = true;
					break;
				}
			}
			if (!merged)
			{
				equations.push_back(equationNo);
				coefficients.push_back(item.coefficient);
			}
		}

		for (size_t j = 0; j < equations.size(); ++j)
		{
			if (fabs(coefficients[j]) < 1.0e-20)
				continue;
			for (size_t i = 0; i <= j; ++i)
			{
				if (fabs(coefficients[i]) < 1.0e-20)
					continue;
				(*StiffnessMatrix)(equations[i], equations[j]) +=
					ConstraintPenalty * coefficients[i] * coefficients[j];
			}
		}
	}
}

