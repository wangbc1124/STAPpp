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

#include "Node.h"
#include "ElementGroup.h"
#include "Outputter.h"
#include "Solver.h"
#include "LoadCaseData.h"
#include "SkylineMatrix.h"
#include <vector>

using namespace std;

//!	Clear an array
template <class type> void clear( type* a, unsigned int N );

//!	Domain class : Define the problem domain
/*!	Only a single instance of Domain class can be created */
class CDomain
{
private:
	struct CDofAlias
	{
		unsigned int master_node = 0;
		unsigned int master_dof = 0;
	};

	struct CMPCEquationTerm
	{
		unsigned int node = 0;
		unsigned int dof = 0;
		double coefficient = 0.0;
	};

	struct CMPCEquation
	{
		std::vector<CMPCEquationTerm> terms;
		double rhs = 0.0;
	};

	struct CElementGroupHeader
	{
		bool pending = false;
		int element_type = 0;
		unsigned int nume = 0;
		unsigned int nummat = 0;
	};

	struct CSparseAssemblyTiming
	{
		double pattern_time = 0.0;
		double assembly_time = 0.0;
		double element_assembly_time = 0.0;
		double mpc_assembly_time = 0.0;
		double element_stiffness_time = 0.0;
		double pattern_insert_time = 0.0;
		double value_insert_time = 0.0;
		double active_dof_pack_time = 0.0;
		double mpc_pattern_time = 0.0;
		double export_upper_csr_time = 0.0;
		unsigned long long element_count = 0;
		unsigned long long element_value_insertions = 0;
		unsigned long long mpc_value_insertions = 0;
	};

//!	The instance of the Domain class
	static CDomain* _instance;

//!	Input file stream for reading data from input data file
	ifstream Input;

//!	Heading information for use in labeling the outpu
	char Title[256]; 

//!	Solution MODEX
/*!		0 : Data check only;
		1 : Execution */
	unsigned int MODEX;

//!	Total number of nodal points
	unsigned int NUMNP;

//!	List of all nodes in the domain
	CNode* NodeList;

//!	Total number of element groups.
/*! An element group consists of a convenient collection of elements with same type */
	unsigned int NUMEG;

//! Element group list
    CElementGroup* EleGrpList;
    
//!	Number of load cases
	unsigned int NLCASE;

//!	List of all load cases
	CLoadCaseData* LoadCases;

//!	Number of concentrated loads applied in each load case
	unsigned int* NLOAD;

//!	Total number of equations in the system
	unsigned int NEQ;

//! Optional equation aliases for tied translational DOFs
	std::vector<CDofAlias> DofAliases;

//! Multi-point constraint equations
	std::vector<CMPCEquation> ConstraintEquations;

//! Penalty factor used for MPC enforcement
	double ConstraintPenalty;

//! Time spent converting simple MPC equations to aliases during input
	double LastMpcAliasTime;

//! Time spent in the last sparse CSR assembly
	CSparseAssemblyTiming LastSparseAssemblyTiming;

//! First element group header read while probing an optional MPC section
	CElementGroupHeader PendingElementGroupHeader;

//!	Banded stiffness matrix
/*! A one-dimensional array storing only the elements below the	skyline of the 
    global stiffness matrix. */
    CSkylineMatrix<double>* StiffnessMatrix;

//!	Global nodal force/displacement vector
	double* Force;

private:

//!	Constructor
	CDomain();

//!	Desconstructor
	~CDomain();

public:

//!	Return pointer to the instance of the Domain class
	static CDomain* GetInstance();

//!	Read domain data from the input data file
	bool ReadData(string FileName, string OutFile);

//!	Read nodal point data
	bool ReadNodalPoints();

//!	Read optional DOF alias data inserted by the preprocessor
	bool ReadDofAliases();

//! Read optional MPC equations inserted by the preprocessor
	bool ReadConstraintEquations();

//! Convert simple two-term MPC equations into equation aliases
	void ConvertSimpleMPCToAliases();

//!	Read load case data
	bool ReadLoadCases();

//!	Read element data
	bool ReadElements();

//!	Calculate global equation numbers corresponding to every degree of freedom of each node
	void CalculateEquationNumber();

//!	Calculate column heights
	void CalculateColumnHeights();

//! Allocate storage for matrices
/*!	Allocate Force, ColumnHeights, DiagonalAddress and StiffnessMatrix and 
    calculate the column heights and address of diagonal elements */
	void AllocateMatrices();

//!	Assemble the banded gloabl stiffness matrix
	void AssembleStiffnessMatrix();

//!	Add penalty-form MPC contributions to the assembled stiffness matrix
	void ApplyConstraintEquations();

//! Assemble stiffness matrix into full CSR storage for iterative solvers
	CCSRMatrix* AssembleSparseStiffnessMatrix(const std::string& backend_name = "standard");

//! Assemble stiffness matrix into symmetric half storage for the PARDISO mainline
	CSparseSymmetricMatrix* AssemblePardisoStiffnessMatrix();

//! Add penalty-form MPC contributions to a sparse stiffness matrix
	void ApplySparseConstraintEquations(CCSRMatrix& matrix, unsigned long long* insertion_count = 0);

//!	Assemble the global nodal force vector for load case LoadCase
	bool AssembleForce(unsigned int LoadCase); 

//! Write nodal displacement results to CSV
	bool WriteDisplacementCSV(const string& FileName);

//! Resolve one nodal DOF displacement including alias backfill
	double GetResolvedDisplacement(unsigned int node_index, unsigned int dof_index) const;

//! Collect node-wise active equation blocks after BCs and aliasing
	std::vector<std::vector<unsigned int> > BuildNodeEquationBlocks() const;

//!	Return solution mode
	inline unsigned int GetMODEX() { return MODEX; }

//!	Return the title of problem
	inline string GetTitle() { return Title; }

//!	Return the total number of equations
	inline unsigned int GetNEQ() { return NEQ; }

//!	Return the total number of nodal points
	inline unsigned int GetNUMNP() { return NUMNP; }

//!	Return the node list
	inline CNode* GetNodeList() { return NodeList; }

//!	Return total number of element groups
	inline unsigned int GetNUMEG() { return NUMEG; }

//! Return element group list
    inline CElementGroup* GetEleGrpList() { return EleGrpList; }

//!	Return pointer to the global nodal force vector
	inline double* GetForce() { return Force; }

//!	Return pointer to the global nodal displacement vector
	inline double* GetDisplacement() { return Force; }

//!	Return the total number of load cases
	inline unsigned int GetNLCASE() { return NLCASE; }

//!	Return the number of concentrated loads applied in each load case
	inline unsigned int* GetNLOAD() { return NLOAD; }

//!	Return the list of load cases
	inline CLoadCaseData* GetLoadCases() { return LoadCases; }

//!	Return pointer to the banded stiffness matrix
	inline CSkylineMatrix<double>* GetStiffnessMatrix() { return StiffnessMatrix; }

//! Return time spent converting simple MPC equations during the last input pass
	inline double GetLastMpcAliasTime() const { return LastMpcAliasTime; }

//! Return time spent building the last CSR sparsity pattern
	inline double GetLastCsrPatternTime() const { return LastSparseAssemblyTiming.pattern_time; }

//! Return time spent filling the last CSR values
	inline double GetLastCsrAssemblyTime() const { return LastSparseAssemblyTiming.assembly_time; }
	inline double GetLastElementCsrAssemblyTime() const { return LastSparseAssemblyTiming.element_assembly_time; }
	inline double GetLastMpcCsrAssemblyTime() const { return LastSparseAssemblyTiming.mpc_assembly_time; }
	inline double GetLastElementStiffnessTime() const { return LastSparseAssemblyTiming.element_stiffness_time; }
	inline double GetLastPatternInsertTime() const { return LastSparseAssemblyTiming.pattern_insert_time; }
	inline double GetLastValueInsertTime() const { return LastSparseAssemblyTiming.value_insert_time; }
	inline double GetLastActiveDofPackTime() const { return LastSparseAssemblyTiming.active_dof_pack_time; }
	inline double GetLastMpcPatternTime() const { return LastSparseAssemblyTiming.mpc_pattern_time; }
	inline double GetLastExportUpperCsrTime() const { return LastSparseAssemblyTiming.export_upper_csr_time; }
	inline unsigned long long GetLastElementCount() const { return LastSparseAssemblyTiming.element_count; }
	inline unsigned long long GetLastElementValueInsertions() const { return LastSparseAssemblyTiming.element_value_insertions; }
	inline unsigned long long GetLastMpcValueInsertions() const { return LastSparseAssemblyTiming.mpc_value_insertions; }

};
