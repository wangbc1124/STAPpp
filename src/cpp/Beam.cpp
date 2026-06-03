#include "Beam.h"

#include <cmath>
#include <iomanip>

using namespace std;

CBeam::CBeam()
{
	NEN_ = 2;
	nodes_ = new CNode*[NEN_];

	ND_ = 6;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

CBeam::~CBeam()
{
}

bool CBeam::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2;

	Input >> N1 >> N2 >> MSet;
	ElementMaterial_ = dynamic_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

void CBeam::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CBeam::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		LocationMatrix_[i++] = nodes_[N]->bcode[0];
		LocationMatrix_[i++] = nodes_[N]->bcode[1];
		LocationMatrix_[i++] = nodes_[N]->bcode[2];
	}
}

void CBeam::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CBeamMaterial* material = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
	double E = material->E;
	double A = material->Area;
	double I = material->I;

	double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1];

	double dx = x2 - x1;
	double dy = y2 - y1;
	double L = sqrt(dx * dx + dy * dy);

	double c = dx / L;
	double s = dy / L;

	// Local stiffness matrix (6x6) for DOFs [u1, v1, theta1, u2, v2, theta2]
	double k[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			k[i][j] = 0.0;

	double EAL = E * A / L;
	double EIL3 = E * I / (L * L * L);

	// Axial stiffness
	k[0][0] =  EAL;
	k[0][3] = -EAL;
	k[3][0] = -EAL;
	k[3][3] =  EAL;

	// Bending stiffness
	k[1][1] =  12.0 * EIL3;
	k[1][2] =  6.0 * L * EIL3;
	k[1][4] = -12.0 * EIL3;
	k[1][5] =  6.0 * L * EIL3;

	k[2][1] =  6.0 * L * EIL3;
	k[2][2] =  4.0 * L * L * EIL3;
	k[2][4] = -6.0 * L * EIL3;
	k[2][5] =  2.0 * L * L * EIL3;

	k[4][1] = -12.0 * EIL3;
	k[4][2] = -6.0 * L * EIL3;
	k[4][4] =  12.0 * EIL3;
	k[4][5] = -6.0 * L * EIL3;

	k[5][1] =  6.0 * L * EIL3;
	k[5][2] =  2.0 * L * L * EIL3;
	k[5][4] = -6.0 * L * EIL3;
	k[5][5] =  4.0 * L * L * EIL3;

	// Rotation matrix: T = [[r, 0], [0, r]] where r is 3x3 per node
	double T[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			T[i][j] = 0.0;

	T[0][0] =  c;  T[0][1] = s;
	T[1][0] = -s;  T[1][1] = c;
	T[2][2] = 1.0;

	T[3][3] =  c;  T[3][4] = s;
	T[4][3] = -s;  T[4][4] = c;
	T[5][5] = 1.0;

	// K_global = T^T * K_local * T
	double K[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			K[i][j] = 0.0;

	// First compute K_local * T
	double KT[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
		{
			KT[i][j] = 0.0;
			for (int m = 0; m < 6; m++)
				KT[i][j] += k[i][m] * T[m][j];
		}

	// Then K = T^T * KT
	for (int i = 0; i < 6; i++)
		for (int j = i; j < 6; j++)
		{
			double val = 0.0;
			for (int m = 0; m < 6; m++)
				val += T[m][i] * KT[m][j];
			K[i][j] = val;
		}

	// Copy upper-triangular K to column-by-column storage
	for (int j = 0; j < 6; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CBeam::ElementStress(double* stress, double* Displacement)
{
	CBeamMaterial* material = dynamic_cast<CBeamMaterial*>(ElementMaterial_);
	double E = material->E;
	double A = material->Area;
	double I = material->I;

	double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1];

	double dx = x2 - x1;
	double dy = y2 - y1;
	double L = sqrt(dx * dx + dy * dy);

	double c = dx / L;
	double s = dy / L;

	// Extract global displacements
	double u_glob[6];
	for (int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			u_glob[i] = Displacement[LocationMatrix_[i] - 1];
		else
			u_glob[i] = 0.0;
	}

	// Transform to local coordinates: u_local = T * u_global
	double t1 =  c * u_glob[0] + s * u_glob[1];
	double t2 = -s * u_glob[0] + c * u_glob[1];
	double r1 = u_glob[2];
	double t3 =  c * u_glob[3] + s * u_glob[4];
	double t4 = -s * u_glob[3] + c * u_glob[4];
	double r2 = u_glob[5];

	// Axial stress
	double axial_stress = E * (t3 - t1) / L;

	// Bending moments at nodes 1 and 2 (from beam stiffness matrix)
	double EIL2 = E * I / (L * L);
	double M1 = EIL2 * ( 6.0 * t2 + 4.0 * L * r1 - 6.0 * t4 + 2.0 * L * r2);
	double M2 = EIL2 * ( 6.0 * t2 + 2.0 * L * r1 - 6.0 * t4 + 4.0 * L * r2);

	// Max bending stress at extreme fiber (y_max depends on section)
	// For rectangular section: y_max = h/2, I = b*h^3/12
	// stress_bend = M * y_max / I
	// We output M directly (user can convert to stress)
	double y_max = 1.0; // Default, user scales by their section
	double bend_stress1 = M1 * y_max / I;
	double bend_stress2 = M2 * y_max / I;

	stress[0] = axial_stress;
	stress[1] = bend_stress1;
	stress[2] = bend_stress2;
}
