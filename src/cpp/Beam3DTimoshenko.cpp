/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Beam3DTimoshenko.h"
#include "Material.h"

#include <cmath>
#include <iomanip>

using namespace std;

namespace
{
	void BuildBeamTriad(const CNode* node1, const CNode* node2, const double n1[3], double R[3][3], double& L)
	{
		double dx = node2->XYZ[0] - node1->XYZ[0];
		double dy = node2->XYZ[1] - node1->XYZ[1];
		double dz = node2->XYZ[2] - node1->XYZ[2];
		L = sqrt(dx*dx + dy*dy + dz*dz);
		if (L < 1.0e-20)
			L = 1.0;

		double ex[3] = {dx / L, dy / L, dz / L};

		// Project user section direction to the plane normal to beam axis.
		double nref[3] = {n1[0], n1[1], n1[2]};
		double proj = nref[0] * ex[0] + nref[1] * ex[1] + nref[2] * ex[2];
		double ey[3] = {
			nref[0] - proj * ex[0],
			nref[1] - proj * ex[1],
			nref[2] - proj * ex[2]
		};

		double ey_norm = sqrt(ey[0] * ey[0] + ey[1] * ey[1] + ey[2] * ey[2]);
		if (ey_norm < 1.0e-10)
		{
			if (fabs(ex[2]) < 0.9)
			{
				ey[0] = -ex[1];
				ey[1] = ex[0];
				ey[2] = 0.0;
			}
			else
			{
				ey[0] = 0.0;
				ey[1] = -ex[2];
				ey[2] = ex[1];
			}
			ey_norm = sqrt(ey[0] * ey[0] + ey[1] * ey[1] + ey[2] * ey[2]);
		}
		for (int i = 0; i < 3; ++i)
			ey[i] /= ey_norm;

		double ez[3] = {
			ex[1] * ey[2] - ex[2] * ey[1],
			ex[2] * ey[0] - ex[0] * ey[2],
			ex[0] * ey[1] - ex[1] * ey[0]
		};

		R[0][0] = ex[0]; R[0][1] = ey[0]; R[0][2] = ez[0];
		R[1][0] = ex[1]; R[1][1] = ey[1]; R[1][2] = ez[1];
		R[2][0] = ex[2]; R[2][1] = ey[2]; R[2][2] = ez[2];
	}

	void BuildTransform(const double R[3][3], double T[12][12])
	{
		for (int i = 0; i < 12; ++i)
			for (int j = 0; j < 12; ++j)
				T[i][j] = 0.0;

		for (int n = 0; n < 2; ++n)
		{
			int base = 6 * n;
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
				{
					T[base + i][base + j] = R[j][i];
					T[base + 3 + i][base + 3 + j] = R[j][i];
				}
		}
	}

	void BuildLocalStiffness(const CBeam3DTimoshenkoMaterial* mat, double L, double k[12][12])
	{
		for (int i = 0; i < 12; ++i)
			for (int j = 0; j < 12; ++j)
				k[i][j] = 0.0;

		double E = mat->E;
		double G = E / (2.0 * (1.0 + mat->Nu));
		double A = mat->Area;
		double Iy = mat->Iy;
		double Iz = mat->Iz;
		double J = mat->J;
		double Asy = (mat->Asy > 1.0e-12) ? mat->Asy : A;
		double Asz = (mat->Asz > 1.0e-12) ? mat->Asz : A;

		double phi_z = 12.0 * E * Iz / (G * Asy * L * L);
		double phi_y = 12.0 * E * Iy / (G * Asz * L * L);
		double psi_z = 1.0 / (1.0 + phi_z);
		double psi_y = 1.0 / (1.0 + phi_y);

		double EAL = E * A / L;
		double GJL = G * J / L;

		k[0][0] = EAL;  k[0][6] = -EAL;
		k[6][0] = -EAL; k[6][6] = EAL;

		k[3][3] = GJL;  k[3][9] = -GJL;
		k[9][3] = -GJL; k[9][9] = GJL;

		double kv = 12.0 * E * Iz * psi_z / (L * L * L);
		double kvt = 6.0 * E * Iz * psi_z / (L * L);
		double ktt = (4.0 + phi_z) * E * Iz * psi_z / L;
		double ktt2 = (2.0 - phi_z) * E * Iz * psi_z / L;

		k[1][1]   = kv;    k[1][5]   = kvt;   k[1][7]   = -kv;   k[1][11]  = kvt;
		k[5][1]   = kvt;   k[5][5]   = ktt;   k[5][7]   = -kvt;  k[5][11]  = ktt2;
		k[7][1]   = -kv;   k[7][5]   = -kvt;  k[7][7]   = kv;    k[7][11]  = -kvt;
		k[11][1]  = kvt;   k[11][5]  = ktt2;  k[11][7]  = -kvt;  k[11][11] = ktt;

		double kw = 12.0 * E * Iy * psi_y / (L * L * L);
		double kwt = 6.0 * E * Iy * psi_y / (L * L);
		double kry = (4.0 + phi_y) * E * Iy * psi_y / L;
		double kry2 = (2.0 - phi_y) * E * Iy * psi_y / L;

		k[2][2]   = kw;    k[2][4]   = -kwt;  k[2][8]   = -kw;   k[2][10]  = -kwt;
		k[4][2]   = -kwt;  k[4][4]   = kry;   k[4][8]   = kwt;   k[4][10]  = kry2;
		k[8][2]   = -kw;   k[8][4]   = kwt;   k[8][8]   = kw;    k[8][10]  = kwt;
		k[10][2]  = -kwt;  k[10][4]  = kry2;  k[10][8]  = kwt;   k[10][10] = kry;
	}

	void Multiply12x12(const double A[12][12], const double B[12][12], double C[12][12])
	{
		for (int i = 0; i < 12; ++i)
			for (int j = 0; j < 12; ++j)
			{
				C[i][j] = 0.0;
				for (int k = 0; k < 12; ++k)
					C[i][j] += A[i][k] * B[k][j];
			}
	}
}

CBeam3DTimoshenko::CBeam3DTimoshenko()
{
	NEN_ = 2;
	nodes_ = new CNode*[NEN_];
	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];
	ElementMaterial_ = nullptr;
}

CBeam3DTimoshenko::~CBeam3DTimoshenko()
{
}

bool CBeam3DTimoshenko::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2;

	Input >> N1 >> N2 >> MSet;
	ElementMaterial_ = dynamic_cast<CBeam3DTimoshenkoMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

void CBeam3DTimoshenko::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CBeam3DTimoshenko::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int d = 0; d < 6; d++)
			LocationMatrix_[i++] = nodes_[N]->bcode[d];
}

void CBeam3DTimoshenko::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CBeam3DTimoshenkoMaterial* mat = dynamic_cast<CBeam3DTimoshenkoMaterial*>(ElementMaterial_);
	double R[3][3], T[12][12], k_local[12][12], TK[12][12], K_global[12][12], L;
	BuildBeamTriad(nodes_[0], nodes_[1], mat->n1, R, L);
	BuildTransform(R, T);
	BuildLocalStiffness(mat, L, k_local);

	for (int i = 0; i < 12; ++i)
		for (int j = 0; j < 12; ++j)
		{
			TK[i][j] = 0.0;
			for (int m = 0; m < 12; ++m)
				TK[i][j] += T[m][i] * k_local[m][j];
		}

	Multiply12x12(TK, T, K_global);

	for (int j = 0; j < 12; ++j)
		for (int i = 0; i <= j; ++i)
			Matrix[(j + 1) * j / 2 + j - i] = K_global[i][j];
}

void CBeam3DTimoshenko::ElementStress(double* stress, double* Displacement)
{
	for (int i = 0; i < 8; ++i)
		stress[i] = 0.0;

	CBeam3DTimoshenkoMaterial* mat = dynamic_cast<CBeam3DTimoshenkoMaterial*>(ElementMaterial_);
	double R[3][3], T[12][12], k_local[12][12], u_global[12], u_local[12], f_local[12], L;
	BuildBeamTriad(nodes_[0], nodes_[1], mat->n1, R, L);
	BuildTransform(R, T);
	BuildLocalStiffness(mat, L, k_local);

	for (int i = 0; i < 12; ++i)
	{
		u_global[i] = (LocationMatrix_[i] ? Displacement[LocationMatrix_[i] - 1] : 0.0);
		u_local[i] = 0.0;
		for (int j = 0; j < 12; ++j)
			u_local[i] += T[i][j] * u_global[j];
	}

	for (int i = 0; i < 12; ++i)
	{
		f_local[i] = 0.0;
		for (int j = 0; j < 12; ++j)
			f_local[i] += k_local[i][j] * u_local[j];
	}

	double axial_force = f_local[6];
	double axial_stress = axial_force / ((mat->Area > 1.0e-20) ? mat->Area : 1.0);

	stress[0] = axial_force;
	stress[1] = f_local[10];
	stress[2] = f_local[11];
	stress[3] = f_local[9];
	stress[4] = axial_stress;
	stress[5] = f_local[4];
	stress[6] = f_local[5];
	stress[7] = f_local[3];
}
