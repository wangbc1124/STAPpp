/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Beam3D.h"
#include "Material.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBeam3D::CBeam3D()
{
	NEN_ = 2;	// Each element has 2 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 12;   // 6 DOF/node x 2 nodes
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CBeam3D::~CBeam3D()
{
}

//	Read element data from stream Input
bool CBeam3D::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2;

	Input >> N1 >> N2 >> MSet;
	ElementMaterial_ = dynamic_cast<CBeam3DMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream
void CBeam3D::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Generate location matrix (all 6 DOFs per node)
void CBeam3D::GenerateLocationMatrix()
{
	// DOF order per node: dx, dy, dz, rx, ry, rz
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		for (unsigned int d = 0; d < 6; d++)
			LocationMatrix_[i++] = nodes_[N]->bcode[d];
	}
}

//	Calculate element stiffness matrix (12x12, upper triangular, col-by-col)
void CBeam3D::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CBeam3DMaterial* mat = dynamic_cast<CBeam3DMaterial*>(ElementMaterial_);
	double E  = mat->E;
	double Nu = mat->Nu;
	double A  = mat->Area;
	double Iy = mat->Iy;
	double Iz = mat->Iz;
	double J  = mat->J;
	double G  = E / (2.0 * (1.0 + Nu));

	double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1], z1 = nodes_[0]->XYZ[2];
	double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1], z2 = nodes_[1]->XYZ[2];

	// --- Element length and direction ---
	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;
	double L = sqrt(dx*dx + dy*dy + dz*dz);

	if (L < 1e-20) L = 1.0;  // avoid division by zero

	// Local x' axis = unit vector along element (node1 -> node2)
	double ex[3] = {dx/L, dy/L, dz/L};

	// Local y' axis: user n1 direction crossed with x', normalized
	double n1x = mat->n1[0], n1y = mat->n1[1], n1z = mat->n1[2];
	double n1_len = sqrt(n1x*n1x + n1y*n1y + n1z*n1z);

	// Cross product: ey' = n1 x ex'
	double ey_tmp[3];
	ey_tmp[0] = n1y * ex[2] - n1z * ex[1];
	ey_tmp[1] = n1z * ex[0] - n1x * ex[2];
	ey_tmp[2] = n1x * ex[1] - n1y * ex[0];
	double ey_len = sqrt(ey_tmp[0]*ey_tmp[0] + ey_tmp[1]*ey_tmp[1] + ey_tmp[2]*ey_tmp[2]);

	double ey[3];
	if (ey_len > 1e-10) {
		ey[0] = ey_tmp[0] / ey_len;
		ey[1] = ey_tmp[1] / ey_len;
		ey[2] = ey_tmp[2] / ey_len;
	} else {
		// n1 is parallel to ex — fallback using global Z or Y
		if (fabs(ex[2]) < 0.999) {
			// Use global Z: ey' = Z x ex'
			ey_tmp[0] = -ex[1];
			ey_tmp[1] =  ex[0];
			ey_tmp[2] =  0.0;
		} else {
			// ex nearly parallel to Z, use global Y: ey' = ex x Y
			ey_tmp[0] = -ex[2];
			ey_tmp[1] =  0.0;
			ey_tmp[2] =  ex[0];
		}
		ey_len = sqrt(ey_tmp[0]*ey_tmp[0] + ey_tmp[1]*ey_tmp[1] + ey_tmp[2]*ey_tmp[2]);
		ey[0] = ey_tmp[0] / ey_len;
		ey[1] = ey_tmp[1] / ey_len;
		ey[2] = ey_tmp[2] / ey_len;
	}

	// Local z' axis: ez' = ex' x ey'
	double ez[3];
	ez[0] = ex[1] * ey[2] - ex[2] * ey[1];
	ez[1] = ex[2] * ey[0] - ex[0] * ey[2];
	ez[2] = ex[0] * ey[1] - ex[1] * ey[0];

	// --- Build rotation matrix R(3x3): columns are local axes in global coords ---
	// R * v_local = v_global,   T * [dx,dy,dz,rx,ry,rz]_local^T = [...]_global^T
	double R[3][3] = {
		{ex[0], ey[0], ez[0]},
		{ex[1], ey[1], ez[1]},
		{ex[2], ey[2], ez[2]}
	};

	// --- Build local stiffness matrix k(12x12) ---
	double k[12][12];
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
			k[i][j] = 0.0;

	// DOF order in local coordinates per node: u=axial(x'), v=lateral(y'), w=lateral(z'), tx=torsion(x'), ty=bending(y'), tz=bending(z')
	// Node 1: indices 0-5, Node 2: indices 6-11

	double EAL = E * A / L;
	double GJL = G * J / L;
	double EIy_L3 = E * Iy / (L*L*L);    // bending in XZ plane (about y')
	double EIz_L3 = E * Iz / (L*L*L);    // bending in XY plane (about z')

	// Axial stiffness: u1(0), u2(6)
	k[0][0] =  EAL;  k[0][6] = -EAL;
	k[6][0] = -EAL;  k[6][6] =  EAL;

	// Torsional stiffness: tx1(3), tx2(9)
	k[3][3] =  GJL;  k[3][9] = -GJL;
	k[9][3] = -GJL;  k[9][9] =  GJL;

	// Bending in XY plane (about z'): v(1), tz(5), v(7), tz(11) — uses Iz
	k[1][1]  =  12.0 * EIz_L3;
	k[1][5]  =   6.0 * L * EIz_L3;
	k[1][7]  = -12.0 * EIz_L3;
	k[1][11] =   6.0 * L * EIz_L3;

	k[5][1]  =   6.0 * L * EIz_L3;
	k[5][5]  =   4.0 * L*L * EIz_L3;
	k[5][7]  =  -6.0 * L * EIz_L3;
	k[5][11] =   2.0 * L*L * EIz_L3;

	k[7][1]  = -12.0 * EIz_L3;
	k[7][5]  =  -6.0 * L * EIz_L3;
	k[7][7]  =  12.0 * EIz_L3;
	k[7][11] =  -6.0 * L * EIz_L3;

	k[11][1]  =   6.0 * L * EIz_L3;
	k[11][5]  =   2.0 * L*L * EIz_L3;
	k[11][7]  =  -6.0 * L * EIz_L3;
	k[11][11] =   4.0 * L*L * EIz_L3;

	// Bending in XZ plane (about y'): w(2), ty(4), w(8), ty(10) — uses Iy
	// Sign convention: consistent with right-hand rule about y'
	k[2][2]   =  12.0 * EIy_L3;
	k[2][4]   =  -6.0 * L * EIy_L3;   // note sign difference
	k[2][8]   = -12.0 * EIy_L3;
	k[2][10]  =  -6.0 * L * EIy_L3;

	k[4][2]   =  -6.0 * L * EIy_L3;
	k[4][4]   =   4.0 * L*L * EIy_L3;
	k[4][8]   =   6.0 * L * EIy_L3;
	k[4][10]  =   2.0 * L*L * EIy_L3;

	k[8][2]   = -12.0 * EIy_L3;
	k[8][4]   =   6.0 * L * EIy_L3;
	k[8][8]   =  12.0 * EIy_L3;
	k[8][10]  =   6.0 * L * EIy_L3;

	k[10][2]  =  -6.0 * L * EIy_L3;
	k[10][4]  =   2.0 * L*L * EIy_L3;
	k[10][8]  =   6.0 * L * EIy_L3;
	k[10][10] =   4.0 * L*L * EIy_L3;

	// --- Transform: K_global = T^T * K_local * T ---
	// T(12x12) = block diagonal of R(3x3) for translation and R(3x3) for rotation at each node
	// T maps: [global_xyz, global_rxyz]_local^T → each node applies R to both translation and rotation blocks

	// Transform: K_global = T^T * K_local * T

	double T[12][12];
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
			T[i][j] = 0.0;

	// T maps global → local: u_local = T * u_global
	// R maps local → global, so T uses R^T = R[j][i]
	// T = block diagonal: each node block is [R^T   0 ]
	//                                        [ 0   R^T]
	for (int n = 0; n < 2; n++) {
		int base = n * 6;
		// Translation block
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				T[base + i][base + j] = R[j][i];
		// Rotation block
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				T[base + 3 + i][base + 3 + j] = R[j][i];
	}

	// Step 1: TK = T^T * k  (T^T is just transpose since R is orthonormal → T^T = T_transpose)
	double TK[12][12];
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++) {
			TK[i][j] = 0.0;
			for (int m = 0; m < 12; m++)
				TK[i][j] += T[m][i] * k[m][j];   // T^T[i][m] = T[m][i]
		}

	// Step 2: K_global = TK * T
	double K_global[12][12];
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++) {
			K_global[i][j] = 0.0;
			for (int m = 0; m < 12; m++)
				K_global[i][j] += TK[i][m] * T[m][j];
		}

	// Copy upper-triangular K_global to column-by-column skyline storage
	// Format: for each column j, store from diagonal element to top of column
	for (int j = 0; j < 12; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K_global[i][j];
}

//	Calculate element stress
void CBeam3D::ElementStress(double* stress, double* Displacement)
{
	CBeam3DMaterial* mat = dynamic_cast<CBeam3DMaterial*>(ElementMaterial_);
	double E  = mat->E;
	double Nu = mat->Nu;
	double A  = mat->Area;
	double Iy = mat->Iy;
	double Iz = mat->Iz;
	double J  = mat->J;
	double G  = E / (2.0 * (1.0 + Nu));

	double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1], z1 = nodes_[0]->XYZ[2];
	double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1], z2 = nodes_[1]->XYZ[2];

	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;
	double L = sqrt(dx*dx + dy*dy + dz*dz);
	if (L < 1e-20) L = 1.0;

	// Recompute local axes (same as stiffness)
	double ex[3] = {dx/L, dy/L, dz/L};

	double n1x = mat->n1[0], n1y = mat->n1[1], n1z = mat->n1[2];
	double ey_tmp[3];
	ey_tmp[0] = n1y * ex[2] - n1z * ex[1];
	ey_tmp[1] = n1z * ex[0] - n1x * ex[2];
	ey_tmp[2] = n1x * ex[1] - n1y * ex[0];
	double ey_len = sqrt(ey_tmp[0]*ey_tmp[0] + ey_tmp[1]*ey_tmp[1] + ey_tmp[2]*ey_tmp[2]);

	double ey[3];
	if (ey_len > 1e-10) {
		for (int i = 0; i < 3; i++) ey[i] = ey_tmp[i] / ey_len;
	} else {
		if (fabs(ex[2]) < 0.999) {
			ey[0] = -ex[1]; ey[1] =  ex[0]; ey[2] = 0.0;
		} else {
			ey[0] = -ex[2]; ey[1] =  0.0;   ey[2] = ex[0];
		}
		ey_len = sqrt(ey[0]*ey[0] + ey[1]*ey[1] + ey[2]*ey[2]);
		for (int i = 0; i < 3; i++) ey[i] /= ey_len;
	}

	// Rotation matrix R
	double R[3][3] = {
		{ex[0], ey[0], ex[1]*ey[2] - ex[2]*ey[1]},
		{ex[1], ey[1], ex[2]*ey[0] - ex[0]*ey[2]},
		{ex[2], ey[2], ex[0]*ey[1] - ex[1]*ey[0]}
	};

	// Extract global displacements
	double u_glob[12];
	for (int i = 0; i < 12; i++) {
		if (LocationMatrix_[i])
			u_glob[i] = Displacement[LocationMatrix_[i] - 1];
		else
			u_glob[i] = 0.0;
	}

	// Transform to local: u_local = R^T * u_global
	double u_loc[12];
	for (int i = 0; i < 12; i++) {
		u_loc[i] = 0.0;
		int n = i / 6;
		int d = i % 6;
		int base = n * 6;
		if (d < 3) {
			for (int m = 0; m < 3; m++)
				u_loc[i] += R[m][d] * u_glob[base + m];
		} else {
			for (int m = 0; m < 3; m++)
				u_loc[i] += R[m][d-3] * u_glob[base + 3 + m];
		}
	}

	// DOF in local coords:
	// u_loc[0-5]:  u1, v1, w1, tx1, ty1, tz1
	// u_loc[6-11]: u2, v2, w2, tx2, ty2, tz2

	// Axial force: N = EA * (u2 - u1) / L
	double N_axial = E * A * (u_loc[6] - u_loc[0]) / L;

	// Torsion: T = GJ * (tx2 - tx1) / L
	double T_torque = G * J * (u_loc[9] - u_loc[3]) / L;

	// Bending moments at node 1 (about z'): Mz1 = EIz * (6*(v1-v2)/L² + 4*L*tz1 + 2*L*tz2) / L
	// ... at node 2
	// Using beam stiffness equations:
	//   M1 = EI/L² * (6v1 + 4L*θ1 - 6v2 + 2L*θ2)
	//   M2 = EI/L² * (6v1 + 2L*θ1 - 6v2 + 4L*θ2)
	double Mz1 = E * Iz / (L*L) * (6.0*u_loc[1] + 4.0*L*u_loc[5] - 6.0*u_loc[7] + 2.0*L*u_loc[11]);
	double Mz2 = E * Iz / (L*L) * (6.0*u_loc[1] + 2.0*L*u_loc[5] - 6.0*u_loc[7] + 4.0*L*u_loc[11]);

	// Bending moments about y' (XZ plane): similar but with sign for w,ty
	double My1 = E * Iy / (L*L) * (-6.0*u_loc[2] + 4.0*L*u_loc[4] + 6.0*u_loc[8] + 2.0*L*u_loc[10]);
	double My2 = E * Iy / (L*L) * (-6.0*u_loc[2] + 2.0*L*u_loc[4] + 6.0*u_loc[8] + 4.0*L*u_loc[10]);

	// Axial stress: sigma = N / A
	double axial_stress = N_axial / A;

	// Output: [N, My1, Mz1, T, sigma, My2, Mz2, T]  (8 values)
	stress[0] = N_axial;
	stress[1] = My1;
	stress[2] = Mz1;
	stress[3] = T_torque;
	stress[4] = axial_stress;
	stress[5] = My2;
	stress[6] = Mz2;
	stress[7] = T_torque;
}
