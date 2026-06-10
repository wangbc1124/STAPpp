/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Shell4.h"
#include "Material.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CShell4::CShell4()
{
	NEN_ = 4;
	nodes_ = new CNode*[NEN_];

	ND_ = 24;   // 6 DOF/node x 4 nodes
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
CShell4::~CShell4()
{
}

//	Read element data from stream Input
bool CShell4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2, N3, N4;

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void CShell4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

//	Generate location matrix (all 6 DOFs per node)
void CShell4::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		for (unsigned int d = 0; d < 6; d++)
			LocationMatrix_[i++] = nodes_[N]->bcode[d];
	}
}

// --- Helper: Q4 shape functions ---
static void ShellShape(double xi, double eta, double N[4], double dNdxi[4], double dNdeta[4])
{
	N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
	N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
	N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
	N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);

	dNdxi[0] = -0.25 * (1.0 - eta);
	dNdxi[1] =  0.25 * (1.0 - eta);
	dNdxi[2] =  0.25 * (1.0 + eta);
	dNdxi[3] = -0.25 * (1.0 + eta);

	dNdeta[0] = -0.25 * (1.0 - xi);
	dNdeta[1] = -0.25 * (1.0 + xi);
	dNdeta[2] =  0.25 * (1.0 + xi);
	dNdeta[3] =  0.25 * (1.0 - xi);
}

//	Calculate element stiffness matrix (24x24)
void CShell4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
	double E  = mat->E;
	double nu = mat->Nu;
	double t  = mat->Thickness;

	double x[4], y[4], z[4];
	for (int i = 0; i < 4; i++) {
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
		z[i] = nodes_[i]->XYZ[2];
	}

	// --- 1. Compute local coordinate system ---
	// Normal: average of cross products of diagonals
	double v13[3] = {x[2]-x[0], y[2]-y[0], z[2]-z[0]};
	double v24[3] = {x[3]-x[1], y[3]-y[1], z[3]-z[1]};
	double ez[3] = {
		v13[1]*v24[2] - v13[2]*v24[1],
		v13[2]*v24[0] - v13[0]*v24[2],
		v13[0]*v24[1] - v13[1]*v24[0]
	};
	double nz = sqrt(ez[0]*ez[0] + ez[1]*ez[1] + ez[2]*ez[2]);
	if (nz < 1e-20) nz = 1.0;
	ez[0] /= nz; ez[1] /= nz; ez[2] /= nz;

	// Local x': direction of side 1->2
	double ex[3] = {x[1]-x[0], y[1]-y[0], z[1]-z[0]};
	double nx = sqrt(ex[0]*ex[0] + ex[1]*ex[1] + ex[2]*ex[2]);
	if (nx < 1e-20) {
		ex[0] = 1.0; ex[1] = 0.0; ex[2] = 0.0;
	} else {
		ex[0] /= nx; ex[1] /= nx; ex[2] /= nx;
	}

	// Local y': z' x x'
	double ey[3] = {
		ez[1]*ex[2] - ez[2]*ex[1],
		ez[2]*ex[0] - ez[0]*ex[2],
		ez[0]*ex[1] - ez[1]*ex[0]
	};
	double ny = sqrt(ey[0]*ey[0] + ey[1]*ey[1] + ey[2]*ey[2]);
	if (ny < 1e-20) {
		ey[0] = 0.0; ey[1] = 1.0; ey[2] = 0.0;
	} else {
		ey[0] /= ny; ey[1] /= ny; ey[2] /= ny;
	}

	// Recompute ex = ey x ez (ensure orthogonality)
	ex[0] = ey[1]*ez[2] - ey[2]*ez[1];
	ex[1] = ey[2]*ez[0] - ey[0]*ez[2];
	ex[2] = ey[0]*ez[1] - ey[1]*ez[0];

	// --- 2. Build rotation matrix R(3x3) ---
	// R * v_local = v_global
	double R[3][3] = {
		{ex[0], ey[0], ez[0]},
		{ex[1], ey[1], ez[1]},
		{ex[2], ey[2], ez[2]}
	};

	// --- 3. Project node coordinates onto local xy plane ---
	double xloc[4], yloc[4];
	for (int i = 0; i < 4; i++) {
		double dx = x[i] - x[0];  // relative to first node
		double dy = y[i] - y[0];
		double dz = z[i] - z[0];
		xloc[i] = ex[0]*dx + ex[1]*dy + ex[2]*dz;
		yloc[i] = ey[0]*dx + ey[1]*dy + ey[2]*dz;
	}

	// --- 4. Assemble local stiffness matrix K_local(24x24) ---
	double K[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++)
			K[i][j] = 0.0;

	// Constitutive matrices
	double Dm_fac = E * t / (1.0 - nu*nu);
	double Dm[3][3] = {
		{Dm_fac, Dm_fac*nu, 0.0},
		{Dm_fac*nu, Dm_fac, 0.0},
		{0.0, 0.0, Dm_fac*(1.0-nu)/2.0}
	};

	double Db_fac = E * t*t*t / (12.0 * (1.0 - nu*nu));
	double Db[3][3] = {
		{Db_fac, Db_fac*nu, 0.0},
		{Db_fac*nu, Db_fac, 0.0},
		{0.0, 0.0, Db_fac*(1.0-nu)/2.0}
	};

	double G = E / (2.0 * (1.0 + nu));
	double k_s = 5.0 / 6.0;
	double Ds = k_s * G * t;

	// --- 4a. Membrane stiffness (2x2 Gauss): DOF u,v per node ---
	// Indices in K_local[24]: node N -> u=6N, v=6N+1
	double gp = 1.0 / sqrt(3.0);
	double gauss_pts[2] = {-gp, gp};

	// Check element quality at center first
	bool valid_element = true;
	{
		double Nc[4], dNc_dxi[4], dNc_deta[4];
		ShellShape(0.0, 0.0, Nc, dNc_dxi, dNc_deta);
		double Jc11 = 0, Jc12 = 0, Jc21 = 0, Jc22 = 0;
		for (int n = 0; n < 4; n++) {
			Jc11 += dNc_dxi[n] * xloc[n];   Jc12 += dNc_dxi[n] * yloc[n];
			Jc21 += dNc_deta[n] * xloc[n];  Jc22 += dNc_deta[n] * yloc[n];
		}
		double detJc = Jc11*Jc22 - Jc12*Jc21;
		if (detJc < 1e-8) {
			valid_element = false;
		}
	}

	if (!valid_element) {
		// Skip invalid element �� set stiffness to zero (will be handled by neighbors)
		for (int j = 0; j < 24; j++)
			for (int i = 0; i <= j; i++)
				Matrix[(j + 1) * j / 2 + j - i] = 0.0;
		return;
	}

	for (int gi = 0; gi < 2; gi++) {
		for (int gj = 0; gj < 2; gj++) {
			double xi_p = gauss_pts[gi];
			double eta_p = gauss_pts[gj];

			double N[4], dNdxi[4], dNdeta[4];
			ShellShape(xi_p, eta_p, N, dNdxi, dNdeta);

			// Jacobian (2x2)
			double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
			for (int n = 0; n < 4; n++) {
				J11 += dNdxi[n] * xloc[n];   J12 += dNdxi[n] * yloc[n];
				J21 += dNdeta[n] * xloc[n];  J22 += dNdeta[n] * yloc[n];
			}
			double detJ = J11*J22 - J12*J21;
			if (detJ <= 0) detJ = 1e-20;
			double inv_det = 1.0 / detJ;

			// Physical derivatives
			double dNdx[4], dNdy[4];
			for (int n = 0; n < 4; n++) {
				dNdx[n] = ( J22*dNdxi[n] - J12*dNdeta[n]) * inv_det;
				dNdy[n] = (-J21*dNdxi[n] + J11*dNdeta[n]) * inv_det;
			}

			// Membrane B-matrix (3x24): strains [eps_x, eps_y, gamma_xy]
			// Node n: eps_x from dNdx * u_n, eps_y from dNdy * v_n, gamma_xy from dNdy*u_n + dNdx*v_n
			double Bm[3][24];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 24; j++)
					Bm[i][j] = 0.0;

			for (int n = 0; n < 4; n++) {
				int c = 6 * n;
				Bm[0][c+0] = dNdx[n];            // eps_x = sum dNdx * u
				Bm[1][c+1] = dNdy[n];            // eps_y = sum dNdy * v
				Bm[2][c+0] = dNdy[n];            // gamma_xy = sum (dNdy*u + dNdx*v)
				Bm[2][c+1] = dNdx[n];
			}

			// Km += Bm^T * Dm * Bm * detJ (weight=1 for 2x2)
			// First: Dm_Bm = Dm * Bm
			double DmBm[3][24];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 24; j++) {
					DmBm[i][j] = 0.0;
					for (int k = 0; k < 3; k++)
						DmBm[i][j] += Dm[i][k] * Bm[k][j];
				}

			for (int i = 0; i < 24; i++)
				for (int j = i; j < 24; j++) {
					double val = 0.0;
					for (int k = 0; k < 3; k++)
						val += Bm[k][i] * DmBm[k][j];
					K[i][j] += val * detJ;
				}
		}
	}

	// --- 4b. Plate bending stiffness (2x2 Gauss): DOF w, θx, θy per node ---
	// Plate DOF per node: θx=6N+3, θy=6N+4, w=6N+2
	for (int gi = 0; gi < 2; gi++) {
		for (int gj = 0; gj < 2; gj++) {
			double xi_p = gauss_pts[gi];
			double eta_p = gauss_pts[gj];

			double N[4], dNdxi[4], dNdeta[4];
			ShellShape(xi_p, eta_p, N, dNdxi, dNdeta);

			// Jacobian (same as membrane)
			double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
			for (int n = 0; n < 4; n++) {
				J11 += dNdxi[n] * xloc[n];   J12 += dNdxi[n] * yloc[n];
				J21 += dNdeta[n] * xloc[n];  J22 += dNdeta[n] * yloc[n];
			}
			double detJ = J11*J22 - J12*J21;
			if (detJ <= 0) detJ = 1e-20;
			double inv_det = 1.0 / detJ;

			double dNdx[4], dNdy[4];
			for (int n = 0; n < 4; n++) {
				dNdx[n] = ( J22*dNdxi[n] - J12*dNdeta[n]) * inv_det;
				dNdy[n] = (-J21*dNdxi[n] + J11*dNdeta[n]) * inv_det;
			}

			// Bending B-matrix Bb(3x24): curvatures [κ_x, κ_y, κ_xy]
			double Bb[3][24];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 24; j++)
					Bb[i][j] = 0.0;

			for (int n = 0; n < 4; n++) {
				int c = 6 * n;
				Bb[0][c+4] =  dNdx[n];          // κ_x = sum dNdx * θy
				Bb[1][c+3] = -dNdy[n];          // κ_y = sum -dNdy * θx
				Bb[2][c+3] = -dNdx[n];          // κ_xy = sum (-dNdx*θx + dNdy*θy)
				Bb[2][c+4] =  dNdy[n];
			}

			double DbBb[3][24];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 24; j++) {
					DbBb[i][j] = 0.0;
					for (int k = 0; k < 3; k++)
						DbBb[i][j] += Db[i][k] * Bb[k][j];
				}

			for (int i = 0; i < 24; i++)
				for (int j = i; j < 24; j++) {
					double val = 0.0;
					for (int k = 0; k < 3; k++)
						val += Bb[k][i] * DbBb[k][j];
					K[i][j] += val * detJ;
				}
		}
	}

	// --- 4c. Shear stiffness (1-point reduced integration) ---
	{
		double N[4], dNdxi[4], dNdeta[4];
		ShellShape(0.0, 0.0, N, dNdxi, dNdeta);

		double J11 = 0, J12 = 0, J21 = 0, J22 = 0;
		for (int n = 0; n < 4; n++) {
			J11 += dNdxi[n] * xloc[n];   J12 += dNdxi[n] * yloc[n];
			J21 += dNdeta[n] * xloc[n];  J22 += dNdeta[n] * yloc[n];
		}
		double detJ = J11*J22 - J12*J21;
		if (detJ <= 0) detJ = 1e-20;
		double inv_det = 1.0 / detJ;

		double dNdx[4], dNdy[4];
		for (int n = 0; n < 4; n++) {
			dNdx[n] = ( J22*dNdxi[n] - J12*dNdeta[n]) * inv_det;
			dNdy[n] = (-J21*dNdxi[n] + J11*dNdeta[n]) * inv_det;
		}

		// Shear B-matrix Bs(2x24): shear strains [γ_xz, γ_yz]
		double Bs[2][24];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 24; j++)
				Bs[i][j] = 0.0;

		for (int n = 0; n < 4; n++) {
			int c = 6 * n;
			Bs[0][c+4] =  N[n];                // γ_xz = sum (N*θy + dNdx*w)
			Bs[0][c+2] =  dNdx[n];
			Bs[1][c+3] = -N[n];                // γ_yz = sum (-N*θx + dNdy*w)
			Bs[1][c+2] =  dNdy[n];
		}

		for (int i = 0; i < 24; i++)
			for (int j = i; j < 24; j++) {
				double val = (Bs[0][i]*Bs[0][j] + Bs[1][i]*Bs[1][j]) * Ds;
				K[i][j] += val * detJ * 4.0;  // weight=4 for 1-pt quadrature over 2x2 area
			}
	}

	// --- 4d. Drilling DOF stabilization ---
	double area = 0.0;
	{
		double vx1 = xloc[1]-xloc[0], vy1 = yloc[1]-yloc[0];
		double vx2 = xloc[2]-xloc[1], vy2 = yloc[2]-yloc[1];
		double vx3 = xloc[3]-xloc[2], vy3 = yloc[3]-yloc[2];
		double vx4 = xloc[0]-xloc[3], vy4 = yloc[0]-yloc[3];
		area = 0.5 * fabs(vx1*vy2 - vy1*vx2 + vx2*vy3 - vy2*vx3 + vx3*vy4 - vy3*vx4 + vx4*vy1 - vy4*vx1);
	}
	// Drilling stiffness proportional to max bending stiffness diagonal
	// to properly constrain θz without dominating the solution
	double Db_diag = E * t*t*t / (12.0 * (1.0 - nu*nu));
	double drill_K = Db_diag * area * 0.1;  // ~10% of bending stiffness scale
	for (int n = 0; n < 4; n++) {
		int c = 6 * n + 5;  // θz DOF
		K[c][c] += drill_K;
	}

	// --- 5. Build rotation matrix T(24x24) ---
	double T[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++)
			T[i][j] = 0.0;

	// T maps global → local: u_local = T * u_global
	// R maps local → global, so T uses R^T = R[j][i]
	for (int n = 0; n < 4; n++) {
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

	// --- 6. K_global = T^T * K_local * T ---
	// Step 1: TK = T^T * K
	double TK[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++) {
			TK[i][j] = 0.0;
			for (int m = 0; m < 24; m++)
				TK[i][j] += T[m][i] * K[m][j];
		}

	// Step 2: K_global = TK * T
	double Kg[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++) {
			Kg[i][j] = 0.0;
			for (int m = 0; m < 24; m++)
				Kg[i][j] += TK[i][m] * T[m][j];
		}

	
// --- 7. Copy upper-triangular to column-by-column skyline storage ---
	for (int j = 0; j < 24; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = Kg[i][j];
}

//	Calculate element stress
void CShell4::ElementStress(double* stress, double* Displacement)
{
	CQ4Material* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
	double E  = mat->E;
	double nu = mat->Nu;
	double t  = mat->Thickness;

	double x[4], y[4], z[4];
	for (int i = 0; i < 4; i++) {
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
		z[i] = nodes_[i]->XYZ[2];
	}

	// Recompute local axes (same as stiffness)
	double v13[3] = {x[2]-x[0], y[2]-y[0], z[2]-z[0]};
	double v24[3] = {x[3]-x[1], y[3]-y[1], z[3]-z[1]};
	double ez[3] = {v13[1]*v24[2]-v13[2]*v24[1], v13[2]*v24[0]-v13[0]*v24[2], v13[0]*v24[1]-v13[1]*v24[0]};
	double nz = sqrt(ez[0]*ez[0]+ez[1]*ez[1]+ez[2]*ez[2]);
	if (nz < 1e-20) nz = 1.0;
	ez[0] /= nz; ez[1] /= nz; ez[2] /= nz;

	double ex[3] = {x[1]-x[0], y[1]-y[0], z[1]-z[0]};
	double nx = sqrt(ex[0]*ex[0]+ex[1]*ex[1]+ex[2]*ex[2]);
	if (nx < 1e-20) {ex[0]=1.0; ex[1]=0.0; ex[2]=0.0;}
	else {ex[0]/=nx; ex[1]/=nx; ex[2]/=nx;}

	double ey[3] = {ez[1]*ex[2]-ez[2]*ex[1], ez[2]*ex[0]-ez[0]*ex[2], ez[0]*ex[1]-ez[1]*ex[0]};
	double ny = sqrt(ey[0]*ey[0]+ey[1]*ey[1]+ey[2]*ey[2]);
	if (ny < 1e-20) {ey[0]=0.0; ey[1]=1.0; ey[2]=0.0;}
	else {ey[0]/=ny; ey[1]/=ny; ey[2]/=ny;}

	ex[0] = ey[1]*ez[2]-ey[2]*ez[1];
	ex[1] = ey[2]*ez[0]-ey[0]*ez[2];
	ex[2] = ey[0]*ez[1]-ey[1]*ez[0];

	double R[3][3] = {{ex[0],ey[0],ez[0]},{ex[1],ey[1],ez[1]},{ex[2],ey[2],ez[2]}};

	// Project to local coordinates
	double xloc[4], yloc[4];
	for (int i = 0; i < 4; i++) {
		double dx = x[i]-x[0], dy = y[i]-y[0], dz = z[i]-z[0];
		xloc[i] = ex[0]*dx+ex[1]*dy+ex[2]*dz;
		yloc[i] = ey[0]*dx+ey[1]*dy+ey[2]*dz;
	}

	// Extract global displacements, transform to local
	double u_glob[24];
	for (int i = 0; i < 24; i++) {
		if (LocationMatrix_[i]) u_glob[i] = Displacement[LocationMatrix_[i]-1];
		else u_glob[i] = 0.0;
	}

	// T^T * u_global = u_local (R^T is R transpose, since R is orthonormal)
	double u_loc[24];
	for (int i = 0; i < 24; i++) u_loc[i] = 0.0;
	for (int n = 0; n < 4; n++) {
		int b = n*6;
		for (int d = 0; d < 3; d++) {
			double val = 0.0;
			for (int m = 0; m < 3; m++) val += R[m][d] * u_glob[b+m];
			u_loc[b+d] = val;
		}
		for (int d = 0; d < 3; d++) {
			double val = 0.0;
			for (int m = 0; m < 3; m++) val += R[m][d] * u_glob[b+3+m];
			u_loc[b+3+d] = val;
		}
	}

	// Compute at element center (xi=eta=0)
	double N[4], dNdxi[4], dNdeta[4];
	ShellShape(0.0, 0.0, N, dNdxi, dNdeta);

	double J11=0, J12=0, J21=0, J22=0;
	for (int n=0; n<4; n++) {
		J11+=dNdxi[n]*xloc[n]; J12+=dNdxi[n]*yloc[n];
		J21+=dNdeta[n]*xloc[n]; J22+=dNdeta[n]*yloc[n];
	}
	double detJ = J11*J22-J12*J21;
	if (detJ <= 0) detJ = 1e-20;
	double inv_det=1.0/detJ;

	double dNdx[4], dNdy[4];
	for (int n=0; n<4; n++) {
		dNdx[n]=(J22*dNdxi[n]-J12*dNdeta[n])*inv_det;
		dNdy[n]=(-J21*dNdxi[n]+J11*dNdeta[n])*inv_det;
	}

	// Membrane strains at center
	double eps_x=0, eps_y=0, gamma_xy=0;
	for (int n=0; n<4; n++) {
		eps_x    += dNdx[n] * u_loc[6*n+0];
		eps_y    += dNdy[n] * u_loc[6*n+1];
		gamma_xy += dNdy[n]*u_loc[6*n+0] + dNdx[n]*u_loc[6*n+1];
	}
	double Dm = E*t/(1.0-nu*nu);
	double sx = Dm*(eps_x+nu*eps_y);
	double sy = Dm*(nu*eps_x+eps_y);
	double txy = Dm*(1.0-nu)/2.0*gamma_xy;

	// Bending curvatures at center
	double kap_x=0, kap_y=0, kap_xy=0;
	for (int n=0; n<4; n++) {
		kap_x  += dNdx[n] * u_loc[6*n+4];
		kap_y  -= dNdy[n] * u_loc[6*n+3];
		kap_xy -= dNdx[n]*u_loc[6*n+3] - dNdy[n]*u_loc[6*n+4];
	}
	double Dbf = E*t*t*t/(12.0*(1.0-nu*nu));
	double Mx  = Dbf*(kap_x+nu*kap_y);
	double My  = Dbf*(nu*kap_x+kap_y);
	double Mxy = Dbf*(1.0-nu)/2.0*kap_xy;

	// Shear forces at center
	double gamma_xz=0, gamma_yz=0;
	for (int n=0; n<4; n++) {
		gamma_xz += N[n]*u_loc[6*n+4] + dNdx[n]*u_loc[6*n+2];
		gamma_yz -= N[n]*u_loc[6*n+3] - dNdy[n]*u_loc[6*n+2];
	}
	double G = E/(2.0*(1.0+nu));
	double Ds_val = 5.0/6.0*G*t;
	double Qx = Ds_val*gamma_xz;
	double Qy = Ds_val*gamma_yz;

	stress[0] = sx;
	stress[1] = sy;
	stress[2] = txy;
	stress[3] = Mx;
	stress[4] = My;
	stress[5] = Mxy;
	stress[6] = Qx;
	stress[7] = Qy;
}
