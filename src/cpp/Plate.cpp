#include "Plate.h"

#include <cmath>
#include <iomanip>

using namespace std;

CPlate::CPlate()
{
	NEN_ = 4;
	nodes_ = new CNode*[NEN_];

	ND_ = 12;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

CPlate::~CPlate()
{
}

bool CPlate::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
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

void CPlate::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CPlate::GenerateLocationMatrix()
{
	// DOF order per node: theta_x, theta_y, w (mapped to bcode[0], [1], [2])
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		LocationMatrix_[i++] = nodes_[N]->bcode[0];
		LocationMatrix_[i++] = nodes_[N]->bcode[1];
		LocationMatrix_[i++] = nodes_[N]->bcode[2];
	}
}

// Shape functions for Q4
static void Q4Shape(double xi, double eta, double N[4], double dNdxi[4], double dNdeta[4])
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

void CPlate::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
	double E = material->E;
	double nu = material->Nu;
	double t = material->Thickness;

	// Node coordinates
	double x[4], y[4];
	for (int i = 0; i < 4; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
	}

	// Bending constitutive matrix Db (3x3)
	double Db_fac = E * t * t * t / (12.0 * (1.0 - nu * nu));
	double Db[3][3];
	Db[0][0] = Db_fac;
	Db[0][1] = Db_fac * nu;
	Db[0][2] = 0.0;
	Db[1][0] = Db_fac * nu;
	Db[1][1] = Db_fac;
	Db[1][2] = 0.0;
	Db[2][0] = 0.0;
	Db[2][1] = 0.0;
	Db[2][2] = Db_fac * (1.0 - nu) / 2.0;

	// Shear constitutive matrix Ds (2x2)
	double G = E / (2.0 * (1.0 + nu));
	double k = 5.0 / 6.0;
	double Ds = k * G * t;

	// Jacobian and stiffness
	double K[12][12];
	for (int i = 0; i < 12; i++)
		for (int j = 0; j < 12; j++)
			K[i][j] = 0.0;

	// === Bending stiffness: 2x2 Gauss quadrature ===
	double gp = 1.0 / sqrt(3.0);
	double gauss_pts[2] = {-gp, gp};

	for (int gi = 0; gi < 2; gi++)
	{
		for (int gj = 0; gj < 2; gj++)
		{
			double xi_p = gauss_pts[gi];
			double eta_p = gauss_pts[gj];

			double N[4], dNdxi[4], dNdeta[4];
			Q4Shape(xi_p, eta_p, N, dNdxi, dNdeta);

			// Jacobian
			double J[2][2];
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					J[i][j] = 0.0;

			for (int n = 0; n < 4; n++)
			{
				J[0][0] += dNdxi[n] * x[n];
				J[0][1] += dNdxi[n] * y[n];
				J[1][0] += dNdeta[n] * x[n];
				J[1][1] += dNdeta[n] * y[n];
			}

			double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
			double invJ[2][2];
			double inv_det = 1.0 / detJ;
			invJ[0][0] =  J[1][1] * inv_det;
			invJ[0][1] = -J[0][1] * inv_det;
			invJ[1][0] = -J[1][0] * inv_det;
			invJ[1][1] =  J[0][0] * inv_det;

			// Physical derivatives
			double dNdx[4], dNdy[4];
			for (int n = 0; n < 4; n++)
			{
				dNdx[n] = invJ[0][0] * dNdxi[n] + invJ[0][1] * dNdeta[n];
				dNdy[n] = invJ[1][0] * dNdxi[n] + invJ[1][1] * dNdeta[n];
			}

			// Bending B-matrix (3x12)
			double Bb[3][12];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 12; j++)
					Bb[i][j] = 0.0;

			for (int n = 0; n < 4; n++)
			{
				int col = 3 * n;
				Bb[0][col+1] =  dNdx[n];
				Bb[1][col]   = -dNdy[n];
				Bb[2][col]   = -dNdx[n];
				Bb[2][col+1] =  dNdy[n];
			}

			// Kb += Bb^T * Db * Bb * detJ (weight = 1 for 2x2 Gauss)
			double DbBb[3][12];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 12; j++)
				{
					DbBb[i][j] = 0.0;
					for (int m = 0; m < 3; m++)
						DbBb[i][j] += Db[i][m] * Bb[m][j];
				}

			for (int i = 0; i < 12; i++)
				for (int j = i; j < 12; j++)
				{
					double val = 0.0;
					for (int m = 0; m < 3; m++)
						val += Bb[m][i] * DbBb[m][j];
					K[i][j] += val * detJ;
				}
		}
	}

	// === Shear stiffness: 1-point Gauss (reduced integration) ===
	{
		double xi_p = 0.0, eta_p = 0.0;

		double N[4], dNdxi[4], dNdeta[4];
		Q4Shape(xi_p, eta_p, N, dNdxi, dNdeta);

		double J[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				J[i][j] = 0.0;

		for (int n = 0; n < 4; n++)
		{
			J[0][0] += dNdxi[n] * x[n];
			J[0][1] += dNdxi[n] * y[n];
			J[1][0] += dNdeta[n] * x[n];
			J[1][1] += dNdeta[n] * y[n];
		}

		double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
		double invJ[2][2];
		double inv_det = 1.0 / detJ;
		invJ[0][0] =  J[1][1] * inv_det;
		invJ[0][1] = -J[0][1] * inv_det;
		invJ[1][0] = -J[1][0] * inv_det;
		invJ[1][1] =  J[0][0] * inv_det;

		double dNdx[4], dNdy[4];
		for (int n = 0; n < 4; n++)
		{
			dNdx[n] = invJ[0][0] * dNdxi[n] + invJ[0][1] * dNdeta[n];
			dNdy[n] = invJ[1][0] * dNdxi[n] + invJ[1][1] * dNdeta[n];
		}

		// Shear B-matrix (2x12)
		double Bs[2][12];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 12; j++)
				Bs[i][j] = 0.0;

		for (int n = 0; n < 4; n++)
		{
			int col = 3 * n;
			Bs[0][col+1] =  N[n];
			Bs[0][col+2] =  dNdx[n];
			Bs[1][col]   = -N[n];
			Bs[1][col+2] =  dNdy[n];
		}

		// Ks += Bs^T * Ds * Bs * detJ * 4 (weight = 4 for center-point quadrature)
		double DsBs[2][12];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 12; j++)
				DsBs[i][j] = Ds * Bs[i][j];

		for (int i = 0; i < 12; i++)
			for (int j = i; j < 12; j++)
			{
				double val = 0.0;
				for (int m = 0; m < 2; m++)
					val += Bs[m][i] * DsBs[m][j];
				K[i][j] += val * detJ * 4.0; // weight = 4 for 1-point quadrature over 2x2 domain
			}
	}

	// Copy upper-triangular K to column-by-column storage
	for (int j = 0; j < 12; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CPlate::ElementStress(double* stress, double* Displacement)
{
	CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
	double E = material->E;
	double nu = material->Nu;
	double t = material->Thickness;

	double x[4], y[4];
	for (int i = 0; i < 4; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
	}

	// Extract displacements
	double u[12];
	for (int i = 0; i < 12; i++)
	{
		if (LocationMatrix_[i])
			u[i] = Displacement[LocationMatrix_[i] - 1];
		else
			u[i] = 0.0;
	}

	// Compute at element center (xi=eta=0)
	double N[4], dNdxi[4], dNdeta[4];
	Q4Shape(0.0, 0.0, N, dNdxi, dNdeta);

	double J[2][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			J[i][j] = 0.0;

	for (int n = 0; n < 4; n++)
	{
		J[0][0] += dNdxi[n] * x[n];
		J[0][1] += dNdxi[n] * y[n];
		J[1][0] += dNdeta[n] * x[n];
		J[1][1] += dNdeta[n] * y[n];
	}

	double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	double invJ[2][2];
	double inv_det = 1.0 / detJ;
	invJ[0][0] =  J[1][1] * inv_det;
	invJ[0][1] = -J[0][1] * inv_det;
	invJ[1][0] = -J[1][0] * inv_det;
	invJ[1][1] =  J[0][0] * inv_det;

	double dNdx[4], dNdy[4];
	for (int n = 0; n < 4; n++)
	{
		dNdx[n] = invJ[0][0] * dNdxi[n] + invJ[0][1] * dNdeta[n];
		dNdy[n] = invJ[1][0] * dNdxi[n] + invJ[1][1] * dNdeta[n];
	}

	// Bending B-matrix (3x12)
	double Bb[3][12];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 12; j++)
			Bb[i][j] = 0.0;

	for (int n = 0; n < 4; n++)
	{
		int col = 3 * n;
		Bb[0][col+1] =  dNdx[n];
		Bb[1][col]   = -dNdy[n];
		Bb[2][col]   = -dNdx[n];
		Bb[2][col+1] =  dNdy[n];
	}

	// Curvatures: kappa = Bb * u
	double kappa[3];
	for (int i = 0; i < 3; i++)
	{
		kappa[i] = 0.0;
		for (int j = 0; j < 12; j++)
			kappa[i] += Bb[i][j] * u[j];
	}

	// Bending moments: M = Db * kappa
	double Db_fac = E * t * t * t / (12.0 * (1.0 - nu * nu));
	double Mx = Db_fac * (kappa[0] + nu * kappa[1]);
	double My = Db_fac * (nu * kappa[0] + kappa[1]);
	double Mxy = Db_fac * (1.0 - nu) / 2.0 * kappa[2];

	// Shear B-matrix (2x12)
	double Bs[2][12];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 12; j++)
			Bs[i][j] = 0.0;

	for (int n = 0; n < 4; n++)
	{
		int col = 3 * n;
		Bs[0][col+1] =  N[n];
		Bs[0][col+2] =  dNdx[n];
		Bs[1][col]   = -N[n];
		Bs[1][col+2] =  dNdy[n];
	}

	// Shear strains: gamma = Bs * u
	double gamma_xz = 0.0, gamma_yz = 0.0;
	for (int j = 0; j < 12; j++)
	{
		gamma_xz += Bs[0][j] * u[j];
		gamma_yz += Bs[1][j] * u[j];
	}

	// Shear forces: Q = Ds * gamma
	double G = E / (2.0 * (1.0 + nu));
	double k = 5.0 / 6.0;
	double Ds = k * G * t;
	double Qx = Ds * gamma_xz;
	double Qy = Ds * gamma_yz;

	stress[0] = Mx;
	stress[1] = My;
	stress[2] = Mxy;
	stress[3] = Qx;
	stress[4] = Qy;
}
