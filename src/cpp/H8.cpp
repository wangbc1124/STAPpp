#include "H8.h"

#include <cmath>
#include <iomanip>

using namespace std;

CH8::CH8()
{
	NEN_ = 8;
	nodes_ = new CNode*[NEN_];

	ND_ = 24;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

CH8::~CH8()
{
}

bool CH8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N[8];

	Input >> N[0] >> N[1] >> N[2] >> N[3] >> N[4] >> N[5] >> N[6] >> N[7] >> MSet;
	ElementMaterial_ = dynamic_cast<CH8Material*>(MaterialSets) + MSet - 1;

	for (unsigned int i = 0; i < 8; i++)
		nodes_[i] = &NodeList[N[i] - 1];

	return true;
}

void CH8::Write(COutputter& output)
{
	for (unsigned int i = 0; i < 8; i++)
		output << setw(9) << nodes_[i]->NodeNumber;
	output << setw(12) << ElementMaterial_->nset << endl;
}

void CH8::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		LocationMatrix_[i++] = nodes_[N]->bcode[0];
		LocationMatrix_[i++] = nodes_[N]->bcode[1];
		LocationMatrix_[i++] = nodes_[N]->bcode[2];
	}
}

// Natural coordinates of 8 nodes
static const double xi[8]   = {-1,  1,  1, -1, -1,  1,  1, -1};
static const double eta[8]  = {-1, -1,  1,  1, -1, -1,  1,  1};
static const double zeta[8] = {-1, -1, -1, -1,  1,  1,  1,  1};

void CH8::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CH8Material* material = dynamic_cast<CH8Material*>(ElementMaterial_);
	double E = material->E;
	double nu = material->Nu;

	// 3D constitutive matrix D (6x6)
	double D[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			D[i][j] = 0.0;

	double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
	D[0][0] = factor * (1.0 - nu);  D[0][1] = factor * nu;        D[0][2] = factor * nu;
	D[1][0] = factor * nu;          D[1][1] = factor * (1.0 - nu);  D[1][2] = factor * nu;
	D[2][0] = factor * nu;          D[2][1] = factor * nu;          D[2][2] = factor * (1.0 - nu);
	D[3][3] = factor * (1.0 - 2.0 * nu) / 2.0;
	D[4][4] = factor * (1.0 - 2.0 * nu) / 2.0;
	D[5][5] = factor * (1.0 - 2.0 * nu) / 2.0;

	// Node coordinates
	double x[8], y[8], z[8];
	for (int i = 0; i < 8; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
		z[i] = nodes_[i]->XYZ[2];
	}

	// 2x2x2 Gauss quadrature
	double gp = 1.0 / sqrt(3.0);
	double gauss_pts[2] = {-gp, gp};

	double K[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++)
			K[i][j] = 0.0;

	for (int gi = 0; gi < 2; gi++)
	{
		for (int gj = 0; gj < 2; gj++)
		{
			for (int gk = 0; gk < 2; gk++)
			{
				double xi_p = gauss_pts[gi];
				double eta_p = gauss_pts[gj];
				double zeta_p = gauss_pts[gk];

				// Shape function derivatives in natural coords
				double dNdxi[8], dNdeta[8], dNdzeta[8];
				for (int n = 0; n < 8; n++)
				{
					dNdxi[n]   = 0.125 * xi[n]   * (1.0 + eta_p * eta[n]) * (1.0 + zeta_p * zeta[n]);
					dNdeta[n]  = 0.125 * eta[n]  * (1.0 + xi_p * xi[n])   * (1.0 + zeta_p * zeta[n]);
					dNdzeta[n] = 0.125 * zeta[n] * (1.0 + xi_p * xi[n])   * (1.0 + eta_p * eta[n]);
				}

				// Jacobian matrix J (3x3)
				double J[3][3];
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						J[i][j] = 0.0;

				for (int n = 0; n < 8; n++)
				{
					J[0][0] += dNdxi[n] * x[n];
					J[0][1] += dNdxi[n] * y[n];
					J[0][2] += dNdxi[n] * z[n];
					J[1][0] += dNdeta[n] * x[n];
					J[1][1] += dNdeta[n] * y[n];
					J[1][2] += dNdeta[n] * z[n];
					J[2][0] += dNdzeta[n] * x[n];
					J[2][1] += dNdzeta[n] * y[n];
					J[2][2] += dNdzeta[n] * z[n];
				}

				// Determinant of J
				double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
				            - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
				            + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

				// Inverse of J (3x3)
				double invJ[3][3];
				double inv_det = 1.0 / detJ;
				invJ[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv_det;
				invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * inv_det;
				invJ[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv_det;
				invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * inv_det;
				invJ[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv_det;
				invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * inv_det;
				invJ[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv_det;
				invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * inv_det;
				invJ[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv_det;

				// Physical derivatives dN/dx, dN/dy, dN/dz
				double dNdx[8], dNdy[8], dNdz[8];
				for (int n = 0; n < 8; n++)
				{
					dNdx[n] = invJ[0][0] * dNdxi[n] + invJ[0][1] * dNdeta[n] + invJ[0][2] * dNdzeta[n];
					dNdy[n] = invJ[1][0] * dNdxi[n] + invJ[1][1] * dNdeta[n] + invJ[1][2] * dNdzeta[n];
					dNdz[n] = invJ[2][0] * dNdxi[n] + invJ[2][1] * dNdeta[n] + invJ[2][2] * dNdzeta[n];
				}

				// B matrix (6x24)
				double B[6][24];
				for (int i = 0; i < 6; i++)
					for (int j = 0; j < 24; j++)
						B[i][j] = 0.0;

				for (int n = 0; n < 8; n++)
				{
					int col = 3 * n;
					B[0][col]   = dNdx[n];
					B[1][col+1] = dNdy[n];
					B[2][col+2] = dNdz[n];
					B[3][col]   = dNdy[n];
					B[3][col+1] = dNdx[n];
					B[4][col+1] = dNdz[n];
					B[4][col+2] = dNdy[n];
					B[5][col]   = dNdz[n];
					B[5][col+2] = dNdx[n];
				}

				// K += B^T * D * B * detJ (weight=1 at each Gauss point)
				double DB[6][24];
				for (int i = 0; i < 6; i++)
					for (int j = 0; j < 24; j++)
					{
						DB[i][j] = 0.0;
						for (int m = 0; m < 6; m++)
							DB[i][j] += D[i][m] * B[m][j];
					}

				for (int i = 0; i < 24; i++)
					for (int j = i; j < 24; j++)
					{
						double val = 0.0;
						for (int m = 0; m < 6; m++)
							val += B[m][i] * DB[m][j];
						K[i][j] += val * detJ;
					}
			}
		}
	}

	// Copy upper-triangular K to column-by-column storage
	int idx = 0;
	for (int j = 0; j < 24; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CH8::ElementStress(double* stress, double* Displacement)
{
	CH8Material* material = dynamic_cast<CH8Material*>(ElementMaterial_);
	double E = material->E;
	double nu = material->Nu;

	// Constitutive matrix
	double D[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			D[i][j] = 0.0;

	double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
	D[0][0] = factor * (1.0 - nu);  D[0][1] = factor * nu;        D[0][2] = factor * nu;
	D[1][0] = factor * nu;          D[1][1] = factor * (1.0 - nu);  D[1][2] = factor * nu;
	D[2][0] = factor * nu;          D[2][1] = factor * nu;          D[2][2] = factor * (1.0 - nu);
	D[3][3] = factor * (1.0 - 2.0 * nu) / 2.0;
	D[4][4] = factor * (1.0 - 2.0 * nu) / 2.0;
	D[5][5] = factor * (1.0 - 2.0 * nu) / 2.0;

	// Node coordinates
	double x[8], y[8], z[8];
	for (int i = 0; i < 8; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
		z[i] = nodes_[i]->XYZ[2];
	}

	// Extract displacements
	double u[24];
	for (int i = 0; i < 24; i++)
	{
		if (LocationMatrix_[i])
			u[i] = Displacement[LocationMatrix_[i] - 1];
		else
			u[i] = 0.0;
	}

	// Compute stress at element center (xi=eta=zeta=0)
	double xi_p = 0.0, eta_p = 0.0, zeta_p = 0.0;

	double dNdxi[8], dNdeta[8], dNdzeta[8];
	for (int n = 0; n < 8; n++)
	{
		dNdxi[n]   = 0.125 * xi[n]   * (1.0 + eta_p * eta[n]) * (1.0 + zeta_p * zeta[n]);
		dNdeta[n]  = 0.125 * eta[n]  * (1.0 + xi_p * xi[n])   * (1.0 + zeta_p * zeta[n]);
		dNdzeta[n] = 0.125 * zeta[n] * (1.0 + xi_p * xi[n])   * (1.0 + eta_p * eta[n]);
	}

	double J[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			J[i][j] = 0.0;

	for (int n = 0; n < 8; n++)
	{
		J[0][0] += dNdxi[n] * x[n];
		J[0][1] += dNdxi[n] * y[n];
		J[0][2] += dNdxi[n] * z[n];
		J[1][0] += dNdeta[n] * x[n];
		J[1][1] += dNdeta[n] * y[n];
		J[1][2] += dNdeta[n] * z[n];
		J[2][0] += dNdzeta[n] * x[n];
		J[2][1] += dNdzeta[n] * y[n];
		J[2][2] += dNdzeta[n] * z[n];
	}

	double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
	            - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
	            + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

	double inv_det = 1.0 / detJ;
	double invJ[3][3];
	invJ[0][0] =  (J[1][1] * J[2][2] - J[1][2] * J[2][1]) * inv_det;
	invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) * inv_det;
	invJ[0][2] =  (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv_det;
	invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) * inv_det;
	invJ[1][1] =  (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv_det;
	invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) * inv_det;
	invJ[2][0] =  (J[1][0] * J[2][1] - J[1][1] * J[2][0]) * inv_det;
	invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) * inv_det;
	invJ[2][2] =  (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv_det;

	double dNdx[8], dNdy[8], dNdz[8];
	for (int n = 0; n < 8; n++)
	{
		dNdx[n] = invJ[0][0] * dNdxi[n] + invJ[0][1] * dNdeta[n] + invJ[0][2] * dNdzeta[n];
		dNdy[n] = invJ[1][0] * dNdxi[n] + invJ[1][1] * dNdeta[n] + invJ[1][2] * dNdzeta[n];
		dNdz[n] = invJ[2][0] * dNdxi[n] + invJ[2][1] * dNdeta[n] + invJ[2][2] * dNdzeta[n];
	}

	// B matrix
	double B[6][24];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 24; j++)
			B[i][j] = 0.0;

	for (int n = 0; n < 8; n++)
	{
		int col = 3 * n;
		B[0][col]   = dNdx[n];
		B[1][col+1] = dNdy[n];
		B[2][col+2] = dNdz[n];
		B[3][col]   = dNdy[n];
		B[3][col+1] = dNdx[n];
		B[4][col+1] = dNdz[n];
		B[4][col+2] = dNdy[n];
		B[5][col]   = dNdz[n];
		B[5][col+2] = dNdx[n];
	}

	// Strain = B * u
	double strain[6];
	for (int i = 0; i < 6; i++)
	{
		strain[i] = 0.0;
		for (int j = 0; j < 24; j++)
			strain[i] += B[i][j] * u[j];
	}

	// Stress = D * strain
	for (int i = 0; i < 6; i++)
	{
		stress[i] = 0.0;
		for (int j = 0; j < 6; j++)
			stress[i] += D[i][j] * strain[j];
	}
}
