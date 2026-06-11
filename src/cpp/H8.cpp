#include "H8.h"

#include <cstdlib>
#include <cmath>
#include <iomanip>

using namespace std;

static double EnvDouble(const char* name, double default_value);
static bool EnvBool(const char* name, bool default_value);

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

double CH8::HourglassAlphaBase() const
{
	return EnvDouble("STAP_H8_ALPHA", 0.009);
}

double CH8::HourglassAlphaMin() const
{
	return EnvDouble("STAP_H8_ALPHA_MIN", 0.0055);
}

double CH8::GlobalYStiffnessScale() const
{
	return EnvDouble("STAP_H8_GLOBAL_Y_SCALE", 1.0);
}

bool CH8::UseSelectiveReducedIntegration() const
{
	return EnvBool("STAP_H8_SRI", false);
}

double CH8::SRIHourglassBlend() const
{
	return EnvDouble("STAP_H8_SRI_HG_BLEND", 0.0);
}

double CH8::FBHourglassScale() const
{
	return EnvDouble("STAP_H8_FB_HG_SCALE", 0.0);
}

bool CH8::UseOrthogonalHourglass() const
{
	return EnvBool("STAP_H8_ORTHO_HG", false);
}

double CH8::OrthogonalHourglassScale() const
{
	return EnvDouble("STAP_H8_ORTHO_HG_SCALE", HourglassAlphaBase());
}

double CH8::OrthogonalHourglassDirectionScale(int direction) const
{
	if (direction == 0)
		return EnvDouble("STAP_H8_ORTHO_HG_X_SCALE", 1.0);
	if (direction == 1)
		return EnvDouble("STAP_H8_ORTHO_HG_Y_SCALE", 1.0);
	return EnvDouble("STAP_H8_ORTHO_HG_Z_SCALE", 1.0);
}

double CH8RPier::HourglassAlphaBase() const
{
	return EnvDouble("STAP_PIER_H8_ALPHA", 0.00904);
}

double CH8RPier::HourglassAlphaMin() const
{
	return EnvDouble("STAP_PIER_H8_ALPHA_MIN", 0.0055);
}

double CH8RPier::GlobalYStiffnessScale() const
{
	return EnvDouble("STAP_PIER_GLOBAL_Y_SCALE", 1.0);
}

bool CH8RPier::UseSelectiveReducedIntegration() const
{
	return EnvBool("STAP_PIER_H8_SRI", false);
}

double CH8RPier::SRIHourglassBlend() const
{
	return EnvDouble("STAP_PIER_H8_SRI_HG_BLEND", 0.02);
}

double CH8RPier::FBHourglassScale() const
{
	return EnvDouble("STAP_PIER_H8_FB_HG_SCALE", 0.0);
}

bool CH8RPier::UseOrthogonalHourglass() const
{
	return EnvBool("STAP_PIER_H8_ORTHO_HG", true);
}

double CH8RPier::OrthogonalHourglassScale() const
{
	return EnvDouble("STAP_PIER_H8_ORTHO_HG_SCALE", 0.015);
}

double CH8RPier::OrthogonalHourglassDirectionScale(int direction) const
{
	if (direction == 0)
		return EnvDouble("STAP_PIER_H8_ORTHO_HG_X_SCALE", 1.0);
	if (direction == 1)
		return EnvDouble("STAP_PIER_H8_ORTHO_HG_Y_SCALE", 1.0);
	return EnvDouble("STAP_PIER_H8_ORTHO_HG_Z_SCALE", 0.25);
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

static double EnvDouble(const char* name, double default_value)
{
	const char* value = getenv(name);
	if (!value || !*value)
		return default_value;
	char* end = nullptr;
	double parsed = strtod(value, &end);
	return (end && end != value) ? parsed : default_value;
}

static bool EnvBool(const char* name, bool default_value)
{
	const char* value = getenv(name);
	if (!value || !*value)
		return default_value;
	if (value[0] == '0' || value[0] == 'n' || value[0] == 'N' || value[0] == 'f' || value[0] == 'F')
		return false;
	return true;
}

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

	double D_vol[6][6], D_dev[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
		{
			D_vol[i][j] = 0.0;
			D_dev[i][j] = D[i][j];
		}
	double bulk = E / (3.0 * (1.0 - 2.0 * nu));
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			D_vol[i][j] = bulk;
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			D_dev[i][j] = D[i][j] - D_vol[i][j];

	// Node coordinates
	double x[8], y[8], z[8];
	for (int i = 0; i < 8; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
		z[i] = nodes_[i]->XYZ[2];
	}

	// Single-point integration at center for reduced integration stiffness
	// B0 at center (xi=eta=zeta=0): shape function derivatives
	double K[24][24];
	double K_sri[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++)
		{
			K[i][j] = 0.0;
			K_sri[i][j] = 0.0;
		}

	double B0[6][24];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 24; j++)
			B0[i][j] = 0.0;

	double detJ0_center = 0.0;
	{
		double dNdxi0[8], dNdeta0[8], dNdzeta0[8];
		for (int n = 0; n < 8; n++)
		{
			dNdxi0[n]   = 0.125 * xi[n];
			dNdeta0[n]  = 0.125 * eta[n];
			dNdzeta0[n] = 0.125 * zeta[n];
		}

		double J0[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				J0[i][j] = 0.0;

		for (int n = 0; n < 8; n++)
		{
			J0[0][0] += dNdxi0[n] * x[n];   J0[0][1] += dNdxi0[n] * y[n];   J0[0][2] += dNdxi0[n] * z[n];
			J0[1][0] += dNdeta0[n] * x[n];  J0[1][1] += dNdeta0[n] * y[n];  J0[1][2] += dNdeta0[n] * z[n];
			J0[2][0] += dNdzeta0[n] * x[n]; J0[2][1] += dNdzeta0[n] * y[n]; J0[2][2] += dNdzeta0[n] * z[n];
		}

		double detJ0 = J0[0][0] * (J0[1][1] * J0[2][2] - J0[1][2] * J0[2][1])
		             - J0[0][1] * (J0[1][0] * J0[2][2] - J0[1][2] * J0[2][0])
		             + J0[0][2] * (J0[1][0] * J0[2][1] - J0[1][1] * J0[2][0]);
		detJ0_center = detJ0;

		// Only use 1-pt stiffness if Jacobian is well-conditioned
		if (detJ0 > 1e-15) {

		double inv_det0 = 1.0 / detJ0;
		double invJ0[3][3];
		invJ0[0][0] =  (J0[1][1] * J0[2][2] - J0[1][2] * J0[2][1]) * inv_det0;
		invJ0[0][1] = -(J0[0][1] * J0[2][2] - J0[0][2] * J0[2][1]) * inv_det0;
		invJ0[0][2] =  (J0[0][1] * J0[1][2] - J0[0][2] * J0[1][1]) * inv_det0;
		invJ0[1][0] = -(J0[1][0] * J0[2][2] - J0[1][2] * J0[2][0]) * inv_det0;
		invJ0[1][1] =  (J0[0][0] * J0[2][2] - J0[0][2] * J0[2][0]) * inv_det0;
		invJ0[1][2] = -(J0[0][0] * J0[1][2] - J0[0][2] * J0[1][0]) * inv_det0;
		invJ0[2][0] =  (J0[1][0] * J0[2][1] - J0[1][1] * J0[2][0]) * inv_det0;
		invJ0[2][1] = -(J0[0][0] * J0[2][1] - J0[0][1] * J0[2][0]) * inv_det0;
		invJ0[2][2] =  (J0[0][0] * J0[1][1] - J0[0][1] * J0[1][0]) * inv_det0;

		double dNdx0[8], dNdy0[8], dNdz0[8];
		for (int n = 0; n < 8; n++)
		{
			dNdx0[n] = invJ0[0][0] * dNdxi0[n] + invJ0[0][1] * dNdeta0[n] + invJ0[0][2] * dNdzeta0[n];
			dNdy0[n] = invJ0[1][0] * dNdxi0[n] + invJ0[1][1] * dNdeta0[n] + invJ0[1][2] * dNdzeta0[n];
			dNdz0[n] = invJ0[2][0] * dNdxi0[n] + invJ0[2][1] * dNdeta0[n] + invJ0[2][2] * dNdzeta0[n];
		}

		for (int n = 0; n < 8; n++)
		{
			int col = 3 * n;
			B0[0][col]   = dNdx0[n];
			B0[1][col+1] = dNdy0[n];
			B0[2][col+2] = dNdz0[n];
			B0[3][col]   = dNdy0[n]; B0[3][col+1] = dNdx0[n];
			B0[4][col+1] = dNdz0[n]; B0[4][col+2] = dNdy0[n];
			B0[5][col]   = dNdz0[n]; B0[5][col+2] = dNdx0[n];
		}

		// 1-pt stiffness: K1 = B0^T * D * B0 * 8*detJ0 (weight=8 for whole element)
		double DB0[6][24];
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 24; j++) {
				DB0[i][j] = 0.0;
				for (int m = 0; m < 6; m++)
					DB0[i][j] += D[i][m] * B0[m][j];
			}

		for (int i = 0; i < 24; i++)
			for (int j = i; j < 24; j++) {
				double val = 0.0;
				for (int m = 0; m < 6; m++)
					val += B0[m][i] * DB0[m][j];
				K[i][j] += val * detJ0 * 8.0;
			}

		double DB0Vol[6][24];
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 24; j++) {
				DB0Vol[i][j] = 0.0;
				for (int m = 0; m < 6; m++)
					DB0Vol[i][j] += D_vol[i][m] * B0[m][j];
			}
		for (int i = 0; i < 24; i++)
			for (int j = i; j < 24; j++) {
				double val = 0.0;
				for (int m = 0; m < 6; m++)
					val += B0[m][i] * DB0Vol[m][j];
				K_sri[i][j] += val * detJ0 * 8.0;
			}
		}  // detJ0 > 1e-15
	}

	// 2x2x2 Gauss quadrature — K_full for hourglass stabilization
	double gp = 1.0 / sqrt(3.0);
	double gauss_pts[2] = {-gp, gp};
	double K_full[24][24];
	for (int i = 0; i < 24; i++)
		for (int j = 0; j < 24; j++)
			K_full[i][j] = 0.0;

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

				// K_full += B^T * D * B * detJ (weight=1 at each Gauss point)
				double DB[6][24], DBDev[6][24];
				for (int i = 0; i < 6; i++)
					for (int j = 0; j < 24; j++)
					{
						DB[i][j] = 0.0;
						DBDev[i][j] = 0.0;
						for (int m = 0; m < 6; m++)
						{
							DB[i][j] += D[i][m] * B[m][j];
							DBDev[i][j] += D_dev[i][m] * B[m][j];
						}
					}

				for (int i = 0; i < 24; i++)
					for (int j = i; j < 24; j++)
					{
						double val = 0.0, val_dev = 0.0;
						for (int m = 0; m < 6; m++)
						{
							val += B[m][i] * DB[m][j];
							val_dev += B[m][i] * DBDev[m][j];
						}
						K_full[i][j] += val * detJ;
						K_sri[i][j] += val_dev * detJ;
					}
			}
		}
	}

	// Combine: K = K_1pt + alpha * (K_full - K_1pt)
	// Only use reduced integration if 1-pt stiffness is well-conditioned
	// Adaptive stabilization: slender elements get reduced alpha
	// Estimate aspect ratio from node coordinate ranges
	double xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
	for (int n = 1; n < 8; n++) {
		if (x[n] < xmin) xmin = x[n]; if (x[n] > xmax) xmax = x[n];
		if (y[n] < ymin) ymin = y[n]; if (y[n] > ymax) ymax = y[n];
		if (z[n] < zmin) zmin = z[n]; if (z[n] > zmax) zmax = z[n];
	}
	double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
	double hmax = fmax(fmax(dx, dy), dz);
	double hmin = fmin(fmin(dx, dy), dz);
	double aspect = (hmin > 1e-20) ? hmax / hmin : 1.0;

	double sum_diag_1pt = 0.0, sum_diag_full = 0.0;
	for (int i = 0; i < 24; i++) {
		sum_diag_1pt += fabs(K[i][i]);
		sum_diag_full += fabs(K_full[i][i]);
	}
	if (sum_diag_1pt > 1e-10 && sum_diag_full > 1e-10) {
		if (UseOrthogonalHourglass())
		{
			double volume = fmax(fabs(dx * dy * dz), 8.0 * fabs(detJ0_center));
			double characteristic = cbrt(volume);
			double shear = E / (2.0 * (1.0 + nu));
			double hg_scale = OrthogonalHourglassScale();
			if (aspect > 3.0)
			{
				hg_scale *= 3.0 / aspect;
				double alpha_min = HourglassAlphaMin();
				if (hg_scale < alpha_min) hg_scale = alpha_min;
			}
			double k_hg = hg_scale * shear * volume / (characteristic * characteristic);
			double xc = 0.0, yc = 0.0, zc = 0.0;
			for (int n = 0; n < 8; n++)
			{
				xc += x[n];
				yc += y[n];
				zc += z[n];
			}
			xc /= 8.0;
			yc /= 8.0;
			zc /= 8.0;

			double phys[4][8];
			for (int n = 0; n < 8; n++)
			{
				phys[0][n] = 1.0;
				phys[1][n] = x[n] - xc;
				phys[2][n] = y[n] - yc;
				phys[3][n] = z[n] - zc;
			}
			double orth_phys[4][8];
			int phys_count = 0;
			for (int m = 0; m < 4; m++)
			{
				double q[8];
				for (int n = 0; n < 8; n++)
					q[n] = phys[m][n];
				for (int p = 0; p < phys_count; p++)
				{
					double dot = 0.0;
					for (int n = 0; n < 8; n++)
						dot += q[n] * orth_phys[p][n];
					for (int n = 0; n < 8; n++)
						q[n] -= dot * orth_phys[p][n];
				}
				double norm = 0.0;
				for (int n = 0; n < 8; n++)
					norm += q[n] * q[n];
				if (norm <= 1.0e-20)
					continue;
				norm = sqrt(norm);
				for (int n = 0; n < 8; n++)
					orth_phys[phys_count][n] = q[n] / norm;
				phys_count++;
			}

			const double hg_seed[4][8] = {
				{ 1, -1,  1, -1,  1, -1,  1, -1},
				{ 1,  1, -1, -1, -1, -1,  1,  1},
				{ 1, -1, -1,  1, -1,  1,  1, -1},
				{-1,  1, -1,  1,  1, -1,  1, -1}
			};
			for (int mode = 0; mode < 4; mode++)
			{
				double q[8];
				for (int n = 0; n < 8; n++)
					q[n] = hg_seed[mode][n];
				for (int p = 0; p < phys_count; p++)
				{
					double dot = 0.0;
					for (int n = 0; n < 8; n++)
						dot += q[n] * orth_phys[p][n];
					for (int n = 0; n < 8; n++)
						q[n] -= dot * orth_phys[p][n];
				}
				double norm = 0.0;
				for (int n = 0; n < 8; n++)
					norm += q[n] * q[n];
				if (norm <= 1.0e-20)
					continue;
				norm = sqrt(norm);
				for (int n = 0; n < 8; n++)
					q[n] /= norm;
				for (int a = 0; a < 8; a++)
					for (int b = a; b < 8; b++)
						for (int d = 0; d < 3; d++)
						{
							double direction_scale = OrthogonalHourglassDirectionScale(d);
							if (direction_scale < 0.0)
								direction_scale = 0.0;
							int i = 3 * a + d;
							int j = 3 * b + d;
							K[i][j] += direction_scale * k_hg * q[a] * q[b];
						}
			}
		}
		else
		{
			double alpha_base = HourglassAlphaBase();
			double alpha = alpha_base;
			if (aspect > 3.0) {
				alpha = alpha_base * 3.0 / aspect;
				double alpha_min = HourglassAlphaMin();
				if (alpha < alpha_min) alpha = alpha_min;
			}
			for (int i = 0; i < 24; i++)
				for (int j = i; j < 24; j++)
					K[i][j] = K[i][j] + alpha * (K_full[i][j] - K[i][j]);
		}
	} else {
		// Fall back to full integration for degenerate elements
		for (int i = 0; i < 24; i++)
			for (int j = i; j < 24; j++)
				K[i][j] = K_full[i][j];
	}

	if (UseSelectiveReducedIntegration())
	{
		double blend = SRIHourglassBlend();
		if (blend < 0.0) blend = 0.0;
		if (blend > 1.0) blend = 1.0;
		for (int i = 0; i < 24; i++)
			for (int j = i; j < 24; j++)
				K[i][j] = K_sri[i][j] + blend * (K_full[i][j] - K_sri[i][j]);
	}

	double fb_scale = FBHourglassScale();
	if (fb_scale > 0.0)
	{
		double volume = fmax(fabs(dx * dy * dz), 8.0 * fabs(detJ0_center));
		double characteristic = cbrt(volume);
		double shear = E / (2.0 * (1.0 + nu));
		double k_hg = fb_scale * shear * volume / (characteristic * characteristic);
		const double gamma[4][8] = {
			{ 1, -1,  1, -1,  1, -1,  1, -1},
			{ 1,  1, -1, -1, -1, -1,  1,  1},
			{ 1, -1, -1,  1, -1,  1,  1, -1},
			{-1,  1, -1,  1,  1, -1,  1, -1}
		};
		for (int mode = 0; mode < 4; mode++)
			for (int a = 0; a < 8; a++)
				for (int b = a; b < 8; b++)
					for (int d = 0; d < 3; d++)
					{
						int i = 3 * a + d;
						int j = 3 * b + d;
						K[i][j] += k_hg * gamma[mode][a] * gamma[mode][b] / 8.0;
					}
	}

	double y_scale = GlobalYStiffnessScale();
	if (y_scale < 1.0e-6)
		y_scale = 1.0e-6;
	if (fabs(y_scale - 1.0) > 1.0e-12)
	{
		double y_factor = sqrt(y_scale);
		for (int i = 0; i < 24; i++)
		{
			double si = (i % 3 == 1) ? y_factor : 1.0;
			for (int j = i; j < 24; j++)
			{
				double sj = (j % 3 == 1) ? y_factor : 1.0;
				K[i][j] *= si * sj;
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
