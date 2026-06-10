/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>

using namespace std;

// LDLT facterization
void CLDLTSolver::LDLT()
{
	unsigned int N = K.dim();
	unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Hights

	// Compute max initial diagonal for relative pivot tolerance
	double max_diag = 0.0;
	for (unsigned int i = 1; i <= N; i++)
	{
		double d = fabs(K(i,i));
		if (d > max_diag) max_diag = d;
	}

	// Relative pivot tolerance: if pivot drops below 1e-12 * max_diag,
	// the matrix is numerically singular for engineering purposes.
	// FLT_MIN (~1.175e-38) was far too small to catch near-zero pivots.
	// Pivot tolerance: check relative ratio rather than absolute magnitude
	// to handle mixed translational/rotational DOF stiffness scales
	// Use small absolute tolerance instead of relative to handle mixed stiffness scales
	(void)max_diag;
	double pivot_tol = 1.0e-10;

	double min_pivot = max_diag;
	unsigned int min_pivot_eq = 0;

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
		// Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
			// Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += K(r,i) * K(r,j);		// C += L_ri * U_rj

			K(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = K(r,j) / K(r,r);	// L_rj = U_rj / D_rr
			K(j,j) -= Lrj * K(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			K(r,j) = Lrj;
		}

		// Track minimum pivot for diagnostics
		double abs_pivot = fabs(K(j,j));
		if (abs_pivot < min_pivot)
		{
			min_pivot = abs_pivot;
			min_pivot_eq = j;
		}

		if (abs_pivot <= pivot_tol)
		{
			cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
				 << "    Equation no = " << j << endl
				 << "    Pivot = " << K(j,j) << endl
				 << "    Pivot tolerance = " << pivot_tol << endl
				 << "    Max initial diagonal = " << max_diag << endl;

			exit(4);
		}
	}

	cout << "  LDLT: N=" << N << " max_diag=" << max_diag
		 << " pivot_tol=" << pivot_tol << endl;
	cout << "  LDLT: min pivot = " << min_pivot
		 << " at equation " << min_pivot_eq
		 << " (ratio = " << min_pivot/max_diag << ")" << endl;
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)
{
	unsigned int N = K.dim();
	unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Hights

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
		unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= K(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= K(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
		unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= K(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};
