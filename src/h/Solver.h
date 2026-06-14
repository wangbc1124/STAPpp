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

#include "SkylineMatrix.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

class CLinearOperator
{
public:
    virtual ~CLinearOperator() {}
    virtual unsigned int dim() const = 0;
    virtual void MatVec(const std::vector<double>& x, std::vector<double>& y) const = 0;
};

class CSparseSymmetricMatrix : public CLinearOperator
{
private:
    unsigned int N_;
    std::vector<std::vector<unsigned int> > pattern_;
    std::vector<unsigned int> column_ptr_;
    std::vector<unsigned int> row_ind_;
    std::vector<double> values_;
    unsigned int max_half_bandwidth_;

public:
    struct SymmetryDiagnostic
    {
        double symmetry_error;
        bool numerically_symmetric;
    };

    explicit CSparseSymmetricMatrix(unsigned int N = 0);

    void AddPattern(unsigned int row, unsigned int col);
    void AddPattern(unsigned int* location_matrix, std::size_t nd);
    void FinalizePattern();
    void AddValue(unsigned int row, unsigned int col, double value);
    unsigned long long Assembly(double* element_matrix, unsigned int* location_matrix, std::size_t nd);
    void MatVec(const std::vector<double>& x, std::vector<double>& y) const;
    std::vector<double> Diagonal() const;
    double DiagonalMaxAbs() const;
    SymmetryDiagnostic CheckSymmetry(double tolerance = 1.0e-10) const;
    void CheckFinite(const std::string& context) const;
    double RelativeResidual(const double* rhs, const double* solution) const;
    void ExportUpperCSR(std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& a) const;
    double AverageColumnNNZ() const;
    std::size_t MaxColumnNNZ() const;

    inline unsigned int dim() const { return N_; }
    inline std::size_t nnz() const { return values_.size(); }
    inline unsigned int GetMaximumHalfBandwidth() const { return max_half_bandwidth_; }
};

//! Sparse matrix in full CSR storage. Indices are zero based.
class CCSRMatrix : public CLinearOperator
{
public:
    enum StorageMode
    {
        kFull,
        kSymmetricUpper
    };

private:
    unsigned int N_;
    StorageMode storage_mode_;
    std::vector<std::vector<unsigned int> > pattern_;
    std::vector<std::size_t> row_ptr_;
    std::vector<unsigned int> col_ind_;
    std::vector<double> values_;
    std::vector<unsigned int> fast_index_rows_;
    std::vector<std::size_t> fast_index_table_offsets_;
    std::vector<std::size_t> fast_index_table_sizes_;
    std::vector<unsigned int> fast_index_keys_;
    std::vector<std::size_t> fast_index_values_;

public:
    struct SymmetryDiagnostic
    {
        double symmetry_error;
        bool numerically_symmetric;
    };

    explicit CCSRMatrix(unsigned int N = 0, StorageMode storage_mode = kFull);

    void AddPattern(unsigned int row, unsigned int col);
    void FinalizePattern();
    void AddValue(unsigned int row, unsigned int col, double value);
    void MatVec(const std::vector<double>& x, std::vector<double>& y) const;
    std::vector<double> Diagonal() const;
    double DiagonalMaxAbs() const;
    SymmetryDiagnostic CheckSymmetry(double tolerance = 1.0e-10) const;
    void CheckFinite(const std::string& context) const;
    double GetValue(unsigned int row, unsigned int col) const;
    double RelativeResidual(const double* rhs, const double* solution) const;
    double AverageRowNNZ() const;
    std::size_t MaxRowNNZ() const;

    inline std::size_t RowBegin(unsigned int row) const { return row_ptr_[row]; }
    inline std::size_t RowEnd(unsigned int row) const { return row_ptr_[row + 1]; }
    inline unsigned int ColumnIndex(std::size_t index) const { return col_ind_[index]; }
    inline double ValueAt(std::size_t index) const { return values_[index]; }

    inline unsigned int dim() const { return N_; }
    inline std::size_t nnz() const { return values_.size(); }
    inline StorageMode storage_mode() const { return storage_mode_; }
};

class CDiagonalScaling
{
public:
    struct Stats
    {
        double min_abs_diag;
        double max_abs_diag;
        double max_scale;
        double min_scale;
        unsigned int number_of_small_or_zero_diagonal_entries;
    };

private:
    bool enabled_;
    std::vector<double> scales_;
    Stats stats_;

public:
    CDiagonalScaling();

    void Disable();
    void Setup(const CCSRMatrix& matrix, double epsD = 1.0e-300);
    void ScaleRightHandSide(const double* rhs, std::vector<double>& scaled_rhs) const;
    void RecoverSolution(const std::vector<double>& scaled_solution, double* solution) const;

    inline bool enabled() const { return enabled_; }
    inline double Scale(unsigned int i) const { return enabled_ ? scales_[i] : 1.0; }
    inline const Stats& GetStats() const { return stats_; }
};

class CScaledOperator : public CLinearOperator
{
private:
    const CCSRMatrix& K_;
    const CDiagonalScaling* Scaling_;

public:
    explicit CScaledOperator(const CCSRMatrix& K, const CDiagonalScaling* scaling = nullptr);

    unsigned int dim() const;
    void MatVec(const std::vector<double>& x, std::vector<double>& y) const;
};

class CPreconditioner
{
public:
    virtual ~CPreconditioner() {}
    virtual void Apply(const std::vector<double>& r, std::vector<double>& z) const = 0;
    virtual const char* Name() const = 0;
};

class CIdentityPreconditioner : public CPreconditioner
{
public:
    void Apply(const std::vector<double>& r, std::vector<double>& z) const;
    const char* Name() const;
};

class CJacobiPreconditioner : public CPreconditioner
{
private:
    std::vector<double> inv_diag_;

public:
    CJacobiPreconditioner(const CCSRMatrix& matrix, const CDiagonalScaling* scaling);

    void Apply(const std::vector<double>& r, std::vector<double>& z) const;
    const char* Name() const;
};

class CSSORPreconditioner : public CPreconditioner
{
private:
    const CCSRMatrix& matrix_;
    const CDiagonalScaling* scaling_;
    std::vector<double> diag_;

public:
    CSSORPreconditioner(const CCSRMatrix& matrix, const CDiagonalScaling* scaling);

    void Apply(const std::vector<double>& r, std::vector<double>& z) const;
    const char* Name() const;
};

class CBlockJacobiPreconditioner : public CPreconditioner
{
public:
    struct Stats
    {
        unsigned int num_blocks;
        unsigned int block_size_histogram[6];
        unsigned int num_cholesky_success;
        unsigned int num_shifted_blocks;
        unsigned int num_diag_fallback_blocks;
    };

private:
    struct Block
    {
        std::vector<unsigned int> equations;
        std::vector<double> factor;
        std::vector<double> inv_diag;
        unsigned int size;
        bool diagonal_fallback;
    };

    std::vector<Block> blocks_;
    Stats stats_;

public:
    CBlockJacobiPreconditioner(const CCSRMatrix& matrix,
                               const CDiagonalScaling* scaling,
                               const std::vector<std::vector<unsigned int> >& equation_blocks);

    void Apply(const std::vector<double>& r, std::vector<double>& z) const;
    const char* Name() const;
    inline const Stats& GetStats() const { return stats_; }
};

class CILU0Preconditioner : public CPreconditioner
{
public:
    struct Stats
    {
        unsigned int shifted_pivots;
        unsigned int zero_pattern_fallbacks;
    };

private:
    unsigned int n_;
    std::vector<std::size_t> row_ptr_;
    std::vector<unsigned int> col_ind_;
    std::vector<double> lu_;
    Stats stats_;

public:
    CILU0Preconditioner(const CCSRMatrix& matrix, const CDiagonalScaling* scaling);

    void Apply(const std::vector<double>& r, std::vector<double>& z) const;
    const char* Name() const;
    inline const Stats& GetStats() const { return stats_; }
};

//! Preconditioned conjugate-gradient solver for symmetric positive definite CSR systems.
class CPCGSolver
{
public:
    struct Result
    {
        bool converged;
        unsigned int iterations;
        double initial_relative_residual;
        double final_relative_residual;
    };

private:
    const CLinearOperator& K_;
    const CPreconditioner& Preconditioner_;
    double tolerance_;
    unsigned int max_iterations_;
    unsigned int log_interval_;

public:
    CPCGSolver(const CLinearOperator& K, const CPreconditioner& preconditioner,
               double tolerance, unsigned int max_iterations, unsigned int log_interval);
    Result Solve(const double* rhs, double* solution) const;
};

//! BiCGSTAB solver for sparse systems that are difficult for CG.
class CBiCGSTABSolver
{
public:
    typedef CPCGSolver::Result Result;

private:
    const CLinearOperator& K_;
    const CPreconditioner& Preconditioner_;
    double tolerance_;
    unsigned int max_iterations_;
    unsigned int log_interval_;

public:
    CBiCGSTABSolver(const CLinearOperator& K, const CPreconditioner& preconditioner,
                    double tolerance, unsigned int max_iterations, unsigned int log_interval);
    Result Solve(const double* rhs, double* solution) const;
};

//! Restarted GMRES solver for sparse systems with stronger convergence control.
class CGMRESSolver
{
public:
    typedef CPCGSolver::Result Result;

private:
    const CLinearOperator& K_;
    const CPreconditioner& Preconditioner_;
    double tolerance_;
    unsigned int max_iterations_;
    unsigned int restart_;
    unsigned int log_interval_;

public:
    CGMRESSolver(const CLinearOperator& K, const CPreconditioner& preconditioner,
                 double tolerance, unsigned int max_iterations,
                 unsigned int restart, unsigned int log_interval);
    Result Solve(const double* rhs, double* solution) const;
};

struct CSparseSolverOptions
{
    std::string backend_name;
    std::string requested_solver;
    std::string requested_preconditioner;
    std::string pardiso_mtype_mode;
    std::string scale_mode;
    double tolerance;
    unsigned int max_iterations;
    unsigned int residual_print_interval;
    unsigned int gmres_restart;

    CSparseSolverOptions();
};

struct CSparseBackendSetupTimings
{
    double symmetry_check_time;
    double scaling_time;
    double preconditioner_setup_time;

    CSparseBackendSetupTimings();
};

struct CSparseBackendInfo
{
    struct PardisoInfo
    {
        bool enabled;
        std::string requested_mtype_mode;
        std::string selected_mtype_mode;
        int selected_mtype;
        std::string matrix_part;
        unsigned int attempt_count;
        bool retried_from_spd_to_sym_indef;
        int phase11_error;
        int phase22_error;
        int phase33_error;
        long long factor_nnz;
        long long peak_memory_kb;
        double export_upper_csr_time;
    };

    bool backend_available;
    std::string backend_name;
    std::string requested_solver;
    std::string actual_solver;
    std::string requested_preconditioner;
    std::string actual_preconditioner;
    std::string scale_mode;
    CCSRMatrix::SymmetryDiagnostic symmetry;
    bool scaling_enabled;
    CDiagonalScaling::Stats scaling_stats;
    bool has_block_jacobi_stats;
    CBlockJacobiPreconditioner::Stats block_jacobi_stats;
    bool has_ilu0_stats;
    CILU0Preconditioner::Stats ilu0_stats;
    PardisoInfo pardiso_info;
    CSparseBackendSetupTimings timings;

    CSparseBackendInfo();
};

struct CSparseSolveResult
{
    bool converged;
    unsigned int iterations;
    double initial_relative_residual;
    double final_relative_residual;

    CSparseSolveResult();
};

//! Sparse solver backend abstraction shared by the in-tree implementation and future external-library backends.
class CSparseSolverBackend
{
public:
    virtual ~CSparseSolverBackend() {}

    virtual void Setup(const CCSRMatrix& matrix,
                       const std::vector<std::vector<unsigned int> >* equation_blocks,
                       const CSparseSolverOptions& options) = 0;
    virtual const CSparseBackendInfo& GetInfo() const = 0;
    virtual CSparseSolveResult Solve(const double* rhs, double* solution) const = 0;
};

std::unique_ptr<CSparseSolverBackend> CreateSparseSolverBackend(const std::string& backend_name);
bool IsSparseBackendAvailable(const std::string& backend_name);
CSparseSolveResult SolvePardisoSystem(const CSparseSymmetricMatrix& matrix,
                                      const double* rhs,
                                      double* solution,
                                      const std::string& requested_mtype_mode,
                                      CSparseBackendInfo& info);

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver
{
private:
    
    CSkylineMatrix<double>& K;

public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K): K(*K) {};

//!	Perform L*D*L(T) factorization of the stiffness matrix
	void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	void BackSubstitution(double* Force); 
};
