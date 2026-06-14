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
#include "Outputter.h"
#include "Clock.h"

#ifdef STAPPP_ENABLE_MKL_PARDISO
#include "mkl.h"
#include "mkl_pardiso.h"
#endif

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace
{
const std::size_t kCsrFastIndexRowThreshold = 256;

class CStandardSparseSolverBackend : public CSparseSolverBackend
{
private:
    const CCSRMatrix* Matrix_;
    const vector<vector<unsigned int> >* EquationBlocks_;
    CSparseSolverOptions Options_;
    CSparseBackendInfo Info_;
    CDiagonalScaling Scaling_;
    std::unique_ptr<CScaledOperator> ScaledOperator_;
    unique_ptr<CPreconditioner> Preconditioner_;

public:
    CStandardSparseSolverBackend();

    void Setup(const CCSRMatrix& matrix,
               const vector<vector<unsigned int> >* equation_blocks,
               const CSparseSolverOptions& options);
    const CSparseBackendInfo& GetInfo() const;
    CSparseSolveResult Solve(const double* rhs, double* solution) const;
};

class CPardisoSparseSolverBackend : public CSparseSolverBackend
{
private:
    const CCSRMatrix* Matrix_;
    CSparseSolverOptions Options_;
    mutable CSparseBackendInfo Info_;

public:
    CPardisoSparseSolverBackend();

    void Setup(const CCSRMatrix& matrix,
               const vector<vector<unsigned int> >* equation_blocks,
               const CSparseSolverOptions& options);
    const CSparseBackendInfo& GetInfo() const;
    CSparseSolveResult Solve(const double* rhs, double* solution) const;
};
}

static COutputter& solver_output()
{
    return *COutputter::GetInstance();
}

static void require_finite_scalar(double value, const string& context)
{
    if (!isfinite(value))
    {
        ostringstream msg;
        msg << "*** Error *** Non-finite number detected in " << context;
        throw runtime_error(msg.str());
    }
}

static void require_finite_vector(const vector<double>& values, const string& context)
{
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (!isfinite(values[i]))
        {
            ostringstream msg;
            msg << "*** Error *** Non-finite number detected in " << context
                << " at index " << i << ": " << values[i];
            throw runtime_error(msg.str());
        }
    }
}

static double dot_product(const vector<double>& a, const vector<double>& b)
{
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        sum += a[i] * b[i];
    require_finite_scalar(sum, "dot product");
    return sum;
}

static void print_stagnation_warning(const string& solver_name, unsigned int iter,
                                     vector<double>& window, double rel)
{
    window.push_back(rel);
    if (window.size() > 501)
        window.erase(window.begin());
    if (window.size() == 501)
    {
        const double previous = window.front();
        if (previous > 0.0 && rel > 0.95 * previous)
        {
            solver_output() << "  " << solver_name
                            << " warning: relative residual decreased less than 5% over the last 500 iterations"
                            << " at iteration " << iter << endl;
            window.clear();
            window.push_back(rel);
        }
    }
}

static bool dense_cholesky_factor(vector<double>& block, unsigned int size)
{
    for (unsigned int i = 0; i < size; ++i)
    {
        for (unsigned int j = 0; j <= i; ++j)
        {
            double sum = block[i * size + j];
            for (unsigned int k = 0; k < j; ++k)
                sum -= block[i * size + k] * block[j * size + k];

            if (i == j)
            {
                if (!isfinite(sum) || sum <= DBL_MIN)
                    return false;
                block[i * size + i] = sqrt(sum);
            }
            else
            {
                block[i * size + j] = sum / block[j * size + j];
            }
        }
        for (unsigned int j = i + 1; j < size; ++j)
            block[i * size + j] = 0.0;
    }
    return true;
}

#ifdef STAPPP_ENABLE_MKL_PARDISO
static void pardiso_log_phase(const char* label, MKL_INT phase, MKL_INT error)
{
    solver_output() << "  PARDISO " << label
                    << " phase=" << phase
                    << " error=" << error << endl;
}

struct PardisoMatrixData
{
    std::vector<MKL_INT> ia;
    std::vector<MKL_INT> ja;
    std::vector<double> a;
};

struct PardisoAttemptResult
{
    int phase11_error;
    int phase22_error;
    int phase33_error;
    long long factor_nnz;
    long long peak_memory_kb;
};

static void build_pardiso_csr(const CCSRMatrix& matrix, bool upper_only,
                              PardisoMatrixData& data)
{
    const unsigned int n = matrix.dim();
    data.ia.assign(static_cast<size_t>(n) + 1, 1);
    data.ja.clear();
    data.a.clear();

    size_t reserve_nnz = matrix.nnz();
    if (upper_only)
        reserve_nnz = reserve_nnz / 2 + n;
    data.ja.reserve(reserve_nnz);
    data.a.reserve(reserve_nnz);

    for (unsigned int row = 0; row < n; ++row)
    {
        data.ia[row] = static_cast<MKL_INT>(data.ja.size()) + 1;
        for (size_t idx = matrix.RowBegin(row); idx < matrix.RowEnd(row); ++idx)
        {
            const unsigned int col = matrix.ColumnIndex(idx);
            if (upper_only && col < row)
                continue;

            const double value = matrix.ValueAt(idx);
            require_finite_scalar(value, "PARDISO CSR matrix values");
            data.ja.push_back(static_cast<MKL_INT>(col) + 1);
            data.a.push_back(value);
        }
    }

    data.ia[n] = static_cast<MKL_INT>(data.ja.size()) + 1;
}

static void build_pardiso_csr(const CSparseSymmetricMatrix& matrix,
                              PardisoMatrixData& data,
                              double* export_time_seconds)
{
    Clock timer;
    timer.Start();
    vector<int> ia_int;
    vector<int> ja_int;
    vector<double> a_double;
    matrix.ExportUpperCSR(ia_int, ja_int, a_double);

    data.ia.assign(ia_int.begin(), ia_int.end());
    data.ja.assign(ja_int.begin(), ja_int.end());
    data.a.swap(a_double);
    for (size_t i = 0; i < data.a.size(); ++i)
        require_finite_scalar(data.a[i], "PARDISO symmetric CSR matrix values");
    if (export_time_seconds)
        *export_time_seconds = timer.ElapsedTime();
}

static void release_pardiso_memory(void* pt[64], MKL_INT maxfct, MKL_INT mnum,
                                   MKL_INT mtype, MKL_INT n, MKL_INT iparm[64],
                                   MKL_INT msglvl, PardisoMatrixData& matrix_data,
                                   std::vector<double>& rhs, std::vector<double>& solution)
{
    MKL_INT phase = -1;
    MKL_INT nrhs = 1;
    MKL_INT error = 0;
    pardiso_log_phase("begin", phase, 0);
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n,
            matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
            &nrhs, iparm, &msglvl, rhs.data(), solution.data(), &error);
    pardiso_log_phase("end", phase, error);
}

static PardisoAttemptResult run_pardiso_attempt(const CCSRMatrix& matrix,
                                                const double* rhs,
                                                int mtype,
                                                bool upper_only,
                                                double* solution)
{
    const MKL_INT n = static_cast<MKL_INT>(matrix.dim());
    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT iparm[64] = {0};
    void* pt[64] = {0};
    PardisoMatrixData matrix_data;
    std::vector<double> rhs_copy(rhs, rhs + matrix.dim());
    std::vector<double> x(matrix.dim(), 0.0);
    PardisoAttemptResult attempt;
    attempt.phase11_error = 0;
    attempt.phase22_error = 0;
    attempt.phase33_error = 0;
    attempt.factor_nnz = -1;
    attempt.peak_memory_kb = -1;

    require_finite_vector(rhs_copy, "PARDISO RHS");
    build_pardiso_csr(matrix, upper_only, matrix_data);

    iparm[0] = 1;
    iparm[1] = 0;
    iparm[7] = 0;
    iparm[9] = 8;
    iparm[17] = -1;
    iparm[34] = 0;

    MKL_INT phase = 11;
    MKL_INT error = 0;
    bool need_release = false;

    try
    {
        pardiso_log_phase("begin", phase, 0);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, const_cast<MKL_INT*>(&n),
                matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
                &nrhs, iparm, &msglvl, rhs_copy.data(), x.data(), &error);
        pardiso_log_phase("end", phase, error);
        attempt.phase11_error = error;
        if (error != 0)
        {
            need_release = true;
            release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                                   msglvl, matrix_data, rhs_copy, x);
            return attempt;
        }
        need_release = true;

        phase = 22;
        error = 0;
        pardiso_log_phase("begin", phase, 0);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, const_cast<MKL_INT*>(&n),
                matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
                &nrhs, iparm, &msglvl, rhs_copy.data(), x.data(), &error);
        pardiso_log_phase("end", phase, error);
        attempt.phase22_error = error;
        attempt.factor_nnz = iparm[17];
        if (error != 0)
        {
            release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                                   msglvl, matrix_data, rhs_copy, x);
            return attempt;
        }

        phase = 33;
        error = 0;
        pardiso_log_phase("begin", phase, 0);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, const_cast<MKL_INT*>(&n),
                matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
                &nrhs, iparm, &msglvl, rhs_copy.data(), x.data(), &error);
        pardiso_log_phase("end", phase, error);
        attempt.phase33_error = error;
        if (error == 0)
        {
            require_finite_vector(x, "PARDISO solution");
            for (unsigned int i = 0; i < matrix.dim(); ++i)
                solution[i] = x[i];
        }

        release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                               msglvl, matrix_data, rhs_copy, x);
        return attempt;
    }
    catch (...)
    {
        if (need_release)
        {
            try
            {
                release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                                       msglvl, matrix_data, rhs_copy, x);
            }
            catch (...)
            {
            }
        }
        throw;
    }
}

static PardisoAttemptResult run_pardiso_attempt(const CSparseSymmetricMatrix& matrix,
                                                const double* rhs,
                                                int mtype,
                                                double* export_time_seconds,
                                                double* solution)
{
    const MKL_INT n = static_cast<MKL_INT>(matrix.dim());
    MKL_INT maxfct = 1;
    MKL_INT mnum = 1;
    MKL_INT nrhs = 1;
    MKL_INT msglvl = 0;
    MKL_INT iparm[64] = {0};
    void* pt[64] = {0};
    PardisoMatrixData matrix_data;
    std::vector<double> rhs_copy(rhs, rhs + matrix.dim());
    std::vector<double> x(matrix.dim(), 0.0);
    PardisoAttemptResult attempt;
    attempt.phase11_error = 0;
    attempt.phase22_error = 0;
    attempt.phase33_error = 0;
    attempt.factor_nnz = -1;
    attempt.peak_memory_kb = -1;

    require_finite_vector(rhs_copy, "PARDISO RHS");
    build_pardiso_csr(matrix, matrix_data, export_time_seconds);

    iparm[0] = 1;
    iparm[1] = 2;
    iparm[7] = 2;
    iparm[9] = 13;
    iparm[10] = 1;
    iparm[12] = 1;
    iparm[17] = -1;
    iparm[18] = -1;
    iparm[26] = 1;
    iparm[34] = 0;

    MKL_INT phase = 11;
    MKL_INT error = 0;
    bool need_release = false;

    try
    {
        pardiso_log_phase("begin", phase, 0);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, const_cast<MKL_INT*>(&n),
                matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
                &nrhs, iparm, &msglvl, rhs_copy.data(), x.data(), &error);
        pardiso_log_phase("end", phase, error);
        attempt.phase11_error = error;
        if (error != 0)
        {
            need_release = true;
            release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                                   msglvl, matrix_data, rhs_copy, x);
            return attempt;
        }
        need_release = true;

        phase = 22;
        error = 0;
        pardiso_log_phase("begin", phase, 0);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, const_cast<MKL_INT*>(&n),
                matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
                &nrhs, iparm, &msglvl, rhs_copy.data(), x.data(), &error);
        pardiso_log_phase("end", phase, error);
        attempt.phase22_error = error;
        attempt.factor_nnz = iparm[17];
        if (error != 0)
        {
            release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                                   msglvl, matrix_data, rhs_copy, x);
            return attempt;
        }

        phase = 33;
        error = 0;
        pardiso_log_phase("begin", phase, 0);
        pardiso(pt, &maxfct, &mnum, &mtype, &phase, const_cast<MKL_INT*>(&n),
                matrix_data.a.data(), matrix_data.ia.data(), matrix_data.ja.data(), 0,
                &nrhs, iparm, &msglvl, rhs_copy.data(), x.data(), &error);
        pardiso_log_phase("end", phase, error);
        attempt.phase33_error = error;
        if (error == 0)
        {
            require_finite_vector(x, "PARDISO solution");
            for (unsigned int i = 0; i < matrix.dim(); ++i)
                solution[i] = x[i];
        }

        release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                               msglvl, matrix_data, rhs_copy, x);
        return attempt;
    }
    catch (...)
    {
        if (need_release)
        {
            try
            {
                release_pardiso_memory(pt, maxfct, mnum, static_cast<MKL_INT>(mtype), n, iparm,
                                       msglvl, matrix_data, rhs_copy, x);
            }
            catch (...)
            {
            }
        }
        throw;
    }
}
#endif

static void dense_cholesky_solve(const vector<double>& factor, unsigned int size,
                                 const vector<double>& rhs, vector<double>& x)
{
    x.assign(size, 0.0);
    vector<double> y(size, 0.0);
    for (unsigned int i = 0; i < size; ++i)
    {
        double sum = rhs[i];
        for (unsigned int k = 0; k < i; ++k)
            sum -= factor[i * size + k] * y[k];
        y[i] = sum / factor[i * size + i];
    }
    for (int ii = static_cast<int>(size) - 1; ii >= 0; --ii)
    {
        const unsigned int i = static_cast<unsigned int>(ii);
        double sum = y[i];
        for (unsigned int k = i + 1; k < size; ++k)
            sum -= factor[k * size + i] * x[k];
        x[i] = sum / factor[i * size + i];
    }
}

CSparseSymmetricMatrix::CSparseSymmetricMatrix(unsigned int N)
    : N_(N), pattern_(N), max_half_bandwidth_(0)
{
}

void CSparseSymmetricMatrix::AddPattern(unsigned int row, unsigned int col)
{
    if (row >= N_ || col >= N_)
        throw out_of_range("symmetric sparse pattern index out of range");
    if (row > col)
        std::swap(row, col);
    pattern_[col].push_back(row);
}

void CSparseSymmetricMatrix::AddPattern(unsigned int* location_matrix, std::size_t nd)
{
    for (std::size_t j = 0; j < nd; ++j)
    {
        const unsigned int Lj = location_matrix[j];
        if (!Lj)
            continue;

        for (std::size_t i = 0; i <= j; ++i)
        {
            const unsigned int Li = location_matrix[i];
            if (!Li)
                continue;

            const unsigned int row = (Li < Lj) ? Li : Lj;
            const unsigned int col = (Li < Lj) ? Lj : Li;
            pattern_[col - 1].push_back(row - 1);
        }
    }
}

void CSparseSymmetricMatrix::FinalizePattern()
{
    column_ptr_.assign(N_ + 1, 0);
    max_half_bandwidth_ = 0;

    for (unsigned int col = 0; col < N_; ++col)
    {
        vector<unsigned int>& rows = pattern_[col];
        rows.push_back(col);
        sort(rows.begin(), rows.end());
        rows.erase(unique(rows.begin(), rows.end()), rows.end());
        column_ptr_[col + 1] = column_ptr_[col] + static_cast<unsigned int>(rows.size());
        if (!rows.empty())
        {
            const unsigned int height = col - rows.front() + 1;
            if (height > max_half_bandwidth_)
                max_half_bandwidth_ = height;
        }
    }

    row_ind_.resize(column_ptr_[N_]);
    values_.assign(column_ptr_[N_], 0.0);
    for (unsigned int col = 0; col < N_; ++col)
    {
        const vector<unsigned int>& rows = pattern_[col];
        copy(rows.begin(), rows.end(), row_ind_.begin() + column_ptr_[col]);
    }

    vector<vector<unsigned int> >().swap(pattern_);
}

void CSparseSymmetricMatrix::AddValue(unsigned int row, unsigned int col, double value)
{
    if (row >= N_ || col >= N_)
        throw out_of_range("symmetric sparse value index out of range");
    if (row > col)
        std::swap(row, col);

    const unsigned int begin = column_ptr_[col];
    const unsigned int end = column_ptr_[col + 1];
    vector<unsigned int>::const_iterator it =
        lower_bound(row_ind_.begin() + begin, row_ind_.begin() + end, row);
    if (it == row_ind_.begin() + end || *it != row)
        throw runtime_error("symmetric sparse value inserted outside finalized pattern");
    values_[static_cast<size_t>(it - row_ind_.begin())] += value;
}

unsigned long long CSparseSymmetricMatrix::Assembly(double* element_matrix, unsigned int* location_matrix, std::size_t nd)
{
    unsigned long long insertions = 0;
    for (std::size_t j = 0; j < nd; ++j)
    {
        const unsigned int Lj = location_matrix[j];
        if (!Lj)
            continue;

        const unsigned int diagj_element = static_cast<unsigned int>((j + 1) * j / 2);
        for (std::size_t i = 0; i <= j; ++i)
        {
            const unsigned int Li = location_matrix[i];
            if (!Li)
                continue;

            const double value = element_matrix[diagj_element + j - i];

            const unsigned int row = (Li < Lj) ? Li : Lj;
            const unsigned int col = (Li < Lj) ? Lj : Li;
            const unsigned int begin = column_ptr_[col - 1];
            const unsigned int end = column_ptr_[col];
            vector<unsigned int>::const_iterator it =
                lower_bound(row_ind_.begin() + begin, row_ind_.begin() + end, row - 1);
            if (it == row_ind_.begin() + end || *it != row - 1)
                throw runtime_error("symmetric sparse value inserted outside finalized pattern");
            values_[static_cast<size_t>(it - row_ind_.begin())] += value;
            ++insertions;
        }
    }
    return insertions;
}

void CSparseSymmetricMatrix::MatVec(const vector<double>& x, vector<double>& y) const
{
    require_finite_vector(x, "symmetric sparse matvec input");
    y.assign(N_, 0.0);
    for (unsigned int col = 0; col < N_; ++col)
    {
        for (unsigned int k = column_ptr_[col]; k < column_ptr_[col + 1]; ++k)
        {
            const unsigned int row = row_ind_[k];
            const double value = values_[k];
            if (row == col)
            {
                y[col] += value * x[col];
            }
            else
            {
                y[row] += value * x[col];
                y[col] += value * x[row];
            }
        }
        require_finite_scalar(y[col], "symmetric sparse matvec output");
    }
}

vector<double> CSparseSymmetricMatrix::Diagonal() const
{
    vector<double> diagonal(N_, 0.0);
    for (unsigned int col = 0; col < N_; ++col)
    {
        const unsigned int begin = column_ptr_[col];
        const unsigned int end = column_ptr_[col + 1];
        vector<unsigned int>::const_iterator it =
            lower_bound(row_ind_.begin() + begin, row_ind_.begin() + end, col);
        if (it != row_ind_.begin() + end && *it == col)
            diagonal[col] = values_[static_cast<size_t>(it - row_ind_.begin())];
    }
    return diagonal;
}

double CSparseSymmetricMatrix::DiagonalMaxAbs() const
{
    vector<double> diagonal = Diagonal();
    double max_diag = 0.0;
    for (size_t i = 0; i < diagonal.size(); ++i)
        max_diag = max(max_diag, fabs(diagonal[i]));
    return max_diag;
}

CSparseSymmetricMatrix::SymmetryDiagnostic CSparseSymmetricMatrix::CheckSymmetry(double tolerance) const
{
    (void)tolerance;
    CheckFinite("symmetric sparse matrix");
    SymmetryDiagnostic diagnostic;
    diagnostic.symmetry_error = 0.0;
    diagnostic.numerically_symmetric = true;
    return diagnostic;
}

void CSparseSymmetricMatrix::CheckFinite(const string& context) const
{
    require_finite_vector(values_, context);
}

double CSparseSymmetricMatrix::RelativeResidual(const double* rhs, const double* solution) const
{
    const unsigned int N = dim();
    vector<double> x(N), b(N), Ax(N);
    for (unsigned int i = 0; i < N; ++i)
    {
        x[i] = solution[i];
        b[i] = rhs[i];
    }
    require_finite_vector(b, "checked residual RHS");
    require_finite_vector(x, "checked residual solution");
    MatVec(x, Ax);
    double norm_r = 0.0;
    double norm_b = 0.0;
    for (unsigned int i = 0; i < N; ++i)
    {
        const double diff = Ax[i] - b[i];
        norm_r += diff * diff;
        norm_b += b[i] * b[i];
    }
    return (norm_b > 0.0) ? sqrt(norm_r / norm_b) : sqrt(norm_r);
}

void CSparseSymmetricMatrix::ExportUpperCSR(vector<int>& ia, vector<int>& ja, vector<double>& a) const
{
    ia.assign(N_ + 1, 1);
    ja.assign(values_.size(), 0);
    a.assign(values_.size(), 0.0);

    for (unsigned int col = 0; col < N_; ++col)
    {
        for (unsigned int k = column_ptr_[col]; k < column_ptr_[col + 1]; ++k)
            ++ia[row_ind_[k] + 1];
    }

    for (unsigned int row = 1; row <= N_; ++row)
        ia[row] += ia[row - 1] - 1;

    vector<int> next = ia;
    for (unsigned int row = 0; row < next.size(); ++row)
        --next[row];

    for (unsigned int col = 0; col < N_; ++col)
    {
        for (unsigned int k = column_ptr_[col]; k < column_ptr_[col + 1]; ++k)
        {
            const unsigned int row = row_ind_[k];
            const int dest = next[row]++;
            ja[dest] = static_cast<int>(col + 1);
            a[dest] = values_[k];
        }
    }
}

double CSparseSymmetricMatrix::AverageColumnNNZ() const
{
    if (N_ == 0)
        return 0.0;
    return static_cast<double>(values_.size()) / static_cast<double>(N_);
}

std::size_t CSparseSymmetricMatrix::MaxColumnNNZ() const
{
    std::size_t max_nnz = 0;
    for (unsigned int col = 0; col < N_; ++col)
    {
        const std::size_t col_nnz = column_ptr_[col + 1] - column_ptr_[col];
        if (col_nnz > max_nnz)
            max_nnz = col_nnz;
    }
    return max_nnz;
}

CCSRMatrix::CCSRMatrix(unsigned int N, StorageMode storage_mode)
    : N_(N), storage_mode_(storage_mode), pattern_(N)
{
}

void CCSRMatrix::AddPattern(unsigned int row, unsigned int col)
{
    if (row >= N_ || col >= N_)
        throw out_of_range("CSR pattern index out of range");
    if (storage_mode_ == kSymmetricUpper && col < row)
        std::swap(row, col);
    pattern_[row].push_back(col);
}

void CCSRMatrix::FinalizePattern()
{
    row_ptr_.assign(static_cast<size_t>(N_) + 1, 0);
    fast_index_rows_.clear();
    fast_index_table_offsets_.clear();
    fast_index_table_sizes_.clear();
    fast_index_keys_.clear();
    fast_index_values_.clear();

    for (unsigned int row = 0; row < N_; ++row)
    {
        vector<unsigned int>& cols = pattern_[row];
        sort(cols.begin(), cols.end());
        cols.erase(unique(cols.begin(), cols.end()), cols.end());
        row_ptr_[row + 1] = row_ptr_[row] + cols.size();
    }

    col_ind_.assign(row_ptr_[N_], 0);
    values_.assign(row_ptr_[N_], 0.0);
    for (unsigned int row = 0; row < N_; ++row)
    {
        const vector<unsigned int>& cols = pattern_[row];
        for (size_t k = 0; k < cols.size(); ++k)
            col_ind_[row_ptr_[row] + k] = cols[k];

        if (cols.size() >= kCsrFastIndexRowThreshold)
        {
            fast_index_rows_.push_back(row);
            const std::size_t table_size = cols.size() * 2 + 1;
            fast_index_table_offsets_.push_back(fast_index_keys_.size());
            fast_index_table_sizes_.push_back(table_size);
            fast_index_keys_.insert(fast_index_keys_.end(), table_size, std::numeric_limits<unsigned int>::max());
            fast_index_values_.insert(fast_index_values_.end(), table_size, static_cast<std::size_t>(-1));
            for (size_t k = 0; k < cols.size(); ++k)
            {
                const unsigned int key = cols[k];
                const std::size_t slot = row_ptr_[row] + k;
                std::size_t probe = static_cast<std::size_t>(key) % table_size;
                while (fast_index_keys_[fast_index_table_offsets_.back() + probe] != std::numeric_limits<unsigned int>::max())
                    probe = (probe + 1) % table_size;
                fast_index_keys_[fast_index_table_offsets_.back() + probe] = key;
                fast_index_values_[fast_index_table_offsets_.back() + probe] = slot;
            }
        }
    }
    vector<vector<unsigned int> >().swap(pattern_);
}

void CCSRMatrix::AddValue(unsigned int row, unsigned int col, double value)
{
    require_finite_scalar(value, "CSR value insertion");
    if (row >= N_ || col >= N_)
        throw out_of_range("CSR value index out of range");
    if (storage_mode_ == kSymmetricUpper && col < row)
        std::swap(row, col);

    vector<unsigned int>::const_iterator fast_row =
        lower_bound(fast_index_rows_.begin(), fast_index_rows_.end(), row);
    if (fast_row != fast_index_rows_.end() && *fast_row == row)
    {
        const size_t fast_pos = static_cast<size_t>(fast_row - fast_index_rows_.begin());
        const size_t table_offset = fast_index_table_offsets_[fast_pos];
        const size_t table_size = fast_index_table_sizes_[fast_pos];
        std::size_t probe = static_cast<std::size_t>(col) % table_size;
        while (fast_index_keys_[table_offset + probe] != std::numeric_limits<unsigned int>::max() &&
               fast_index_keys_[table_offset + probe] != col)
        {
            probe = (probe + 1) % table_size;
        }
        if (fast_index_keys_[table_offset + probe] != col)
            throw runtime_error("CSR fast index inconsistent with finalized pattern");
        const size_t fast_slot = fast_index_values_[table_offset + probe];
        values_[fast_slot] += value;
        return;
    }

    const size_t begin = row_ptr_[row];
    const size_t end = row_ptr_[row + 1];
    vector<unsigned int>::const_iterator it =
        lower_bound(col_ind_.begin() + begin, col_ind_.begin() + end, col);
    if (it == col_ind_.begin() + end || *it != col)
        throw runtime_error("CSR value inserted outside finalized pattern");
    values_[static_cast<size_t>(it - col_ind_.begin())] += value;
}

void CCSRMatrix::MatVec(const vector<double>& x, vector<double>& y) const
{
    require_finite_vector(x, "CSR matvec input");
    y.assign(N_, 0.0);
    if (storage_mode_ == kSymmetricUpper)
    {
        for (unsigned int row = 0; row < N_; ++row)
        {
            for (size_t k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
            {
                const unsigned int col = col_ind_[k];
                const double value = values_[k];
                y[row] += value * x[col];
                if (col != row)
                    y[col] += value * x[row];
            }
            require_finite_scalar(y[row], "CSR matvec output");
        }
        return;
    }

    for (unsigned int row = 0; row < N_; ++row)
    {
        double sum = 0.0;
        for (size_t k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
            sum += values_[k] * x[col_ind_[k]];
        require_finite_scalar(sum, "CSR matvec output");
        y[row] = sum;
    }
}

vector<double> CCSRMatrix::Diagonal() const
{
    vector<double> diag(N_, 0.0);
    for (unsigned int row = 0; row < N_; ++row)
    {
        for (size_t k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
        {
            if (col_ind_[k] == row)
            {
                diag[row] = values_[k];
                break;
            }
        }
    }
    return diag;
}

double CCSRMatrix::DiagonalMaxAbs() const
{
    double max_diag = 0.0;
    vector<double> diag = Diagonal();
    for (size_t i = 0; i < diag.size(); ++i)
        max_diag = max(max_diag, fabs(diag[i]));
    return max_diag;
}

double CCSRMatrix::AverageRowNNZ() const
{
    if (N_ == 0)
        return 0.0;
    return static_cast<double>(values_.size()) / static_cast<double>(N_);
}

std::size_t CCSRMatrix::MaxRowNNZ() const
{
    std::size_t max_row_nnz = 0;
    for (unsigned int row = 0; row < N_; ++row)
    {
        const std::size_t row_nnz = row_ptr_[row + 1] - row_ptr_[row];
        if (row_nnz > max_row_nnz)
            max_row_nnz = row_nnz;
    }
    return max_row_nnz;
}

double CCSRMatrix::GetValue(unsigned int row, unsigned int col) const
{
    if (row >= N_ || col >= N_)
        throw out_of_range("CSR value index out of range");
    const size_t begin = row_ptr_[row];
    const size_t end = row_ptr_[row + 1];
    vector<unsigned int>::const_iterator it =
        lower_bound(col_ind_.begin() + begin, col_ind_.begin() + end, col);
    if (it == col_ind_.begin() + end || *it != col)
        return 0.0;
    return values_[static_cast<size_t>(it - col_ind_.begin())];
}

CCSRMatrix::SymmetryDiagnostic CCSRMatrix::CheckSymmetry(double tolerance) const
{
    CheckFinite("CSR matrix");
    if (storage_mode_ == kSymmetricUpper)
    {
        SymmetryDiagnostic diagnostic;
        diagnostic.symmetry_error = 0.0;
        diagnostic.numerically_symmetric = true;
        return diagnostic;
    }

    double max_abs = 0.0;
    double max_diff = 0.0;

    for (unsigned int row = 0; row < N_; ++row)
    {
        for (size_t k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
        {
            const unsigned int col = col_ind_[k];
            const double aij = values_[k];
            max_abs = max(max_abs, fabs(aij));

            const size_t begin = row_ptr_[col];
            const size_t end = row_ptr_[col + 1];
            vector<unsigned int>::const_iterator it =
                lower_bound(col_ind_.begin() + begin, col_ind_.begin() + end, row);
            const double aji = (it != col_ind_.begin() + end && *it == row)
                ? values_[static_cast<size_t>(it - col_ind_.begin())]
                : 0.0;
            max_diff = max(max_diff, fabs(aij - aji));
        }
    }

    SymmetryDiagnostic diagnostic;
    diagnostic.symmetry_error = (max_abs > 0.0) ? max_diff / max_abs : 0.0;
    diagnostic.numerically_symmetric = diagnostic.symmetry_error <= tolerance;
    return diagnostic;
}

void CCSRMatrix::CheckFinite(const string& context) const
{
    require_finite_vector(values_, context);
}

double CCSRMatrix::RelativeResidual(const double* rhs, const double* solution) const
{
    const unsigned int N = dim();
    vector<double> x(N), b(N), Ax(N);
    for (unsigned int i = 0; i < N; ++i)
    {
        x[i] = solution[i];
        b[i] = rhs[i];
    }
    require_finite_vector(b, "checked residual RHS");
    require_finite_vector(x, "checked residual solution");
    MatVec(x, Ax);
    double norm_r = 0.0;
    double norm_b = 0.0;
    for (unsigned int i = 0; i < N; ++i)
    {
        const double diff = Ax[i] - b[i];
        norm_r += diff * diff;
        norm_b += b[i] * b[i];
    }
    require_finite_scalar(norm_r, "checked residual norm");
    return sqrt(norm_r) / (sqrt(norm_b) > 0.0 ? sqrt(norm_b) : 1.0);
}

CDiagonalScaling::CDiagonalScaling()
{
    Disable();
}

void CDiagonalScaling::Disable()
{
    enabled_ = false;
    scales_.clear();
    stats_.min_abs_diag = 0.0;
    stats_.max_abs_diag = 0.0;
    stats_.max_scale = 1.0;
    stats_.min_scale = 1.0;
    stats_.number_of_small_or_zero_diagonal_entries = 0;
}

void CDiagonalScaling::Setup(const CCSRMatrix& matrix, double epsD)
{
    enabled_ = true;
    const vector<double> diag = matrix.Diagonal();
    require_finite_vector(diag, "diagonal scaling diagonal");
    scales_.assign(diag.size(), 1.0);
    stats_.min_abs_diag = numeric_limits<double>::max();
    stats_.max_abs_diag = 0.0;
    stats_.max_scale = 0.0;
    stats_.min_scale = numeric_limits<double>::max();
    stats_.number_of_small_or_zero_diagonal_entries = 0;

    for (unsigned int i = 0; i < diag.size(); ++i)
    {
        const double abs_diag = fabs(diag[i]);
        stats_.min_abs_diag = min(stats_.min_abs_diag, abs_diag);
        stats_.max_abs_diag = max(stats_.max_abs_diag, abs_diag);
        if (abs_diag <= epsD)
            stats_.number_of_small_or_zero_diagonal_entries++;
        const double denom = max(abs_diag, epsD);
        scales_[i] = 1.0 / sqrt(denom);
        stats_.max_scale = max(stats_.max_scale, scales_[i]);
        stats_.min_scale = min(stats_.min_scale, scales_[i]);
    }

    if (diag.empty())
    {
        stats_.min_abs_diag = 0.0;
        stats_.min_scale = 1.0;
        stats_.max_scale = 1.0;
    }
}

void CDiagonalScaling::ScaleRightHandSide(const double* rhs, vector<double>& scaled_rhs) const
{
    if (!enabled_)
    {
        scaled_rhs.assign(rhs, rhs + scaled_rhs.size());
        return;
    }
    scaled_rhs.assign(scales_.size(), 0.0);
    for (unsigned int i = 0; i < scales_.size(); ++i)
        scaled_rhs[i] = scales_[i] * rhs[i];
}

void CDiagonalScaling::RecoverSolution(const vector<double>& scaled_solution, double* solution) const
{
    if (!enabled_)
    {
        for (unsigned int i = 0; i < scaled_solution.size(); ++i)
            solution[i] = scaled_solution[i];
        return;
    }
    for (unsigned int i = 0; i < scaled_solution.size(); ++i)
        solution[i] = scales_[i] * scaled_solution[i];
}

CScaledOperator::CScaledOperator(const CCSRMatrix& K, const CDiagonalScaling* scaling)
    : K_(K), Scaling_(scaling)
{
}

unsigned int CScaledOperator::dim() const
{
    return K_.dim();
}

void CScaledOperator::MatVec(const vector<double>& x, vector<double>& y) const
{
    if (!Scaling_ || !Scaling_->enabled())
    {
        K_.MatVec(x, y);
        return;
    }

    vector<double> tmp(x.size(), 0.0);
    for (unsigned int i = 0; i < x.size(); ++i)
        tmp[i] = Scaling_->Scale(i) * x[i];
    K_.MatVec(tmp, y);
    for (unsigned int i = 0; i < y.size(); ++i)
        y[i] *= Scaling_->Scale(i);
}

void CIdentityPreconditioner::Apply(const vector<double>& r, vector<double>& z) const
{
    z = r;
}

const char* CIdentityPreconditioner::Name() const
{
    return "none";
}

CJacobiPreconditioner::CJacobiPreconditioner(const CCSRMatrix& matrix, const CDiagonalScaling* scaling)
{
    vector<double> diag = matrix.Diagonal();
    require_finite_vector(diag, "Jacobi preconditioner diagonal");
    inv_diag_.assign(diag.size(), 1.0);
    for (unsigned int i = 0; i < diag.size(); ++i)
    {
        double value = diag[i];
        if (scaling && scaling->enabled())
        {
            const double scale = scaling->Scale(i);
            value *= scale * scale;
        }
        const double denom = max(fabs(value), 1.0e-300);
        inv_diag_[i] = 1.0 / denom;
    }
}

void CJacobiPreconditioner::Apply(const vector<double>& r, vector<double>& z) const
{
    require_finite_vector(r, "Jacobi preconditioner input");
    z.assign(inv_diag_.size(), 0.0);
    for (unsigned int i = 0; i < inv_diag_.size(); ++i)
        z[i] = r[i] * inv_diag_[i];
    require_finite_vector(z, "Jacobi preconditioner output");
}

const char* CJacobiPreconditioner::Name() const
{
    return "jacobi";
}

CSSORPreconditioner::CSSORPreconditioner(const CCSRMatrix& matrix, const CDiagonalScaling* scaling)
    : matrix_(matrix), scaling_(scaling), diag_(matrix.Diagonal())
{
    require_finite_vector(diag_, "SSOR preconditioner diagonal");
}

void CSSORPreconditioner::Apply(const vector<double>& r, vector<double>& z) const
{
    require_finite_vector(r, "SSOR preconditioner input");
    vector<double> y(matrix_.dim(), 0.0);
    z.assign(matrix_.dim(), 0.0);

    for (unsigned int i = 0; i < matrix_.dim(); ++i)
    {
        double sum = r[i];
        for (size_t k = matrix_.RowBegin(i); k < matrix_.RowEnd(i); ++k)
        {
            unsigned int col = matrix_.ColumnIndex(k);
            double value = matrix_.ValueAt(k);
            if (scaling_ && scaling_->enabled())
                value *= scaling_->Scale(i) * scaling_->Scale(col);
            if (col < i)
                sum -= value * y[col];
        }
        double diag = diag_[i];
        if (scaling_ && scaling_->enabled())
        {
            const double scale = scaling_->Scale(i);
            diag *= scale * scale;
        }
        y[i] = (fabs(diag) <= DBL_MIN) ? sum : sum / diag;
    }

    for (unsigned int i = 0; i < matrix_.dim(); ++i)
    {
        double diag = diag_[i];
        if (scaling_ && scaling_->enabled())
        {
            const double scale = scaling_->Scale(i);
            diag *= scale * scale;
        }
        y[i] *= diag;
    }

    for (int ii = static_cast<int>(matrix_.dim()) - 1; ii >= 0; --ii)
    {
        unsigned int i = static_cast<unsigned int>(ii);
        double sum = y[i];
        for (size_t k = matrix_.RowBegin(i); k < matrix_.RowEnd(i); ++k)
        {
            unsigned int col = matrix_.ColumnIndex(k);
            double value = matrix_.ValueAt(k);
            if (scaling_ && scaling_->enabled())
                value *= scaling_->Scale(i) * scaling_->Scale(col);
            if (col > i)
                sum -= value * z[col];
        }
        double diag = diag_[i];
        if (scaling_ && scaling_->enabled())
        {
            const double scale = scaling_->Scale(i);
            diag *= scale * scale;
        }
        z[i] = (fabs(diag) <= DBL_MIN) ? sum : sum / diag;
    }
    require_finite_vector(z, "SSOR preconditioner output");
}

const char* CSSORPreconditioner::Name() const
{
    return "ssor";
}

CBlockJacobiPreconditioner::CBlockJacobiPreconditioner(
    const CCSRMatrix& matrix, const CDiagonalScaling* scaling,
    const vector<vector<unsigned int> >& equation_blocks)
{
    stats_.num_blocks = 0;
    stats_.num_cholesky_success = 0;
    stats_.num_shifted_blocks = 0;
    stats_.num_diag_fallback_blocks = 0;
    for (unsigned int i = 0; i < 6; ++i)
        stats_.block_size_histogram[i] = 0;

    const vector<double> diag = matrix.Diagonal();
    require_finite_vector(diag, "block Jacobi diagonal");

    for (size_t block_index = 0; block_index < equation_blocks.size(); ++block_index)
    {
        const vector<unsigned int>& eqs = equation_blocks[block_index];
        if (eqs.empty())
            continue;

        Block block;
        block.equations = eqs;
        block.size = static_cast<unsigned int>(eqs.size());
        block.diagonal_fallback = false;
        block.factor.assign(block.size * block.size, 0.0);
        block.inv_diag.assign(block.size, 0.0);

        stats_.num_blocks++;
        stats_.block_size_histogram[min(block.size, 6u) - 1]++;

        double max_abs_diag = 0.0;
        for (unsigned int i = 0; i < block.size; ++i)
        {
            const unsigned int row = eqs[i];
            for (unsigned int j = 0; j < block.size; ++j)
            {
                const unsigned int col = eqs[j];
                double value = matrix.GetValue(row, col);
                if (scaling && scaling->enabled())
                    value *= scaling->Scale(row) * scaling->Scale(col);
                block.factor[i * block.size + j] = value;
            }
            max_abs_diag = max(max_abs_diag, fabs(block.factor[i * block.size + i]));
        }
        if (max_abs_diag <= 0.0)
            max_abs_diag = 1.0;

        for (unsigned int i = 0; i < block.size; ++i)
        {
            for (unsigned int j = i + 1; j < block.size; ++j)
            {
                const double sym = 0.5 * (block.factor[i * block.size + j] +
                                          block.factor[j * block.size + i]);
                block.factor[i * block.size + j] = sym;
                block.factor[j * block.size + i] = sym;
            }
        }

        vector<double> factor = block.factor;
        bool success = dense_cholesky_factor(factor, block.size);
        bool shifted = false;
        double shift = 1.0e-12 * max_abs_diag;
        while (!success && shift <= 1.0e-6 * max_abs_diag)
        {
            factor = block.factor;
            for (unsigned int i = 0; i < block.size; ++i)
                factor[i * block.size + i] += shift;
            success = dense_cholesky_factor(factor, block.size);
            if (success)
                shifted = true;
            else
                shift *= 10.0;
        }

        if (success)
        {
            block.factor.swap(factor);
            stats_.num_cholesky_success++;
            if (shifted)
                stats_.num_shifted_blocks++;
        }
        else
        {
            block.diagonal_fallback = true;
            for (unsigned int i = 0; i < block.size; ++i)
            {
                double value = diag[eqs[i]];
                if (scaling && scaling->enabled())
                {
                    const double scale = scaling->Scale(eqs[i]);
                    value *= scale * scale;
                }
                block.inv_diag[i] = 1.0 / max(fabs(value), 1.0e-300);
            }
            stats_.num_diag_fallback_blocks++;
        }

        blocks_.push_back(block);
    }
}

void CBlockJacobiPreconditioner::Apply(const vector<double>& r, vector<double>& z) const
{
    require_finite_vector(r, "block Jacobi input");
    z.assign(r.size(), 0.0);
    vector<double> local_rhs;
    vector<double> local_x;
    for (size_t block_index = 0; block_index < blocks_.size(); ++block_index)
    {
        const Block& block = blocks_[block_index];
        local_rhs.assign(block.size, 0.0);
        for (unsigned int i = 0; i < block.size; ++i)
            local_rhs[i] = r[block.equations[i]];

        if (block.diagonal_fallback)
        {
            for (unsigned int i = 0; i < block.size; ++i)
                z[block.equations[i]] = local_rhs[i] * block.inv_diag[i];
        }
        else
        {
            dense_cholesky_solve(block.factor, block.size, local_rhs, local_x);
            for (unsigned int i = 0; i < block.size; ++i)
                z[block.equations[i]] = local_x[i];
        }
    }
    require_finite_vector(z, "block Jacobi output");
}

const char* CBlockJacobiPreconditioner::Name() const
{
    return "block-jacobi";
}

CILU0Preconditioner::CILU0Preconditioner(const CCSRMatrix& matrix, const CDiagonalScaling* scaling)
    : n_(matrix.dim()), row_ptr_(n_ + 1, 0)
{
    stats_.shifted_pivots = 0;
    stats_.zero_pattern_fallbacks = 0;

    for (unsigned int row = 0; row < n_; ++row)
        row_ptr_[row + 1] = row_ptr_[row] + (matrix.RowEnd(row) - matrix.RowBegin(row));
    col_ind_.assign(row_ptr_[n_], 0);
    lu_.assign(row_ptr_[n_], 0.0);

    for (unsigned int row = 0; row < n_; ++row)
    {
        size_t out = row_ptr_[row];
        for (size_t k = matrix.RowBegin(row); k < matrix.RowEnd(row); ++k, ++out)
        {
            const unsigned int col = matrix.ColumnIndex(k);
            col_ind_[out] = col;
            double value = matrix.ValueAt(k);
            if (scaling && scaling->enabled())
                value *= scaling->Scale(row) * scaling->Scale(col);
            lu_[out] = value;
        }
    }

    vector<size_t> diag_index(n_, static_cast<size_t>(-1));
    for (unsigned int row = 0; row < n_; ++row)
    {
        for (size_t k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
        {
            if (col_ind_[k] == row)
            {
                diag_index[row] = k;
                break;
            }
        }
        if (diag_index[row] == static_cast<size_t>(-1))
            throw runtime_error("*** Error *** ILU(0) requires an explicit diagonal entry in every row.");
    }

    for (unsigned int i = 0; i < n_; ++i)
    {
        for (size_t kk = row_ptr_[i]; kk < row_ptr_[i + 1]; ++kk)
        {
            const unsigned int kcol = col_ind_[kk];
            if (kcol >= i)
                break;

            const double pivot = lu_[diag_index[kcol]];
            if (!isfinite(pivot))
                throw runtime_error("*** Error *** Non-finite ILU(0) pivot.");
            if (fabs(pivot) <= 1.0e-30)
            {
                lu_[diag_index[kcol]] = (pivot >= 0.0 ? 1.0 : -1.0) * 1.0e-12;
                stats_.shifted_pivots++;
            }

            lu_[kk] /= lu_[diag_index[kcol]];

            size_t ij = kk + 1;
            size_t kj = diag_index[kcol] + 1;
            while (ij < row_ptr_[i + 1] && kj < row_ptr_[kcol + 1])
            {
                const unsigned int col_i = col_ind_[ij];
                const unsigned int col_k = col_ind_[kj];
                if (col_i == col_k)
                {
                    lu_[ij] -= lu_[kk] * lu_[kj];
                    ++ij;
                    ++kj;
                }
                else if (col_i < col_k)
                {
                    ++ij;
                }
                else
                {
                    ++kj;
                }
            }
        }

        double& diag = lu_[diag_index[i]];
        if (!isfinite(diag))
            throw runtime_error("*** Error *** Non-finite ILU(0) diagonal.");
        if (fabs(diag) <= 1.0e-30)
        {
            diag = (diag >= 0.0 ? 1.0 : -1.0) * 1.0e-12;
            stats_.shifted_pivots++;
        }
    }
}

void CILU0Preconditioner::Apply(const vector<double>& r, vector<double>& z) const
{
    require_finite_vector(r, "ILU0 preconditioner input");
    z.assign(n_, 0.0);
    vector<double> y(n_, 0.0);

    for (unsigned int i = 0; i < n_; ++i)
    {
        double sum = r[i];
        double diag = 1.0;
        for (size_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k)
        {
            const unsigned int col = col_ind_[k];
            if (col < i)
                sum -= lu_[k] * y[col];
            else if (col == i)
            {
                diag = lu_[k];
                break;
            }
        }
        y[i] = sum;
        require_finite_scalar(diag, "ILU0 lower diagonal");
    }

    for (int ii = static_cast<int>(n_) - 1; ii >= 0; --ii)
    {
        const unsigned int i = static_cast<unsigned int>(ii);
        double sum = y[i];
        double diag = 0.0;
        for (size_t k = row_ptr_[i]; k < row_ptr_[i + 1]; ++k)
        {
            const unsigned int col = col_ind_[k];
            if (col > i)
                sum -= lu_[k] * z[col];
            else if (col == i)
                diag = lu_[k];
        }
        if (fabs(diag) <= 1.0e-30)
            throw runtime_error("*** Error *** ILU(0) encountered a zero pivot during apply.");
        z[i] = sum / diag;
    }

    require_finite_vector(z, "ILU0 preconditioner output");
}

const char* CILU0Preconditioner::Name() const
{
    return "ilu0";
}

CPCGSolver::CPCGSolver(const CLinearOperator& K, const CPreconditioner& preconditioner,
                       double tolerance, unsigned int max_iterations, unsigned int log_interval)
    : K_(K), Preconditioner_(preconditioner), tolerance_(tolerance),
      max_iterations_(max_iterations), log_interval_(log_interval)
{
}

CPCGSolver::Result CPCGSolver::Solve(const double* rhs, double* solution) const
{
    const unsigned int N = K_.dim();
    vector<double> b(N), x(N, 0.0), r(N), z(N), p(N), Ap(N);
    for (unsigned int i = 0; i < N; ++i)
        b[i] = rhs[i];
    require_finite_vector(b, "PCG RHS");

    const double norm_b = sqrt(dot_product(b, b));
    const double denom = norm_b > 0.0 ? norm_b : 1.0;
    r = b;

    Preconditioner_.Apply(r, z);
    p = z;

    double rz_old = dot_product(r, z);
    Result result;
    result.converged = false;
    result.iterations = 0;
    result.initial_relative_residual = sqrt(dot_product(r, r)) / denom;
    result.final_relative_residual = result.initial_relative_residual;
    solver_output() << "  PCG iter 0 relative_residual=" << result.initial_relative_residual << endl;

    if (result.initial_relative_residual <= tolerance_)
    {
        for (unsigned int i = 0; i < N; ++i)
            solution[i] = 0.0;
        result.converged = true;
        return result;
    }

    vector<double> residual_window;
    residual_window.push_back(result.initial_relative_residual);
    for (unsigned int iter = 1; iter <= max_iterations_; ++iter)
    {
        K_.MatVec(p, Ap);
        double denom_alpha = dot_product(p, Ap);
        if (!isfinite(denom_alpha) || denom_alpha <= DBL_MIN)
        {
            solver_output() << "  PCG stopped: non-positive or invalid pAp=" << denom_alpha
                            << " at iteration " << iter
                            << ". Suggest using sparse-bicgstab." << endl;
            break;
        }
        double alpha = rz_old / denom_alpha;
        if (!isfinite(alpha))
        {
            solver_output() << "  PCG stopped: invalid alpha at iteration " << iter
                            << ". Suggest using sparse-bicgstab." << endl;
            break;
        }

        for (unsigned int i = 0; i < N; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double rr = dot_product(r, r);
        double rel = sqrt(rr) / denom;
        result.iterations = iter;
        result.final_relative_residual = rel;
        if (log_interval_ && (iter == 1 || iter % log_interval_ == 0))
            solver_output() << "  PCG iter " << iter << " relative_residual=" << rel << endl;
        print_stagnation_warning("PCG", iter, residual_window, rel);
        if (rel <= tolerance_)
        {
            result.converged = true;
            break;
        }

        Preconditioner_.Apply(r, z);
        double rz_new = dot_product(r, z);
        if (!isfinite(rz_new) || fabs(rz_old) <= DBL_MIN)
        {
            solver_output() << "  PCG stopped: invalid preconditioned residual at iteration "
                            << iter << ". Suggest using sparse-bicgstab." << endl;
            break;
        }
        double beta = rz_new / rz_old;
        if (!isfinite(beta))
        {
            solver_output() << "  PCG stopped: invalid beta at iteration " << iter
                            << ". Suggest using sparse-bicgstab." << endl;
            break;
        }
        for (unsigned int i = 0; i < N; ++i)
            p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
    }

    for (unsigned int i = 0; i < N; ++i)
        solution[i] = x[i];
    return result;
}

CBiCGSTABSolver::CBiCGSTABSolver(const CLinearOperator& K, const CPreconditioner& preconditioner,
                                 double tolerance, unsigned int max_iterations,
                                 unsigned int log_interval)
    : K_(K), Preconditioner_(preconditioner), tolerance_(tolerance),
      max_iterations_(max_iterations), log_interval_(log_interval)
{
}

CBiCGSTABSolver::Result CBiCGSTABSolver::Solve(const double* rhs, double* solution) const
{
    const unsigned int N = K_.dim();
    vector<double> b(N), x(N, 0.0), r(N), r_hat(N), p(N, 0.0), v(N, 0.0);
    vector<double> s(N), t(N), phat(N), shat(N);
    for (unsigned int i = 0; i < N; ++i)
        b[i] = rhs[i];
    require_finite_vector(b, "BiCGSTAB RHS");

    const double norm_b = sqrt(dot_product(b, b));
    const double denom = norm_b > 0.0 ? norm_b : 1.0;
    r = b;
    r_hat = r;

    Result result;
    result.converged = false;
    result.iterations = 0;
    result.initial_relative_residual = sqrt(dot_product(r, r)) / denom;
    result.final_relative_residual = result.initial_relative_residual;
    solver_output() << "  BiCGSTAB iter 0 relative_residual=" << result.initial_relative_residual << endl;
    if (result.initial_relative_residual <= tolerance_)
    {
        for (unsigned int i = 0; i < N; ++i)
            solution[i] = 0.0;
        result.converged = true;
        return result;
    }

    vector<double> residual_window;
    residual_window.push_back(result.initial_relative_residual);
    double rho_old = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    for (unsigned int iter = 1; iter <= max_iterations_; ++iter)
    {
        double rho_new = dot_product(r_hat, r);
        if (!isfinite(rho_new) || fabs(rho_new) <= DBL_MIN)
        {
            solver_output() << "  BiCGSTAB stopped: invalid rho at iteration " << iter
                            << ". Suggest using sparse-gmres for diagnostics." << endl;
            break;
        }
        double beta = (rho_new / rho_old) * (alpha / omega);
        for (unsigned int i = 0; i < N; ++i)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);

        Preconditioner_.Apply(p, phat);
        K_.MatVec(phat, v);
        double rhat_v = dot_product(r_hat, v);
        if (!isfinite(rhat_v) || fabs(rhat_v) <= DBL_MIN)
        {
            solver_output() << "  BiCGSTAB stopped: invalid alpha denominator at iteration "
                            << iter << ". Suggest using sparse-gmres for diagnostics." << endl;
            break;
        }
        alpha = rho_new / rhat_v;
        if (!isfinite(alpha))
        {
            solver_output() << "  BiCGSTAB stopped: invalid alpha at iteration " << iter
                            << ". Suggest using sparse-gmres for diagnostics." << endl;
            break;
        }

        for (unsigned int i = 0; i < N; ++i)
            s[i] = r[i] - alpha * v[i];
        double rel_s = sqrt(dot_product(s, s)) / denom;
        if (rel_s <= tolerance_)
        {
            for (unsigned int i = 0; i < N; ++i)
                x[i] += alpha * phat[i];
            result.converged = true;
            result.iterations = iter;
            result.final_relative_residual = rel_s;
            print_stagnation_warning("BiCGSTAB", iter, residual_window, rel_s);
            break;
        }

        Preconditioner_.Apply(s, shat);
        K_.MatVec(shat, t);
        double tt = dot_product(t, t);
        if (!isfinite(tt) || tt <= DBL_MIN)
        {
            solver_output() << "  BiCGSTAB stopped: invalid omega denominator at iteration "
                            << iter << ". Suggest using sparse-gmres for diagnostics." << endl;
            break;
        }
        omega = dot_product(t, s) / tt;
        if (!isfinite(omega) || fabs(omega) <= DBL_MIN)
        {
            solver_output() << "  BiCGSTAB stopped: invalid omega at iteration " << iter
                            << ". Suggest using sparse-gmres for diagnostics." << endl;
            break;
        }

        for (unsigned int i = 0; i < N; ++i)
        {
            x[i] += alpha * phat[i] + omega * shat[i];
            r[i] = s[i] - omega * t[i];
        }

        double rr = dot_product(r, r);
        double rel = sqrt(rr) / denom;
        result.iterations = iter;
        result.final_relative_residual = rel;
        if (log_interval_ && (iter == 1 || iter % log_interval_ == 0))
            solver_output() << "  BiCGSTAB iter " << iter << " relative_residual=" << rel << endl;
        print_stagnation_warning("BiCGSTAB", iter, residual_window, rel);
        if (rel <= tolerance_)
        {
            result.converged = true;
            break;
        }

        rho_old = rho_new;
    }

    for (unsigned int i = 0; i < N; ++i)
        solution[i] = x[i];
    return result;
}

CGMRESSolver::CGMRESSolver(const CLinearOperator& K, const CPreconditioner& preconditioner,
                           double tolerance, unsigned int max_iterations,
                           unsigned int restart, unsigned int log_interval)
    : K_(K), Preconditioner_(preconditioner), tolerance_(tolerance),
      max_iterations_(max_iterations), restart_(restart ? restart : 50),
      log_interval_(log_interval)
{
}

CGMRESSolver::Result CGMRESSolver::Solve(const double* rhs, double* solution) const
{
    const unsigned int N = K_.dim();
    const unsigned int m = restart_;
    vector<double> b(N), x(N, 0.0), Ax(N), r(N), z(N), w(N), y(m + 1, 0.0);
    for (unsigned int i = 0; i < N; ++i)
        b[i] = rhs[i];
    require_finite_vector(b, "GMRES RHS");

    const double norm_b = sqrt(dot_product(b, b));
    const double denom = norm_b > 0.0 ? norm_b : 1.0;

    Result result;
    result.converged = false;
    result.iterations = 0;
    result.initial_relative_residual = norm_b / denom;
    result.final_relative_residual = result.initial_relative_residual;
    solver_output() << "  GMRES iter 0 relative_residual=" << result.initial_relative_residual << endl;
    if (result.initial_relative_residual <= tolerance_)
    {
        for (unsigned int i = 0; i < N; ++i)
            solution[i] = 0.0;
        result.converged = true;
        return result;
    }

    vector<vector<double> > V(m + 1, vector<double>(N, 0.0));
    vector<vector<double> > Z(m, vector<double>(N, 0.0));
    vector<vector<double> > H(m + 1, vector<double>(m, 0.0));
    vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0);
    vector<double> residual_window;
    residual_window.push_back(result.initial_relative_residual);

    while (result.iterations < max_iterations_ && !result.converged)
    {
        K_.MatVec(x, Ax);
        double rr = 0.0;
        for (unsigned int i = 0; i < N; ++i)
        {
            r[i] = b[i] - Ax[i];
            rr += r[i] * r[i];
        }
        double beta = sqrt(rr);
        result.final_relative_residual = beta / denom;
        require_finite_scalar(result.final_relative_residual, "GMRES restart residual");
        if (result.final_relative_residual <= tolerance_)
        {
            result.converged = true;
            break;
        }

        for (unsigned int i = 0; i < N; ++i)
            V[0][i] = r[i] / beta;
        fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        unsigned int inner = 0;
        for (; inner < m && result.iterations < max_iterations_; ++inner)
        {
            Preconditioner_.Apply(V[inner], Z[inner]);
            K_.MatVec(Z[inner], w);
            for (unsigned int k = 0; k <= inner; ++k)
            {
                H[k][inner] = dot_product(w, V[k]);
                for (unsigned int i = 0; i < N; ++i)
                    w[i] -= H[k][inner] * V[k][i];
            }
            H[inner + 1][inner] = sqrt(dot_product(w, w));
            if (H[inner + 1][inner] > DBL_MIN)
            {
                for (unsigned int i = 0; i < N; ++i)
                    V[inner + 1][i] = w[i] / H[inner + 1][inner];
            }

            for (unsigned int k = 0; k < inner; ++k)
            {
                double temp = cs[k] * H[k][inner] + sn[k] * H[k + 1][inner];
                H[k + 1][inner] = -sn[k] * H[k][inner] + cs[k] * H[k + 1][inner];
                H[k][inner] = temp;
            }

            double h0 = H[inner][inner];
            double h1 = H[inner + 1][inner];
            double norm = sqrt(h0 * h0 + h1 * h1);
            if (!isfinite(norm) || norm <= DBL_MIN)
            {
                solver_output() << "  GMRES stopped: invalid Hessenberg norm at iteration "
                                << result.iterations + 1 << endl;
                break;
            }
            cs[inner] = h0 / norm;
            sn[inner] = h1 / norm;
            H[inner][inner] = cs[inner] * h0 + sn[inner] * h1;
            H[inner + 1][inner] = 0.0;

            double gtmp = cs[inner] * g[inner] + sn[inner] * g[inner + 1];
            g[inner + 1] = -sn[inner] * g[inner] + cs[inner] * g[inner + 1];
            g[inner] = gtmp;

            result.iterations++;
            result.final_relative_residual = fabs(g[inner + 1]) / denom;
            require_finite_scalar(result.final_relative_residual, "GMRES residual");
            if (log_interval_ && (result.iterations == 1 || result.iterations % log_interval_ == 0))
                solver_output() << "  GMRES iter " << result.iterations
                                << " relative_residual=" << result.final_relative_residual << endl;
            print_stagnation_warning("GMRES", result.iterations, residual_window,
                                     result.final_relative_residual);
            if (result.final_relative_residual <= tolerance_)
            {
                inner++;
                result.converged = true;
                break;
            }
        }

        if (inner == 0)
            break;

        y.assign(m + 1, 0.0);
        for (int row = static_cast<int>(inner) - 1; row >= 0; --row)
        {
            double sum = g[static_cast<unsigned int>(row)];
            for (unsigned int col = static_cast<unsigned int>(row) + 1; col < inner; ++col)
                sum -= H[static_cast<unsigned int>(row)][col] * y[col];
            if (fabs(H[static_cast<unsigned int>(row)][static_cast<unsigned int>(row)]) <= DBL_MIN)
            {
                solver_output() << "  GMRES stopped: singular least-squares system." << endl;
                inner = 0;
                break;
            }
            y[static_cast<unsigned int>(row)] =
                sum / H[static_cast<unsigned int>(row)][static_cast<unsigned int>(row)];
        }
        if (inner == 0)
            break;

        for (unsigned int col = 0; col < inner; ++col)
        {
            for (unsigned int i = 0; i < N; ++i)
                x[i] += Z[col][i] * y[col];
        }
    }

    for (unsigned int i = 0; i < N; ++i)
        solution[i] = x[i];
    return result;
}

CSparseSolverOptions::CSparseSolverOptions()
    : backend_name("standard"),
      requested_solver("sparse-auto"),
      requested_preconditioner("ssor"),
      pardiso_mtype_mode("auto"),
      scale_mode("diag"),
      tolerance(1.0e-6),
      max_iterations(5000),
      residual_print_interval(50),
      gmres_restart(80)
{
}

CSparseBackendSetupTimings::CSparseBackendSetupTimings()
    : symmetry_check_time(0.0),
      scaling_time(0.0),
      preconditioner_setup_time(0.0)
{
}

CSparseBackendInfo::CSparseBackendInfo()
    : backend_available(true),
      backend_name("standard"),
      requested_solver("sparse-auto"),
      actual_solver("sparse-bicgstab"),
      requested_preconditioner("ssor"),
      actual_preconditioner("block-jacobi"),
      scale_mode("diag"),
      scaling_enabled(false),
      has_block_jacobi_stats(false),
      has_ilu0_stats(false)
{
    symmetry.symmetry_error = 0.0;
    symmetry.numerically_symmetric = false;
    scaling_stats.min_abs_diag = 0.0;
    scaling_stats.max_abs_diag = 0.0;
    scaling_stats.max_scale = 1.0;
    scaling_stats.min_scale = 1.0;
    scaling_stats.number_of_small_or_zero_diagonal_entries = 0;
    for (unsigned int i = 0; i < 6; ++i)
        block_jacobi_stats.block_size_histogram[i] = 0;
    block_jacobi_stats.num_blocks = 0;
    block_jacobi_stats.num_cholesky_success = 0;
    block_jacobi_stats.num_shifted_blocks = 0;
    block_jacobi_stats.num_diag_fallback_blocks = 0;
    ilu0_stats.shifted_pivots = 0;
    ilu0_stats.zero_pattern_fallbacks = 0;
    pardiso_info.enabled = false;
    pardiso_info.requested_mtype_mode = "auto";
    pardiso_info.selected_mtype_mode = "";
    pardiso_info.selected_mtype = 0;
    pardiso_info.matrix_part = "";
    pardiso_info.attempt_count = 0;
    pardiso_info.retried_from_spd_to_sym_indef = false;
    pardiso_info.phase11_error = 0;
    pardiso_info.phase22_error = 0;
    pardiso_info.phase33_error = 0;
    pardiso_info.factor_nnz = -1;
    pardiso_info.peak_memory_kb = -1;
}

CSparseSolveResult::CSparseSolveResult()
    : converged(false),
      iterations(0),
      initial_relative_residual(0.0),
      final_relative_residual(0.0)
{
}

CStandardSparseSolverBackend::CStandardSparseSolverBackend()
    : Matrix_(0),
      EquationBlocks_(0)
{
}

void CStandardSparseSolverBackend::Setup(const CCSRMatrix& matrix,
                                         const vector<vector<unsigned int> >* equation_blocks,
                                         const CSparseSolverOptions& options)
{
    Matrix_ = &matrix;
    EquationBlocks_ = equation_blocks;
    Options_ = options;
    Info_ = CSparseBackendInfo();
    Info_.backend_name = "standard";
    Info_.requested_solver = options.requested_solver;
    Info_.requested_preconditioner = options.requested_preconditioner;
    Info_.scale_mode = options.scale_mode;

    Clock symmetryTimer;
    symmetryTimer.Start();
    Info_.symmetry = matrix.CheckSymmetry();
    Info_.timings.symmetry_check_time = symmetryTimer.ElapsedTime();

    Info_.actual_solver = options.requested_solver;
    if (options.requested_solver == "sparse-auto")
        Info_.actual_solver = "sparse-bicgstab";

    Info_.actual_preconditioner = options.requested_preconditioner;
    if (options.requested_solver == "sparse-auto" && options.requested_preconditioner == "ssor")
        Info_.actual_preconditioner = "block-jacobi";

    Clock scalingTimer;
    scalingTimer.Start();
    if (options.scale_mode == "diag")
        Scaling_.Setup(matrix);
    else
        Scaling_.Disable();
    Info_.timings.scaling_time = scalingTimer.ElapsedTime();
    Info_.scaling_enabled = Scaling_.enabled();
    Info_.scaling_stats = Scaling_.GetStats();
    ScaledOperator_.reset(new CScaledOperator(matrix, Scaling_.enabled() ? &Scaling_ : 0));

    Clock preconditionerTimer;
    preconditionerTimer.Start();
    if (Info_.actual_preconditioner == "none")
    {
        Preconditioner_.reset(new CIdentityPreconditioner());
    }
    else if (Info_.actual_preconditioner == "jacobi")
    {
        Preconditioner_.reset(new CJacobiPreconditioner(matrix, Scaling_.enabled() ? &Scaling_ : 0));
    }
    else if (Info_.actual_preconditioner == "ssor")
    {
        Preconditioner_.reset(new CSSORPreconditioner(matrix, Scaling_.enabled() ? &Scaling_ : 0));
    }
    else if (Info_.actual_preconditioner == "ilu0")
    {
        Preconditioner_.reset(new CILU0Preconditioner(matrix, Scaling_.enabled() ? &Scaling_ : 0));
        const CILU0Preconditioner* ilu =
            dynamic_cast<const CILU0Preconditioner*>(Preconditioner_.get());
        if (ilu)
        {
            Info_.has_ilu0_stats = true;
            Info_.ilu0_stats = ilu->GetStats();
        }
    }
    else
    {
        if (!EquationBlocks_)
            throw runtime_error("*** Error *** Block Jacobi preconditioner requires equation blocks.");
        Preconditioner_.reset(
            new CBlockJacobiPreconditioner(matrix, Scaling_.enabled() ? &Scaling_ : 0, *EquationBlocks_));
        const CBlockJacobiPreconditioner* block =
            dynamic_cast<const CBlockJacobiPreconditioner*>(Preconditioner_.get());
        if (block)
        {
            Info_.has_block_jacobi_stats = true;
            Info_.block_jacobi_stats = block->GetStats();
        }
    }
    Info_.timings.preconditioner_setup_time = preconditionerTimer.ElapsedTime();
}

const CSparseBackendInfo& CStandardSparseSolverBackend::GetInfo() const
{
    return Info_;
}

CSparseSolveResult CStandardSparseSolverBackend::Solve(const double* rhs, double* solution) const
{
    if (!Matrix_ || !Preconditioner_.get())
        throw runtime_error("*** Error *** Sparse solver backend is not initialized.");

    vector<double> scaled_rhs(Matrix_->dim(), 0.0);
    if (Scaling_.enabled())
        Scaling_.ScaleRightHandSide(rhs, scaled_rhs);
    else
    {
        for (unsigned int i = 0; i < Matrix_->dim(); ++i)
            scaled_rhs[i] = rhs[i];
    }

    vector<double> scaled_solution(Matrix_->dim(), 0.0);
    CPCGSolver::Result result;
    if (Info_.actual_solver == "sparse-gmres")
    {
        CGMRESSolver solver(*ScaledOperator_, *Preconditioner_, Options_.tolerance,
                            Options_.max_iterations, Options_.gmres_restart,
                            Options_.residual_print_interval);
        result = solver.Solve(scaled_rhs.data(), scaled_solution.data());
    }
    else if (Info_.actual_solver == "sparse-bicgstab")
    {
        CBiCGSTABSolver solver(*ScaledOperator_, *Preconditioner_, Options_.tolerance,
                               Options_.max_iterations, Options_.residual_print_interval);
        result = solver.Solve(scaled_rhs.data(), scaled_solution.data());
    }
    else
    {
        CPCGSolver solver(*ScaledOperator_, *Preconditioner_, Options_.tolerance,
                          Options_.max_iterations, Options_.residual_print_interval);
        result = solver.Solve(scaled_rhs.data(), scaled_solution.data());
    }

    if (Scaling_.enabled())
        Scaling_.RecoverSolution(scaled_solution, solution);
    else
    {
        for (unsigned int i = 0; i < Matrix_->dim(); ++i)
            solution[i] = scaled_solution[i];
    }

    CSparseSolveResult backend_result;
    backend_result.converged = result.converged;
    backend_result.iterations = result.iterations;
    backend_result.initial_relative_residual = result.initial_relative_residual;
    backend_result.final_relative_residual = result.final_relative_residual;
    return backend_result;
}

CPardisoSparseSolverBackend::CPardisoSparseSolverBackend()
    : Matrix_(0)
{
}

void CPardisoSparseSolverBackend::Setup(const CCSRMatrix& matrix,
                                        const vector<vector<unsigned int> >* equation_blocks,
                                        const CSparseSolverOptions& options)
{
    (void)equation_blocks;
    Matrix_ = &matrix;
    Options_ = options;
    Info_ = CSparseBackendInfo();
    Info_.backend_name = "pardiso";
    Info_.requested_solver = options.requested_solver;
    Info_.requested_preconditioner = options.requested_preconditioner;
    Info_.scale_mode = options.scale_mode;
    Info_.actual_solver = "pardiso";
    Info_.actual_preconditioner = "direct";
    Info_.scaling_enabled = false;
    Info_.pardiso_info.enabled = true;
    Info_.pardiso_info.requested_mtype_mode = options.pardiso_mtype_mode;

    Clock symmetryTimer;
    symmetryTimer.Start();
    Info_.symmetry = matrix.CheckSymmetry();
    Info_.timings.symmetry_check_time = symmetryTimer.ElapsedTime();

#ifdef STAPPP_ENABLE_MKL_PARDISO
    Info_.backend_available = true;
#else
    Info_.backend_available = false;
#endif
}

const CSparseBackendInfo& CPardisoSparseSolverBackend::GetInfo() const
{
    return Info_;
}

CSparseSolveResult CPardisoSparseSolverBackend::Solve(const double* rhs, double* solution) const
{
    if (!Matrix_)
        throw runtime_error("*** Error *** PARDISO backend is not initialized.");

#ifdef STAPPP_ENABLE_MKL_PARDISO
    Info_.pardiso_info.enabled = true;
    Info_.pardiso_info.requested_mtype_mode = Options_.pardiso_mtype_mode;
    Info_.pardiso_info.selected_mtype_mode.clear();
    Info_.pardiso_info.selected_mtype = 0;
    Info_.pardiso_info.matrix_part.clear();
    Info_.pardiso_info.attempt_count = 0;
    Info_.pardiso_info.retried_from_spd_to_sym_indef = false;
    Info_.pardiso_info.phase11_error = 0;
    Info_.pardiso_info.phase22_error = 0;
    Info_.pardiso_info.phase33_error = 0;
    Info_.pardiso_info.factor_nnz = -1;
    Info_.pardiso_info.peak_memory_kb = -1;

    struct AttemptConfig
    {
        const char* mode;
        int mtype;
        bool upper_only;
    };

    vector<AttemptConfig> attempts;
    if (Options_.pardiso_mtype_mode == "spd")
        attempts.push_back(AttemptConfig{"spd", 2, true});
    else if (Options_.pardiso_mtype_mode == "sym-indef")
        attempts.push_back(AttemptConfig{"sym-indef", -2, true});
    else if (Options_.pardiso_mtype_mode == "unsym")
        attempts.push_back(AttemptConfig{"unsym", 11, false});
    else
    {
        if (Info_.symmetry.symmetry_error <= 1.0e-10)
            attempts.push_back(AttemptConfig{"spd", 2, true});
        else
            attempts.push_back(AttemptConfig{"sym-indef", -2, true});
    }

    for (size_t attempt_index = 0; attempt_index < attempts.size(); ++attempt_index)
    {
        const AttemptConfig& config = attempts[attempt_index];
        ++Info_.pardiso_info.attempt_count;
        Info_.pardiso_info.selected_mtype_mode = config.mode;
        Info_.pardiso_info.selected_mtype = config.mtype;
        Info_.pardiso_info.matrix_part = config.upper_only ? "upper" : "full";

        PardisoAttemptResult attempt = run_pardiso_attempt(*Matrix_, rhs, config.mtype,
                                                           config.upper_only, solution);
        Info_.pardiso_info.phase11_error = attempt.phase11_error;
        Info_.pardiso_info.phase22_error = attempt.phase22_error;
        Info_.pardiso_info.phase33_error = attempt.phase33_error;
        Info_.pardiso_info.factor_nnz = attempt.factor_nnz;
        Info_.pardiso_info.peak_memory_kb = attempt.peak_memory_kb;

        if (attempt.phase11_error == 0 &&
            attempt.phase22_error == 0 &&
            attempt.phase33_error == 0)
        {
            CSparseSolveResult result;
            result.converged = true;
            result.iterations = 0;
            vector<double> zero_solution(Matrix_->dim(), 0.0);
            result.initial_relative_residual = Matrix_->RelativeResidual(rhs, zero_solution.data());
            result.final_relative_residual = Matrix_->RelativeResidual(rhs, solution);
            return result;
        }

        const bool should_retry_from_spd =
            Options_.pardiso_mtype_mode == "auto" &&
            string(config.mode) == "spd" &&
            attempt.phase11_error == 0 &&
            attempt.phase22_error == -4;
        if (should_retry_from_spd)
        {
            Info_.pardiso_info.retried_from_spd_to_sym_indef = true;
            attempts.push_back(AttemptConfig{"sym-indef", -2, true});
            continue;
        }

        ostringstream msg;
        if (attempt.phase11_error != 0)
            msg << "*** Error *** PARDISO symbolic factorization failed with error code " << attempt.phase11_error;
        else if (attempt.phase22_error != 0)
            msg << "*** Error *** PARDISO numeric factorization failed with error code " << attempt.phase22_error;
        else
            msg << "*** Error *** PARDISO solve phase failed with error code " << attempt.phase33_error;
        msg << " (mtype=" << config.mode
            << ", matrix_part=" << (config.upper_only ? "upper" : "full") << ")";
        throw runtime_error(msg.str());
    }

    throw runtime_error("*** Error *** PARDISO failed without producing a valid attempt result.");
#else
    (void)rhs;
    (void)solution;
    throw runtime_error("*** Error *** PARDISO backend requested, but STAPPP_ENABLE_MKL_PARDISO is OFF. Reconfigure with MKL support on a machine that has Intel oneMKL installed.");
#endif
}

unique_ptr<CSparseSolverBackend> CreateSparseSolverBackend(const string& backend_name)
{
    if (backend_name == "standard")
        return unique_ptr<CSparseSolverBackend>(new CStandardSparseSolverBackend());
    if (backend_name == "pardiso")
        return unique_ptr<CSparseSolverBackend>(new CPardisoSparseSolverBackend());

    throw runtime_error(string("*** Error *** Unsupported sparse solver backend: ") + backend_name);
}

bool IsSparseBackendAvailable(const string& backend_name)
{
    if (backend_name == "standard")
        return true;
    if (backend_name == "pardiso")
    {
#ifdef STAPPP_ENABLE_MKL_PARDISO
        return true;
#else
        return false;
#endif
    }
    return false;
}

CSparseSolveResult SolvePardisoSystem(const CSparseSymmetricMatrix& matrix,
                                      const double* rhs,
                                      double* solution,
                                      const string& requested_mtype_mode,
                                      CSparseBackendInfo& info)
{
    info = CSparseBackendInfo();
    info.backend_name = "pardiso";
    info.backend_available = IsSparseBackendAvailable("pardiso");
    info.requested_solver = "sparse-auto";
    info.actual_solver = "pardiso";
    info.requested_preconditioner = "direct";
    info.actual_preconditioner = "direct";
    info.scale_mode = "none";
    info.scaling_enabled = false;
    info.symmetry.symmetry_error = 0.0;
    info.symmetry.numerically_symmetric = true;
    info.pardiso_info.enabled = true;
    info.pardiso_info.requested_mtype_mode = requested_mtype_mode;
    info.pardiso_info.selected_mtype_mode.clear();
    info.pardiso_info.selected_mtype = 0;
    info.pardiso_info.matrix_part = "upper";
    info.pardiso_info.attempt_count = 0;
    info.pardiso_info.retried_from_spd_to_sym_indef = false;
    info.pardiso_info.phase11_error = 0;
    info.pardiso_info.phase22_error = 0;
    info.pardiso_info.phase33_error = 0;
    info.pardiso_info.factor_nnz = -1;
    info.pardiso_info.peak_memory_kb = -1;
    info.pardiso_info.export_upper_csr_time = 0.0;

#ifdef STAPPP_ENABLE_MKL_PARDISO
    struct AttemptConfig
    {
        const char* mode;
        int mtype;
    };

    vector<AttemptConfig> attempts;
    if (requested_mtype_mode == "spd")
        attempts.push_back(AttemptConfig{"spd", 2});
    else if (requested_mtype_mode == "sym-indef")
        attempts.push_back(AttemptConfig{"sym-indef", -2});
    else if (requested_mtype_mode == "unsym")
        throw runtime_error("*** Error *** Unsymmetric PARDISO mode is not supported by the symmetric half-storage mainline.");
    else
        attempts.push_back(AttemptConfig{"spd", 2});

    for (size_t attempt_index = 0; attempt_index < attempts.size(); ++attempt_index)
    {
        const AttemptConfig& config = attempts[attempt_index];
        ++info.pardiso_info.attempt_count;
        info.pardiso_info.selected_mtype_mode = config.mode;
        info.pardiso_info.selected_mtype = config.mtype;

        double export_upper_csr_time = 0.0;
        PardisoAttemptResult attempt =
            run_pardiso_attempt(matrix, rhs, config.mtype, &export_upper_csr_time, solution);
        info.pardiso_info.phase11_error = attempt.phase11_error;
        info.pardiso_info.phase22_error = attempt.phase22_error;
        info.pardiso_info.phase33_error = attempt.phase33_error;
        info.pardiso_info.factor_nnz = attempt.factor_nnz;
        info.pardiso_info.peak_memory_kb = attempt.peak_memory_kb;
        info.pardiso_info.export_upper_csr_time = export_upper_csr_time;

        if (attempt.phase11_error == 0 &&
            attempt.phase22_error == 0 &&
            attempt.phase33_error == 0)
        {
            CSparseSolveResult result;
            result.converged = true;
            result.iterations = 0;
            vector<double> zero_solution(matrix.dim(), 0.0);
            result.initial_relative_residual = matrix.RelativeResidual(rhs, zero_solution.data());
            result.final_relative_residual = matrix.RelativeResidual(rhs, solution);
            return result;
        }

        const bool should_retry_from_spd =
            requested_mtype_mode == "auto" &&
            string(config.mode) == "spd" &&
            attempt.phase11_error == 0 &&
            attempt.phase22_error == -4;
        if (should_retry_from_spd)
        {
            info.pardiso_info.retried_from_spd_to_sym_indef = true;
            attempts.push_back(AttemptConfig{"sym-indef", -2});
            continue;
        }

        ostringstream msg;
        if (attempt.phase11_error != 0)
            msg << "*** Error *** PARDISO symbolic factorization failed with error code " << attempt.phase11_error;
        else if (attempt.phase22_error != 0)
            msg << "*** Error *** PARDISO numeric factorization failed with error code " << attempt.phase22_error;
        else
            msg << "*** Error *** PARDISO solve phase failed with error code " << attempt.phase33_error;
        msg << " (mtype=" << config.mode << ", matrix_part=upper)";
        throw runtime_error(msg.str());
    }

    throw runtime_error("*** Error *** PARDISO failed without producing a valid attempt result.");
#else
    (void)matrix;
    (void)rhs;
    (void)solution;
    throw runtime_error("*** Error *** PARDISO backend requested, but STAPPP_ENABLE_MKL_PARDISO is OFF. Reconfigure with MKL support on a machine that has Intel oneMKL installed.");
#endif
}

// LDLT facterization
void CLDLTSolver::LDLT()
{
	unsigned int N = K.dim();
	unsigned int* ColumnHeights = K.GetColumnHeights();   // Column Hights
	unsigned int* DiagonalAddress = K.GetDiagonalAddress();
	double* data = K.GetData();

	// Compute max initial diagonal for relative pivot tolerance
	double max_diag = 0.0;
	for (unsigned int i = 1; i <= N; i++)
	{
		double d = fabs(K(i,i));
		if (d > max_diag) max_diag = d;
	}

	(void)max_diag;
	double pivot_tol = 1.0e-10;

	double min_pivot = max_diag;
	unsigned int min_pivot_eq = 0;

	for (unsigned int j = 2; j <= N; j++)
	{
		unsigned int mj = j - ColumnHeights[j-1];
		const unsigned int diag_j = DiagonalAddress[j - 1] - 1;

		for (unsigned int i = mj+1; i <= j-1; i++)
		{
			unsigned int mi = i - ColumnHeights[i-1];
			unsigned int first = max(mi, mj);
			if (first > i - 1)
				continue;

			double C = 0.0;
			unsigned int off_i = DiagonalAddress[i - 1] + (i - first) - 1;
			unsigned int off_j = DiagonalAddress[j - 1] + (j - first) - 1;
			for (unsigned int r = first; r <= i-1; r++, off_i--, off_j--)
				C += data[off_i] * data[off_j];

			data[DiagonalAddress[j - 1] + (j - i) - 1] -= C;
		}

		for (unsigned int r = mj; r <= j-1; r++)
		{
			unsigned int off_rj = DiagonalAddress[j - 1] + (j - r) - 1;
			double Lrj = data[off_rj] / data[DiagonalAddress[r - 1] - 1];
			data[diag_j] -= Lrj * data[off_rj];
			data[off_rj] = Lrj;
		}

		double abs_pivot = fabs(data[diag_j]);
		if (abs_pivot < min_pivot)
		{
			min_pivot = abs_pivot;
			min_pivot_eq = j;
		}

		if (abs_pivot <= pivot_tol)
		{
			cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
				 << "    Equation no = " << j << endl
				 << "    Pivot = " << data[diag_j] << endl
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
	unsigned int* ColumnHeights = K.GetColumnHeights();
	unsigned int* DiagonalAddress = K.GetDiagonalAddress();
	double* data = K.GetData();

	for (unsigned int i = 2; i <= N; i++)
	{
		unsigned int mi = i - ColumnHeights[i-1];

		unsigned int off = DiagonalAddress[i - 1] + (i - mi) - 1;
		for (unsigned int j = mi; j <= i-1; j++)
			Force[i-1] -= data[off--] * Force[j-1];
	}

	for (unsigned int i = 1; i <= N; i++)
		Force[i-1] /= data[DiagonalAddress[i - 1] - 1];

	for (unsigned int j = N; j >= 2; j--)
	{
		unsigned int mj = j - ColumnHeights[j-1];

		unsigned int off = DiagonalAddress[j - 1] + (j - mj) - 1;
		for (unsigned int i = mj; i <= j-1; i++)
			Force[i-1] -= data[off--] * Force[j-1];
	}
};
