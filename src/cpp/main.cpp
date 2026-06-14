/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"
#include <memory>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <exception>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName [--solver skyline|sparse-cg|sparse-bicgstab|sparse-gmres|sparse-auto] [--backend standard|pardiso] [--pardiso-mtype auto|spd|sym-indef|unsym] [--output full|summary] [--precond none|jacobi|ssor|block-jacobi|ilu0] [--scale none|diag] [--tol value] [--max-iter value] [--residual-print-interval value] [--csv file]\n";
		exit(1);
	}

    string solver = "sparse-auto";
    string backend = IsSparseBackendAvailable("pardiso") ? "pardiso" : "standard";
    string pardiso_mtype = "auto";
    string output_mode = "full";
    double pcg_tolerance = 1.0e-6;
    unsigned int pcg_max_iterations = 5000;
    unsigned int residual_print_interval = 50;
    string preconditioner = "ssor";
    string scale_mode = "diag";
    string csv_file;

    const char* env_solver = getenv("STAP_SOLVER");
    if (env_solver)
        solver = env_solver;
    const char* env_backend = getenv("STAP_BACKEND");
    if (env_backend)
        backend = env_backend;
    const char* env_output = getenv("STAP_OUTPUT_MODE");
    if (env_output)
        output_mode = env_output;
    const char* env_tol = getenv("STAP_PCG_TOL");
    if (env_tol)
        pcg_tolerance = atof(env_tol);
    const char* env_max_iter = getenv("STAP_PCG_MAX_ITER");
    if (env_max_iter)
        pcg_max_iterations = static_cast<unsigned int>(atoi(env_max_iter));
    const char* env_precond = getenv("STAP_PCG_PRECOND");
    if (env_precond)
        preconditioner = env_precond;
    const char* env_scale = getenv("STAP_PCG_SCALE");
    if (env_scale)
        scale_mode = env_scale;

    for (int arg = 2; arg < argc; ++arg)
    {
        string opt(argv[arg]);
        if (opt == "--solver" && arg + 1 < argc)
            solver = argv[++arg];
        else if (opt == "--backend" && arg + 1 < argc)
            backend = argv[++arg];
        else if (opt == "--pardiso-mtype" && arg + 1 < argc)
            pardiso_mtype = argv[++arg];
        else if (opt == "--output" && arg + 1 < argc)
            output_mode = argv[++arg];
        else if (opt == "--precond" && arg + 1 < argc)
            preconditioner = argv[++arg];
        else if (opt == "--scale" && arg + 1 < argc)
            scale_mode = argv[++arg];
        else if (opt == "--tol" && arg + 1 < argc)
            pcg_tolerance = atof(argv[++arg]);
        else if (opt == "--max-iter" && arg + 1 < argc)
            pcg_max_iterations = static_cast<unsigned int>(atoi(argv[++arg]));
        else if (opt == "--residual-print-interval" && arg + 1 < argc)
            residual_print_interval = static_cast<unsigned int>(atoi(argv[++arg]));
        else if (opt == "--csv" && arg + 1 < argc)
            csv_file = argv[++arg];
        else
        {
            cout << "*** Error *** Unknown or incomplete option: " << opt << endl;
            exit(1);
        }
    }

    if (solver != "skyline" && solver != "sparse-cg" && solver != "sparse-bicgstab" &&
        solver != "sparse-gmres" && solver != "sparse-auto")
    {
        cout << "*** Error *** Unsupported solver: " << solver << endl;
        exit(1);
    }
    if (backend != "standard" && backend != "pardiso")
    {
        cout << "*** Error *** Unsupported sparse backend: " << backend << endl;
        exit(1);
    }
    if (pardiso_mtype != "auto" && pardiso_mtype != "spd" &&
        pardiso_mtype != "sym-indef" && pardiso_mtype != "unsym")
    {
        cout << "*** Error *** Unsupported PARDISO mtype mode: " << pardiso_mtype << endl;
        exit(1);
    }
    if (output_mode != "full" && output_mode != "summary")
    {
        cout << "*** Error *** Unsupported output mode: " << output_mode << endl;
        exit(1);
    }
    if (preconditioner != "none" && preconditioner != "jacobi" &&
        preconditioner != "ssor" && preconditioner != "block-jacobi" &&
        preconditioner != "ilu0")
    {
        cout << "*** Error *** Unsupported preconditioner: " << preconditioner << endl;
        exit(1);
    }
    if (scale_mode != "none" && scale_mode != "diag")
    {
        cout << "*** Error *** Unsupported scale mode: " << scale_mode << endl;
        exit(1);
    }
    if (residual_print_interval == 0)
    {
        cout << "*** Error *** --residual-print-interval must be positive." << endl;
        exit(1);
    }

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";
    if (csv_file.empty())
        csv_file = filename + ".displacements.csv";

	CDomain* FEMData = CDomain::GetInstance();

    Clock timer;
    timer.Start();

    COutputter* Output = COutputter::GetInstance(OutFile);
    Output->SetSummaryMode(output_mode == "summary");

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();
    double time_mpc_alias = FEMData->GetLastMpcAliasTime();

    if (!FEMData->GetMODEX())
    {
        *Output << "Data check completed !" << endl << endl;
        return 0;
    }

    double time_assemble = 0.0;
    double time_csr_pattern = 0.0;
    double time_csr_assembly = 0.0;
    double time_csr_element_assembly = 0.0;
    double time_csr_mpc_assembly = 0.0;
    double time_element_stiffness = 0.0;
    double time_pattern_insert = 0.0;
    double time_value_insert = 0.0;
    double time_active_dof_pack = 0.0;
    double time_mpc_pattern = 0.0;
    double time_export_upper_csr = 0.0;
    double time_symmetry_check = 0.0;
    double time_scaling = 0.0;
    double time_precond_setup = 0.0;
    double time_iter_solve = 0.0;
    double time_residual_check = 0.0;
    double time_csv_write = 0.0;

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif

    if (solver == "sparse-cg" || solver == "sparse-bicgstab" || solver == "sparse-gmres" || solver == "sparse-auto")
    {
        if (backend == "pardiso")
        {
            try
            {
                CSparseSymmetricMatrix* SparseMatrix = FEMData->AssemblePardisoStiffnessMatrix();
                time_csr_pattern = FEMData->GetLastCsrPatternTime();
                time_csr_assembly = FEMData->GetLastCsrAssemblyTime();
                time_csr_element_assembly = FEMData->GetLastElementCsrAssemblyTime();
                time_csr_mpc_assembly = FEMData->GetLastMpcCsrAssemblyTime();
                time_element_stiffness = FEMData->GetLastElementStiffnessTime();
                time_pattern_insert = FEMData->GetLastPatternInsertTime();
                time_value_insert = FEMData->GetLastValueInsertTime();
                time_active_dof_pack = FEMData->GetLastActiveDofPackTime();
                time_mpc_pattern = FEMData->GetLastMpcPatternTime();
                time_assemble = timer.ElapsedTime();

                CSparseBackendInfo backend_info;
                *Output << " S P A R S E   I T E R A T I V E   S O L V E R" << endl << endl
                        << "     BACKEND = pardiso" << endl
                        << "     BACKEND_AVAILABLE = " << (IsSparseBackendAvailable("pardiso") ? "YES" : "NO") << endl
                        << "     SOLVER = " << solver << endl
                        << "     ACTUAL SOLVER = pardiso" << endl
                        << "     TOLERANCE = " << pcg_tolerance << endl
                        << "     MAX ITERATIONS = " << pcg_max_iterations << endl
                        << "     RESIDUAL PRINT INTERVAL = " << residual_print_interval << endl
                        << "     PRECONDITIONER = direct" << endl
                        << "     SCALE = none" << endl
                        << "     SYMMETRY_ERROR = 0" << endl
                        << "     NUMERICALLY SYMMETRIC = YES" << endl
                        << "     INITIAL GUESS = zero vector" << endl << endl;

                for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); ++lcase)
                {
                    FEMData->AssembleForce(lcase + 1);
                    vector<double> rhs(FEMData->GetNEQ(), 0.0);
                    for (unsigned int i = 0; i < FEMData->GetNEQ(); ++i)
                        rhs[i] = FEMData->GetForce()[i];

                    Clock solveTimer;
                    solveTimer.Start();
                    CSparseSolveResult result =
                        SolvePardisoSystem(*SparseMatrix, rhs.data(), FEMData->GetForce(), pardiso_mtype, backend_info);
                    time_iter_solve += solveTimer.ElapsedTime();

                    Clock residualTimer;
                    residualTimer.Start();
                    const double checked_residual = SparseMatrix->RelativeResidual(rhs.data(), FEMData->GetForce());
                    time_residual_check += residualTimer.ElapsedTime();
                    time_export_upper_csr = backend_info.pardiso_info.export_upper_csr_time;

                    *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl
                            << "     ITERATIVE SOLVER CONVERGED = " << (result.converged ? "YES" : "NO") << endl
                            << "     ITERATIONS = " << result.iterations << endl
                            << "     INITIAL RELATIVE RESIDUAL = " << result.initial_relative_residual << endl
                            << "     FINAL RELATIVE RESIDUAL = " << result.final_relative_residual << endl
                            << "     CHECKED RELATIVE RESIDUAL = " << checked_residual << endl << endl
                            << "     PARDISO_REQUESTED_MTYPE = " << backend_info.pardiso_info.requested_mtype_mode << endl
                            << "     PARDISO_SELECTED_MTYPE = " << backend_info.pardiso_info.selected_mtype_mode << endl
                            << "     PARDISO_SELECTED_MTYPE_CODE = " << backend_info.pardiso_info.selected_mtype << endl
                            << "     PARDISO_MATRIX_PART = " << backend_info.pardiso_info.matrix_part << endl
                            << "     PARDISO_ATTEMPT_COUNT = " << backend_info.pardiso_info.attempt_count << endl
                            << "     PARDISO_RETRY_FROM_SPD_TO_SYM_INDEF = "
                            << (backend_info.pardiso_info.retried_from_spd_to_sym_indef ? "YES" : "NO") << endl
                            << "     PARDISO_PHASE11_ERROR = " << backend_info.pardiso_info.phase11_error << endl
                            << "     PARDISO_PHASE22_ERROR = " << backend_info.pardiso_info.phase22_error << endl
                            << "     PARDISO_PHASE33_ERROR = " << backend_info.pardiso_info.phase33_error << endl
                            << "     PARDISO_FACT_NNZ = " << backend_info.pardiso_info.factor_nnz << endl
                            << "     PARDISO_PEAK_MEMORY_KB = " << backend_info.pardiso_info.peak_memory_kb << endl
                            << endl;

                    if (!result.converged)
                    {
                        cerr << "*** Error *** PARDISO did not produce a converged result." << endl;
                        delete SparseMatrix;
                        exit(5);
                    }
                    if (!isfinite(checked_residual) || checked_residual > pcg_tolerance)
                    {
                        cerr << "*** Error *** Checked relative residual exceeds tolerance: "
                             << checked_residual << " > " << pcg_tolerance << endl;
                        delete SparseMatrix;
                        exit(5);
                    }

                    Clock csvTimer;
                    csvTimer.Start();
                    if (!FEMData->WriteDisplacementCSV(csv_file))
                    {
                        cerr << "*** Error *** Failed to write displacement CSV: " << csv_file << endl;
                        delete SparseMatrix;
                        exit(6);
                    }
                    time_csv_write += csvTimer.ElapsedTime();
                    *Output << "     DISPLACEMENT CSV = " << csv_file << endl << endl;
                }

                delete SparseMatrix;
            }
            catch (const exception& error)
            {
                cerr << error.what() << endl;
                exit(5);
            }
        }
        else
        {
            CCSRMatrix* SparseMatrix = FEMData->AssembleSparseStiffnessMatrix(backend);
            time_csr_pattern = FEMData->GetLastCsrPatternTime();
            time_csr_assembly = FEMData->GetLastCsrAssemblyTime();
            time_csr_element_assembly = FEMData->GetLastElementCsrAssemblyTime();
            time_csr_mpc_assembly = FEMData->GetLastMpcCsrAssemblyTime();
            time_assemble = timer.ElapsedTime();
            vector<vector<unsigned int> > blocks;
            const vector<vector<unsigned int> >* equation_blocks = 0;
            if (preconditioner == "block-jacobi" || (solver == "sparse-auto" && preconditioner == "ssor"))
            {
                blocks = FEMData->BuildNodeEquationBlocks();
                equation_blocks = &blocks;
            }

            CSparseSolverOptions sparse_options;
            sparse_options.requested_solver = solver;
            sparse_options.backend_name = backend;
            sparse_options.requested_preconditioner = preconditioner;
            sparse_options.pardiso_mtype_mode = pardiso_mtype;
            sparse_options.scale_mode = scale_mode;
            sparse_options.tolerance = pcg_tolerance;
            sparse_options.max_iterations = pcg_max_iterations;
            sparse_options.residual_print_interval = residual_print_interval;

            unique_ptr<CSparseSolverBackend> sparse_backend;
            CSparseBackendInfo backend_info;
            try
            {
                sparse_backend = CreateSparseSolverBackend(backend);
                sparse_backend->Setup(*SparseMatrix, equation_blocks, sparse_options);
                backend_info = sparse_backend->GetInfo();
            }
            catch (const exception& error)
            {
                cerr << error.what() << endl;
                delete SparseMatrix;
                exit(5);
            }
            time_symmetry_check = backend_info.timings.symmetry_check_time;
            time_scaling = backend_info.timings.scaling_time;
            time_precond_setup = backend_info.timings.preconditioner_setup_time;

            *Output << " S P A R S E   I T E R A T I V E   S O L V E R" << endl << endl
                    << "     BACKEND = " << backend_info.backend_name << endl
                    << "     BACKEND_AVAILABLE = " << (backend_info.backend_available ? "YES" : "NO") << endl
                    << "     SOLVER = " << solver << endl
                    << "     ACTUAL SOLVER = " << backend_info.actual_solver << endl
                    << "     TOLERANCE = " << pcg_tolerance << endl
                    << "     MAX ITERATIONS = " << pcg_max_iterations << endl
                    << "     RESIDUAL PRINT INTERVAL = " << residual_print_interval << endl
                    << "     PRECONDITIONER = " << backend_info.actual_preconditioner << endl
                    << "     SCALE = " << scale_mode << endl
                    << "     SYMMETRY_ERROR = " << backend_info.symmetry.symmetry_error << endl
                    << "     NUMERICALLY SYMMETRIC = "
                    << (backend_info.symmetry.numerically_symmetric ? "YES" : "NO") << endl
                    << "     INITIAL GUESS = zero vector" << endl << endl;

            if (backend_info.scaling_enabled)
            {
                const CDiagonalScaling::Stats& stats = backend_info.scaling_stats;
                *Output << "     MIN_ABS_DIAG = " << stats.min_abs_diag << endl
                        << "     MAX_ABS_DIAG = " << stats.max_abs_diag << endl
                        << "     MAX_SCALE = " << stats.max_scale << endl
                        << "     MIN_SCALE = " << stats.min_scale << endl
                        << "     SMALL_OR_ZERO_DIAGONAL_COUNT = " << stats.number_of_small_or_zero_diagonal_entries << endl
                        << endl;
            }

            if (backend_info.has_block_jacobi_stats)
            {
                const CBlockJacobiPreconditioner::Stats& stats = backend_info.block_jacobi_stats;
                *Output << "     BLOCK_JACOBI_NUM_BLOCKS = " << stats.num_blocks << endl
                        << "     BLOCK_SIZE_HISTOGRAM = ["
                        << stats.block_size_histogram[0] << ", "
                        << stats.block_size_histogram[1] << ", "
                        << stats.block_size_histogram[2] << ", "
                        << stats.block_size_histogram[3] << ", "
                        << stats.block_size_histogram[4] << ", "
                        << stats.block_size_histogram[5] << "]" << endl
                        << "     BLOCK_JACOBI_CHOLESKY_SUCCESS = " << stats.num_cholesky_success << endl
                        << "     BLOCK_JACOBI_SHIFTED_BLOCKS = " << stats.num_shifted_blocks << endl
                        << "     BLOCK_JACOBI_DIAG_FALLBACK_BLOCKS = " << stats.num_diag_fallback_blocks << endl
                        << endl;
            }
            else if (backend_info.has_ilu0_stats)
            {
                const CILU0Preconditioner::Stats& stats = backend_info.ilu0_stats;
                *Output << "     ILU0_SHIFTED_PIVOTS = " << stats.shifted_pivots << endl
                        << "     ILU0_ZERO_PATTERN_FALLBACKS = " << stats.zero_pattern_fallbacks << endl
                        << endl;
            }

            for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
            {
                FEMData->AssembleForce(lcase + 1);
                vector<double> rhs(FEMData->GetNEQ(), 0.0);
                for (unsigned int i = 0; i < FEMData->GetNEQ(); ++i)
                    rhs[i] = FEMData->GetForce()[i];

                CSparseSolveResult result;
                double checked_residual = 0.0;
                try
                {
                    Clock solveTimer;
                    solveTimer.Start();
                    result = sparse_backend->Solve(rhs.data(), FEMData->GetForce());
                    time_iter_solve += solveTimer.ElapsedTime();
                    backend_info = sparse_backend->GetInfo();
                }
                catch (const exception& error)
                {
                    cerr << error.what() << endl;
                    delete SparseMatrix;
                    exit(5);
                }

                Clock residualTimer;
                residualTimer.Start();
                checked_residual = SparseMatrix->RelativeResidual(rhs.data(), FEMData->GetForce());
                time_residual_check += residualTimer.ElapsedTime();

                *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl
                        << "     ITERATIVE SOLVER CONVERGED = " << (result.converged ? "YES" : "NO") << endl
                        << "     ITERATIONS = " << result.iterations << endl
                        << "     INITIAL RELATIVE RESIDUAL = " << result.initial_relative_residual << endl
                        << "     FINAL RELATIVE RESIDUAL = " << result.final_relative_residual << endl
                        << "     CHECKED RELATIVE RESIDUAL = " << checked_residual << endl << endl;

                if (!result.converged)
                {
                    cerr << "*** Error *** Iterative solver did not converge within max iterations." << endl;
                    delete SparseMatrix;
                    exit(5);
                }
                if (!isfinite(checked_residual) || checked_residual > pcg_tolerance)
                {
                    cerr << "*** Error *** Checked relative residual exceeds tolerance: "
                         << checked_residual << " > " << pcg_tolerance << endl;
                    delete SparseMatrix;
                    exit(5);
                }

                Clock csvTimer;
                csvTimer.Start();
                if (!FEMData->WriteDisplacementCSV(csv_file))
                {
                    cerr << "*** Error *** Failed to write displacement CSV: " << csv_file << endl;
                    delete SparseMatrix;
                    exit(6);
                }
                time_csv_write += csvTimer.ElapsedTime();
                *Output << "     DISPLACEMENT CSV = " << csv_file << endl << endl;
            }

            delete SparseMatrix;
        }
    }
    else
    {
//      Allocate global vectors and matrices, such as the Force, ColumnHeights,
//      DiagonalAddress and StiffnessMatrix, and calculate the column heights
//      and address of diagonal elements
        FEMData->AllocateMatrices();

//      Assemble the banded gloabl stiffness matrix
        FEMData->AssembleStiffnessMatrix();

        time_assemble = timer.ElapsedTime();

//      Solve the linear equilibrium equations for displacements
        CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());

//      Perform L*D*L(T) factorization of stiffness matrix
        Solver->LDLT();

//      Loop over for all load cases
        for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
        {
//          Assemble righ-hand-side vector (force vector)
            FEMData->AssembleForce(lcase + 1);

//          Reduce right-hand-side force vector and back substitute
            Solver->BackSubstitution(FEMData->GetForce());

            *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;

#ifdef _DEBUG_
            Output->PrintDisplacement();
#endif

            Output->OutputNodalDisplacement();

//          Calculate and output stresses of all elements
            Output->OutputElementStress();

            Clock csvTimer;
            csvTimer.Start();
            if (!FEMData->WriteDisplacementCSV(csv_file))
            {
                cerr << "*** Error *** Failed to write displacement CSV: " << csv_file << endl;
                delete Solver;
                exit(6);
            }
            time_csv_write += csvTimer.ElapsedTime();
            *Output << "     DISPLACEMENT CSV = " << csv_file << endl << endl;
        }

        delete Solver;
    }

    double time_solution = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     READ_TIME = " << time_input - time_mpc_alias << endl
            << "     MPC_ALIAS_TIME = " << time_mpc_alias << endl
            << "     CSR_PATTERN_TIME = " << time_csr_pattern << endl
            << "     CSR_ASSEMBLY_TIME = " << time_csr_assembly << endl
            << "     ELEMENT_VALUE_ASSEMBLY_TIME = " << time_csr_element_assembly << endl
            << "     MPC_VALUE_ASSEMBLY_TIME = " << time_csr_mpc_assembly << endl
            << "     ELEMENT_STIFFNESS_TIME = " << time_element_stiffness << endl
            << "     PATTERN_INSERT_TIME = " << time_pattern_insert << endl
            << "     VALUE_INSERT_TIME = " << time_value_insert << endl
            << "     ACTIVE_DOF_PACK_TIME = " << time_active_dof_pack << endl
            << "     MPC_PATTERN_TIME = " << time_mpc_pattern << endl
            << "     EXPORT_UPPER_CSR_TIME = " << time_export_upper_csr << endl
            << "     SYMMETRY_CHECK_TIME = " << time_symmetry_check << endl
            << "     SCALING_TIME = " << time_scaling << endl
            << "     PRECOND_SETUP_TIME = " << time_precond_setup << endl
            << "     ITER_SOLVE_TIME = " << time_iter_solve << endl
            << "     RESIDUAL_CHECK_TIME = " << time_residual_check << endl
            << "     CSV_WRITE_TIME = " << time_csv_write << endl
            << "     TOTAL_TIME = " << time_solution << endl << endl;

	return 0;
}
