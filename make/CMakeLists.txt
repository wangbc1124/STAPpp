cmake_minimum_required(VERSION 3.10)
project(stappp_bridge_fast CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(STAPPP_DEBUG "Output detailed results for debugging." OFF)
option(STAPPP_ENABLE_MKL_PARDISO "Enable Intel oneMKL PARDISO sparse backend." OFF)
set(STAPPP_MKL_ROOT "$ENV{STAPPP_MKL_ROOT}" CACHE PATH "Path to Intel oneMKL root directory")
if(NOT STAPPP_MKL_ROOT AND DEFINED ENV{MKLROOT})
    set(STAPPP_MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "Path to Intel oneMKL root directory" FORCE)
endif()
set(STAPPP_MKL_INCLUDE_DIR "$ENV{STAPPP_MKL_INCLUDE_DIR}" CACHE PATH "Path to Intel oneMKL include directory")
set(STAPPP_MKL_LIBRARY "$ENV{STAPPP_MKL_LIBRARY}" CACHE FILEPATH "Path to Intel oneMKL runtime import library (for example mkl_rt.lib)")
set(STAPPP_MKL_RUNTIME_DIR "$ENV{STAPPP_MKL_RUNTIME_DIR}" CACHE PATH "Path to Intel oneMKL runtime DLL directory")
set(STAPPP_COMPILER_RUNTIME_DIR "$ENV{STAPPP_COMPILER_RUNTIME_DIR}" CACHE PATH "Path to compiler runtime DLL directory")
set(STAPPP_EXTRA_RUNTIME_DIRS "$ENV{STAPPP_EXTRA_RUNTIME_DIRS}" CACHE STRING "Additional runtime directories to add before running Bridge targets")
if(STAPPP_DEBUG)
    add_compile_definitions(_DEBUG_)
endif()

file(GLOB STAPPP_SOURCES CONFIGURE_DEPENDS src/cpp/*.cpp)
file(GLOB STAPPP_HEADERS CONFIGURE_DEPENDS src/h/*.h)

add_executable(stap++ ${STAPPP_SOURCES} ${STAPPP_HEADERS})
target_include_directories(stap++ PRIVATE src/h)

if(MSVC)
    target_compile_options(stap++ PRIVATE /O2)
else()
    target_compile_options(stap++ PRIVATE -O2)
endif()

if(STAPPP_ENABLE_MKL_PARDISO)
    if(NOT STAPPP_MKL_ROOT AND STAPPP_MKL_INCLUDE_DIR)
        get_filename_component(_stappp_mkl_root "${STAPPP_MKL_INCLUDE_DIR}" DIRECTORY)
        set(STAPPP_MKL_ROOT "${_stappp_mkl_root}" CACHE PATH "Path to Intel oneMKL root directory" FORCE)
    endif()
    if(NOT STAPPP_MKL_INCLUDE_DIR AND STAPPP_MKL_ROOT)
        set(STAPPP_MKL_INCLUDE_DIR "${STAPPP_MKL_ROOT}/include" CACHE PATH "Path to Intel oneMKL include directory" FORCE)
    endif()
    if(NOT STAPPP_MKL_LIBRARY AND STAPPP_MKL_ROOT)
        set(STAPPP_MKL_LIBRARY "${STAPPP_MKL_ROOT}/lib/mkl_rt.lib" CACHE FILEPATH "Path to Intel oneMKL runtime import library" FORCE)
    endif()
    if(NOT STAPPP_MKL_INCLUDE_DIR OR NOT EXISTS "${STAPPP_MKL_INCLUDE_DIR}")
        message(FATAL_ERROR "STAPPP_ENABLE_MKL_PARDISO=ON requires STAPPP_MKL_INCLUDE_DIR to point to a valid oneMKL include directory.")
    endif()
    if(NOT STAPPP_MKL_LIBRARY OR NOT EXISTS "${STAPPP_MKL_LIBRARY}")
        message(FATAL_ERROR "STAPPP_ENABLE_MKL_PARDISO=ON requires STAPPP_MKL_LIBRARY to point to a valid oneMKL runtime import library.")
    endif()
    if(NOT STAPPP_MKL_RUNTIME_DIR)
        if(WIN32 AND STAPPP_MKL_ROOT AND EXISTS "${STAPPP_MKL_ROOT}/redist/intel64")
            set(STAPPP_MKL_RUNTIME_DIR "${STAPPP_MKL_ROOT}/redist/intel64" CACHE PATH "Path to Intel oneMKL runtime DLL directory" FORCE)
        elseif(STAPPP_MKL_ROOT AND EXISTS "${STAPPP_MKL_ROOT}/bin")
            set(STAPPP_MKL_RUNTIME_DIR "${STAPPP_MKL_ROOT}/bin" CACHE PATH "Path to Intel oneMKL runtime DLL directory" FORCE)
        else()
            get_filename_component(_stappp_mkl_library_dir "${STAPPP_MKL_LIBRARY}" DIRECTORY)
            set(STAPPP_MKL_RUNTIME_DIR "${_stappp_mkl_library_dir}" CACHE PATH "Path to Intel oneMKL runtime DLL directory" FORCE)
        endif()
    endif()
    if(NOT EXISTS "${STAPPP_MKL_RUNTIME_DIR}")
        message(FATAL_ERROR "STAPPP_ENABLE_MKL_PARDISO=ON requires STAPPP_MKL_RUNTIME_DIR to point to a valid oneMKL runtime DLL directory.")
    endif()
    target_compile_definitions(stap++ PRIVATE STAPPP_ENABLE_MKL_PARDISO=1)
    target_include_directories(stap++ PRIVATE "${STAPPP_MKL_INCLUDE_DIR}")
    target_link_libraries(stap++ PRIVATE "${STAPPP_MKL_LIBRARY}")
endif()

set(STAPPP_RUN_ENV)
if(NOT STAPPP_COMPILER_RUNTIME_DIR)
    get_filename_component(_stappp_compiler_runtime_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
    set(STAPPP_COMPILER_RUNTIME_DIR "${_stappp_compiler_runtime_dir}" CACHE PATH "Path to compiler runtime DLL directory" FORCE)
endif()

if(STAPPP_ENABLE_MKL_PARDISO)
    set(_stappp_runtime_dirs "${STAPPP_MKL_RUNTIME_DIR}" "${STAPPP_COMPILER_RUNTIME_DIR}")
    if(STAPPP_EXTRA_RUNTIME_DIRS)
        list(APPEND _stappp_runtime_dirs ${STAPPP_EXTRA_RUNTIME_DIRS})
    endif()
    if(WIN32)
        string(REPLACE ";" "\\;" STAPPP_RUN_PATH "${_stappp_runtime_dirs}")
        set(STAPPP_RUN_ENV "PATH=${STAPPP_RUN_PATH}")
    else()
        string(REPLACE ";" ":" STAPPP_RUN_LIBRARY_PATH "${_stappp_runtime_dirs}")
        set(STAPPP_RUN_ENV "LD_LIBRARY_PATH=${STAPPP_RUN_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
    endif()
endif()

file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/results")

function(stappp_add_bridge_run_target target_name input_file csv_file pardiso_mtype)
    add_custom_target(${target_name}
        COMMAND ${CMAKE_COMMAND} -E env ${STAPPP_RUN_ENV}
            $<TARGET_FILE:stap++>
            ${input_file}
            --solver sparse-auto
            --backend pardiso
            --pardiso-mtype ${pardiso_mtype}
            --output summary
            --tol 1e-6
            --max-iter 5000
            --csv ${csv_file}
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        DEPENDS stap++
        USES_TERMINAL
        VERBATIM
    )
endfunction()

stappp_add_bridge_run_target(run-bridge-1 "Bridge-1/Bridge-1.dat" "results/Bridge-1.displacements.csv" "auto")
stappp_add_bridge_run_target(run-bridge-2 "Bridge-2/Bridge-2.dat" "results/Bridge-2.displacements.csv" "sym-indef")
stappp_add_bridge_run_target(run-bridge-3 "Bridge-3/Bridge-3.dat" "results/Bridge-3.displacements.csv" "auto")
stappp_add_bridge_run_target(run-bridge-4 "Bridge-4.dat" "results/Bridge-4.displacements.csv" "auto")
