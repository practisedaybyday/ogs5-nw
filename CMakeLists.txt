#### OGS-5 Project ####

# Specify minimum CMake version
cmake_minimum_required(VERSION 2.6)	# or set to 2.4 for older systems (like JUGENE...)

# Project name is set later
PROJECT(tmp)

### CMake setup ###
INCLUDE(CMakeConfiguration/CMakeSetup.cmake)

### Find directories and libraries ###
INCLUDE (CMakeConfiguration/Find.cmake)

###############
### Options ###
###############

# Configurations
OPTION(OGS_FEM "Default configuration" OFF)
OPTION(OGS_FEM_SP "Sparse matrix configuration" OFF)
OPTION(OGS_FEM_JFNK "Jacobi free Newton-Krylov method for H2M problem" OFF)
OPTION(OGS_FEM_GEMS "GEMS configuration" OFF)
OPTION(OGS_FEM_PQC "Phreeqc configuration" OFF)
OPTION(OGS_FEM_BRNS "BRNS configuration" OFF)
OPTION(OGS_FEM_CHEMAPP "CHEMAPP configuration" OFF)
OPTION(OGS_USE_QT "User interface configuration" OFF)
OPTION(OGS_FEM_LIS "Library of Iterative Solvers for Linear Systems configuration" OFF)
OPTION(OGS_FEM_MKL "Math kernel library configuration" OFF)

IF (UNIX) # Only supported on Linux
	OPTION(OGS_FEM_MPI "Message passing interface configuration" OFF)
ENDIF (UNIX)

# Set default configuration when no other config is given
IF (NOT OGS_FEM AND NOT OGS_FEM_SP AND NOT OGS_FEM_GEMS AND NOT OGS_FEM_PQC AND NOT OGS_FEM_BRNS AND NOT OGS_FEM_CHEMAPP AND NOT OGS_FEM_LIS AND NOT OGS_FEM_MKL AND NOT OGS_FEM_MPI AND NOT OGS_USE_QT)
	MESSAGE (STATUS "No configuration specified. Default confuguration is used.")
	SET (OGS_FEM ON)
ENDIF (NOT OGS_FEM AND NOT OGS_FEM_SP AND NOT OGS_FEM_GEMS AND NOT OGS_FEM_PQC AND NOT OGS_FEM_BRNS AND NOT OGS_FEM_CHEMAPP AND NOT OGS_FEM_LIS AND NOT OGS_FEM_MKL AND NOT OGS_FEM_MPI AND NOT OGS_USE_QT)

# Check if a valid OGS configuration is given (see Macros.cmake)
CHECK_CONFIG()

# Additional options
OPTION(OGS_PACKAGING "Creating installers / packages" OFF)
IF (OGS_PACKAGING)
	OPTION (OGS_PACKAGING_ZIP "Do you want to package as zip?" OFF)
ENDIF (OGS_PACKAGING)
OPTION(OGS_PYTHON "Enable python scripting interface" OFF)
OPTION(OGS_BUILD_UTILITIES "Build additional utility programs?" OFF)
OPTION(OGS_CMAKE_DEBUG "Show additional cmake debug information" OFF)
OPTION(OGS_ONE_BIG_GTEST "Should all gtests be collected to one ctest" ON)
OPTION(OGS_USE_NETCDF "NetCDF interface" OFF)
OPTION(OGS_USE_OPENSG "VTK to OpenSG converter" OFF)
OPTION(OGS_USE_VRPN "Use VRPN for head tracking (experimental)" OFF)
OPTION(OGS_VRED_PLUGIN "Enable the Vred plugin (far from experimental)" OFF)
OPTION(OGS_BUILD_INFO "Should build information be generated" ON)
OPTION (CMAKE_CMD_ARGS "Build ogs5 storing the CMake command line args (hint: use cmake.ogs.sh" OFF )

MARK_AS_ADVANCED(FORCE OGS_PYTHON OGS_CMAKE_DEBUG OGS_ONE_BIG_GTEST OGS_USE_OPENSG OGS_USE_VRPN OGS_VRED_PLUGIN OGS_BUILD_INFO CMAKE_CMD_ARGS)

###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
SET (OGS_VERSION_MAJOR 5)
SET (OGS_VERSION_MINOR 1)
SET (OGS_VERSION_PATCH 07)
SET (OGS_RELEASE_PERSONS "WX")
SET (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}(${OGS_RELEASE_PERSONS})")
# should we use ISO dates: 2010-12-31 ?
SET (OGS_DATE "02.07.2011")

### Compiler setup ###
INCLUDE(CMakeConfiguration/CompilerSetup.cmake)

# Check for number of processors
CHECK_PROCESSOR_COUNT()

### General project setup ###
INCLUDE(CMakeConfiguration/GeneralProjectSetup.cmake)

IF(OGS_FEM)
   MESSAGE (STATUS  "Configuring for FEM command line" )	
ENDIF(OGS_FEM)

IF(OGS_FEM_SP)
	MESSAGE (STATUS  "Configuring for FEM command line with new sparse matrix" )	
	SET( PRJ_EXT "-SP" )
	ADD_DEFINITIONS(-DNEW_EQS)
ENDIF(OGS_FEM_SP)

IF(OGS_FEM_JFNK)
	MESSAGE (STATUS  "Configuring for Jacobi free Newton-Krylov method for H2M problem" )	
	SET( PRJ_EXT "-SP_JFNK" )
	ADD_DEFINITIONS(-DNEW_EQS)
	ADD_DEFINITIONS(-DJFNK_H2M)
ENDIF(OGS_FEM_JFNK)

IF(OGS_FEM_MPI)
	MESSAGE (STATUS  "Configuring for FEM command line with MPI" )
	STRING (REGEX MATCH mpi MPI_COMPILER_SET ${CMAKE_CXX_COMPILER})
	IF (NOT MPI_COMPILER_SET)
		MESSAGE (WARNING "Make sure that you specified the mpi c and c++ compiler on the command line. E.g. cmake -DOGS_FEM_MPI=ON -DCMAKE_C_COMPILER=/opt/openmpi-1.4.1/bin/mpicc -DCMAKE_CXX_COMPILER=/opt/openmpi-1.4.1/bin/mpic++ ..")
	ENDIF ()
	SET( PRJ_EXT "-MPI" )
	SET(PARALLEL_USE_MPI ON CACHE BOOL "Use Message Passing Interface")
	MARK_AS_ADVANCED(PARALLEL_USE_MPI)
	ADD_DEFINITIONS(-DUSE_MPI)
	ADD_DEFINITIONS(-DNEW_EQS)
ENDIF(OGS_FEM_MPI)

IF(OGS_FEM_GEMS)
	MESSAGE (STATUS  "Configuring for FEM command line with GEMS" )
	SET( PRJ_EXT "-GEMS" )
	ADD_DEFINITIONS(-DNO_ERROR_CONTROL)
	ADD_DEFINITIONS(-DGEM_REACT)
	ADD_DEFINITIONS(-DGEMlib)
	ADD_DEFINITIONS(-DIPMGEMPLUGIN)
	ADD_DEFINITIONS(-DNEW_EQS)
	IF (MSVC)
		SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
		SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	ENDIF (MSVC)

	IF(PARALLEL_USE_MPI)
		ADD_DEFINITIONS(-DUSE_MPI_GEMS)
	ENDIF(PARALLEL_USE_MPI)
ENDIF(OGS_FEM_GEMS)

IF(OGS_FEM_PQC)
	MESSAGE (STATUS  "Configuring for FEM command line with PQC" )	
	SET( PRJ_EXT "-PQC" )
	MESSAGE (STATUS "Configuration: g++ and LIBPHREEQC")
	ADD_DEFINITIONS(-D_POSIX_SOURCE)
	ADD_DEFINITIONS(-DNO_ERROR_CONTROL)
	ADD_DEFINITIONS(-DNEW_EQS)

	IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)  #WW
		IF (LIBPQC_DEBUG)
			MESSAGE (STATUS  "Compiling ogs and LIBPHREEQC with debug flags" )
			SET(CMAKE_CXX_FLAGS "-O0 -ggdb -Wall")
		ELSE (LIBPQC_DEBUG)
			ADD_DEFINITIONS ( "-O3" )
		ENDIF (LIBPQC_DEBUG)
	ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	SET (LIBPQC ON)
    ADD_SUBDIRECTORY ( PQC )
ENDIF(OGS_FEM_PQC)

IF(OGS_FEM_BRNS)
	MESSAGE (STATUS  "Configuring for FEM command line with BRNS" )
	SET( PRJ_EXT "-BRNS" )
	ADD_DEFINITIONS(-DNEW_EQS)
	ADD_DEFINITIONS(-DBRNS)
	ADD_DEFINITIONS(-DNO_ERROR_CONTROL)
	ADD_DEFINITIONS(-D_POSIX_SOURCE)
	IF(PARALLEL_USE_MPI)
			ADD_DEFINITIONS(-DUSE_MPI_BRNS)
	ENDIF(PARALLEL_USE_MPI)
ENDIF(OGS_FEM_BRNS)

IF(OGS_FEM_LIS)
	MESSAGE (STATUS  "Configuring for FEM command line with LIS" )
	SET( PRJ_EXT "-LIS" )
	ADD_DEFINITIONS(-DNEW_EQS)
	SET(PARALLEL_USE_OPENMP ON CACHE BOOL "Use Library of Iterative Solvers")
	MARK_AS_ADVANCED(PARALLEL_USE_OPENMP)
	SET(LIS ON CACHE BOOL "Use Library of Iterative Solvers")
	MARK_AS_ADVANCED(LIS)
	ADD_DEFINITIONS( -DIPMGEMPLUGIN )
ENDIF(OGS_FEM_LIS)

IF(OGS_FEM_MKL)
	MESSAGE (STATUS  "Configuring for FEM command line with MKL" )
	SET( PRJ_EXT "-MKL" )
	ADD_DEFINITIONS(-DNO_ERROR_CONTROL)
	ADD_DEFINITIONS(-DNEW_EQS)
	SET(PARALLEL_USE_OPENMP ON CACHE BOOL "Use Library of Iterative Solvers")
	MARK_AS_ADVANCED(PARALLEL_USE_OPENMP)
	SET(LIS ON CACHE BOOL "Use Library of Iterative Solvers")
	MARK_AS_ADVANCED(LIS)
	SET(MKL ON CACHE BOOL "Use Math Kernel Library")
	MARK_AS_ADVANCED(MKL)
ENDIF(OGS_FEM_MKL)

IF(OGS_FEM_CHEMAPP)
	MESSAGE (STATUS  "Configuring for FEM command line with CHEMAPP" )
	SET( PRJ_EXT "-CHEMAPP" )
	ADD_DEFINITIONS(-DCHEMAPP)
ENDIF(OGS_FEM_CHEMAPP)

# Add subdirectories with the projects
ADD_SUBDIRECTORY( Base )
ADD_SUBDIRECTORY( MathLib )
ADD_SUBDIRECTORY( GEO )
ADD_SUBDIRECTORY( MSH )
ADD_SUBDIRECTORY( FEM )
IF(OGS_FEM_GEMS)
	ADD_SUBDIRECTORY( GEM )
ENDIF(OGS_FEM_GEMS)
IF(OGS_FEM_CHEMAPP)
	ADD_SUBDIRECTORY( EQL )
	LINK_DIRECTORIES( ${CMAKE_SOURCE_DIR}/EQL )
ENDIF(OGS_FEM_CHEMAPP)
ADD_SUBDIRECTORY( FileIO )
ADD_SUBDIRECTORY( OGSProject )
ADD_SUBDIRECTORY( OGS )

IF(OGS_USE_QT)
   MESSAGE (STATUS  "Configuring for OGS with GUI" )	
   SET(OGS_FEM ON )
	IF(OGS_USE_OPENSG AND CMAKE_COMPILER_IS_GNUCXX)
		ADD_DEFINITIONS (-D__STDC_FORMAT_MACROS)	# for OpenSG compiling
	ENDIF()
   ADD_SUBDIRECTORY( Qt )
ENDIF(OGS_USE_QT)

# Project name
IF(OGS_USE_QT)
	PROJECT( OGS-${OGS_VERSION_MAJOR}-GUI )
ELSE(OGS_USE_QT)
	PROJECT( OGS-FEM-${OGS_VERSION_MAJOR}${PRJ_EXT} )
ENDIF(OGS_USE_QT)

# Create Configure.h
CONFIGURE_FILE (Base/Configure.h.in ${PROJECT_BINARY_DIR}/Base/Configure.h)
CONFIGURE_FILE (Base/BuildInfo.h.in ${PROJECT_BINARY_DIR}/Base/BuildInfo.h)

IF (MSVC)
	ADD_SUBDIRECTORY (CMakeConfiguration)
ENDIF (MSVC)

IF (BENCHMARK_DIR_FOUND OR GTEST_FOUND)
	ENABLE_TESTING()
ENDIF (BENCHMARK_DIR_FOUND OR GTEST_FOUND)

IF (EXAMPLEDATA_DIR_FOUND)
	INSTALL (DIRECTORY ${EXAMPLEDATA_DIR_FOUND} DESTINATION . PATTERN .svn EXCLUDE)
ENDIF (EXAMPLEDATA_DIR_FOUND)

IF (GTEST_FOUND)
	include_directories(${GTEST_INCLUDE_DIRS})
	add_subdirectory(tests)
ENDIF (GTEST_FOUND)

IF (OGS_PACKAGING)
	# Add packaging configuration
	MESSAGE (STATUS "Enabling packaging")
	INCLUDE( CMakeConfiguration/Pack.cmake)
ENDIF (OGS_PACKAGING)

IF (OGS_CMAKE_DEBUG)
	INCLUDE (CMakeConfiguration/PrintVariables.cmake)
ENDIF (OGS_CMAKE_DEBUG)	

IF (BENCHMARK_DIR_FOUND)
	STRING (COMPARE EQUAL "$ENV{HOSTNAME}" "dev2.intern.ufz.de" HOST_IS_DEV2)
	MESSAGE (STATUS "Benchmarking through CTest enabled")
	IF (${HOST_IS_DEV2})
		MESSAGE (STATUS "On dev2: file comparing enabled")
	ENDIF (${HOST_IS_DEV2})
	
	FILE (GLOB BENCHMARK_CONFIGS "${PROJECT_SOURCE_DIR}/../benchmarks/*.cmake")
	FOREACH (BENCHMARK_CONFIG ${BENCHMARK_CONFIGS})
		INCLUDE ("${BENCHMARK_CONFIG}")
	ENDFOREACH(BENCHMARK_CONFIG)
	
	UNSET (COPY_BENCHMARKS_TO_REF CACHE)
ENDIF (BENCHMARK_DIR_FOUND)

IF (OGS_BUILD_UTILITIES)
	ADD_SUBDIRECTORY (UTL/MSHGEOTOOLS/)
	ADD_SUBDIRECTORY (UTL/FileConverter/)
	ADD_SUBDIRECTORY (tests/GEO/)
ENDIF (OGS_BUILD_UTILITIES)
