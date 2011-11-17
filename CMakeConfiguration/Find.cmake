############################
### Find OGS directories ###
############################

IF(DEFINED BENCHMARK_DIR)
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${BENCHMARK_DIR})
ELSE()
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${PROJECT_SOURCE_DIR}/../benchmarks)
ENDIF()

IF(DEFINED EXAMPLEDATA_DIR)
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${EXAMPLEDATA_DIR})
ELSE()
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${PROJECT_SOURCE_DIR}/../ExampleData)
ENDIF()

FIND_PATH (OGS_LIBS_DIR_FOUND geotiff.lib ${PROJECT_SOURCE_DIR}/../Libs/libgeotiff)

# Find precompiled libraries (for BRNS GEMS LIS)
FIND_PATH (OGS_PRECOMPILED_LIBS_DIR_FOUND GEMS3_rl.lib ${PROJECT_SOURCE_DIR}/../Libs/precompiled)
IF (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
	LINK_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
ELSE (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	IF (WIN32)
		IF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also checked out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
	ELSE (WIN32)
		IF (OGS_FEM_LIS)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also checked out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_LIS)	
	ENDIF (WIN32)
ENDIF (OGS_PRECOMPILED_LIBS_DIR_FOUND)

FIND_PATH (LIB_DIR_FOUND GEMS3_rl.lib ${PROJECT_SOURCE_DIR}/LIB)
IF (LIB_DIR_FOUND)
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/LIB )
	LINK_DIRECTORIES( ${PROJECT_SOURCE_DIR}/LIB )
ENDIF (LIB_DIR_FOUND)

######################
### Find libraries ###
######################
FIND_PACKAGE (PythonInterp)
FIND_PACKAGE( Shapelib )
IF(Shapelib_FOUND)
	ADD_DEFINITIONS(-DShapelib_FOUND)
ENDIF() # Shapelib_FOUND

## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
	SET (HAVE_PTHREADS TRUE)
	MESSAGE (STATUS "pthread library found." )
ADD_DEFINITIONS(-DHAVE_PTHREADS)
ENDIF (CMAKE_USE_PTHREADS_INIT )

## boost (see FindBoost.cmake for more options) ##
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
FIND_PACKAGE( Boost COMPONENTS filesystem system regex)

## VTK ##
IF (OGS_LIBS_DIR_FOUND)
	SET (VTK_DIR ${PROJECT_SOURCE_DIR}/../Libs/VTK/build)
ENDIF () # OGS_LIBS_DIR_FOUND
FIND_PACKAGE( VTK )

IF (OGS_PYTHON)
	FIND_PACKAGE (PythonLibs 2.5 REQUIRED)
ENDIF (OGS_PYTHON)

## Qt4 library ##
IF(NOT OGS_DONT_USE_QT)
	FIND_PACKAGE( Qt4 4.5)
ENDIF(NOT OGS_DONT_USE_QT)

IF ( QT4_FOUND )
	# OPTION(OGS_GUI OFF )
	# this is needed to correctly link the qt libraries through target_link_libraries
	# By default only QtCore and QtGui modules are enabled
	# other modules must be enabled like this:
	SET( QT_USE_QTOPENGL TRUE )
	SET( QT_USE_QTSQL TRUE )
	SET( QT_USE_QTTEST TRUE )
	SET( QT_USE_QTXML TRUE )
	IF (QT_QTXMLPATTERNS_FOUND)
		set( QT_USE_QTXMLPATTERNS TRUE )
	ENDIF (QT_QTXMLPATTERNS_FOUND)
	INCLUDE( ${QT_USE_FILE} )
	IF (OGS_USE_QT)
		SET(OGS_USE_NETCDF ON FORCE)
	ENDIF (OGS_USE_QT)
ENDIF (QT4_FOUND )

## VRPN ##
FIND_PACKAGE( VRPN )

## VRED ##
FIND_PATH (VRED_DIR_FOUND vrNodePtr.h ${VRED_DIR}/include/vred)

## OpenSG ##
FIND_PACKAGE( OpenSG COMPONENTS OSGBase OSGSystem)

## gtest ##
IF(CMAKE_MINOR_VERSION GREATER 6)
	IF (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" AND NOT OGS_FEM_GEMS )
		#SET (GTEST_MSVC_SEARCH MT)
		SET (GTEST_ROOT ${PROJECT_SOURCE_DIR}/../Libs/gtest)
	ENDIF (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" AND NOT OGS_FEM_GEMS )
	FIND_PACKAGE (GTest)
ENDIF(CMAKE_MINOR_VERSION GREATER 6)


IF(MKL)
	# Find MKLlib
	FIND_PACKAGE( MKL REQUIRED )
	INCLUDE_DIRECTORIES (${MKL_INCLUDE_DIR})
ENDIF(MKL)

IF(LIS)
	# Find LISlib
	FIND_PACKAGE( LIS REQUIRED )
	set (NEW_EQS ON)
	add_definitions(
		-o3
		-DIPMGEMPLUGIN
	)
ENDIF(LIS)

# Find OpenMP
IF(PARALLEL_USE_OPENMP)
	FIND_PACKAGE( OpenMP REQUIRED )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
ENDIF(PARALLEL_USE_OPENMP)

IF(PARALLEL_USE_MPI)
	IF (WIN32)
#		MESSAGE (FATAL_ERROR "Aborting: MPI is only supported under UNIX/LINUX!")			
		#ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
		FIND_PACKAGE(MPI REQUIRED)
	ENDIF(WIN32)
	IF(UNIX)

# If there is an mpi compiler find it and interogate (farther below) it for the include
# and lib dirs otherwise we will continue to search from ${_MPI_BASE_DIR}.

		IF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
			find_program(MPI_COMPILER
 			 NAMES mpic++ mpicxx mpiCC mpicc
			 HINTS "${_MPI_BASE_DIR}"
 			 PATH_SUFFIXES bin
			  DOC "MPI compiler. Used only to detect MPI compilation flags.")
			IF(MPI_COMPILER)

			MESSAGE (STATUS  "CMake version is less than 2.8, MPI compiler is set directly" )	
			mark_as_advanced(MPI_COMPILER)
				SET(CMAKE_C_COMPILER ${MPI_COMPILER})
				SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})   		
			ENDIF(MPI_COMPILER)
		ELSE( ${CMAKE_MAJOR_VERSION}  EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
			FIND_PACKAGE(MPI REQUIRED)		
		ENDIF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
	ENDIF(UNIX)
ENDIF(PARALLEL_USE_MPI)

######################
### Find tools     ###
######################

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")

FIND_PACKAGE(cppcheck)
