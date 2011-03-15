## Programmed by WW
##           20.04.2010.  WW	


# Project name
IF(NOT OGS_USE_QT)
	PROJECT( OGS-FEM-${OGS_VERSION_MAJOR}${PRJ_EXT} )
ENDIF(NOT OGS_USE_QT)

IF(MKL)
	IF (UNIX)
		# Find MKLlib
		FIND_PACKAGE( MKL REQUIRED )
		SET(PARALLEL_USE_OPENMP ON)

		INCLUDE_DIRECTORIES (${MKL_INCLUDE_DIR})
	ELSE(UNIX)
		MESSAGE (FATAL_ERROR "MKL is only supported under LINUX/UNIX" )	
	ENDIF (UNIX)
ENDIF(MKL)

IF(LIS)
	# Find LISlib
	FIND_PACKAGE( LIS REQUIRED )
	#set (LIS ON)
	set (NEW_EQS ON)
	add_definitions(
		-o3
		-DIPMGEMPLUGIN
	)	
#	IF (UNIX)
#		FIND_PACKAGE( LIS REQUIRED )
#		#set (LIS ON)
#		set (NEW_EQS ON)
#		add_definitions(
#			-o3
#			-DIPMGEMPLUGIN
#		)	
#	ELSE(UNIX)
#		MESSAGE (FATAL_ERROR  "LIS is only supported under LINUX/UNIX" )	
#	ENDIF (UNIX)
ENDIF(LIS)

# Find OpenMP
IF(PARALLEL_USE_OPENMP)
	FIND_PACKAGE( OpenMP )
 
## Might be needed 
#	IF (MSVC)
#		SET(CMAKE_CXX_FLAGS_DEBUG "USE_OPENMP" /openmp)
#		SET(CMAKE_CXX_FLAGS_RELEASE "USE_OPENMP" /openmp)
#	ENDIF(MSVC)

#	IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
#		add_definitions(
#			-O3 
#			-fopenmp 
#			-lpthread 
#		)
#	ENDIF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)

	IF(OPENMP_FOUND)
		SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )	
	ENDIF(OPENMP_FOUND)
ENDIF(PARALLEL_USE_OPENMP)

IF(PARALLEL_USE_MPI)
	IF (WIN32)
		MESSAGE (FATAL_ERROR "Aborting: MPI is only supported under UNIX/LINUX!")
#		FIND_PACKAGE(MPI)
#		IF(MPI_FOUND)		
##			SET(CMAKE_C_COMPILER mpicc)
##			SET(CMAKE_CXX_COMPILER mpicxx)
#		ELSE(MPI_FOUND)
#			MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
#		ENDIF(MPI_FOUND)			
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
				ADD_DEFINITIONS(-DNEW_EQS -DUSE_MPI)    		
			ENDIF(MPI_COMPILER)
		ELSE( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
		FIND_PACKAGE(MPI)
		IF(MPI_FOUND)		
			SET(CMAKE_C_COMPILER ${MPI_COMPILER})
			SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
			ADD_DEFINITIONS(-DNEW_EQS -DUSE_MPI)    			
		ELSE(MPI_FOUND)
			MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
		ENDIF(MPI_FOUND)			
		ENDIF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
	ENDIF(UNIX)
ENDIF(PARALLEL_USE_MPI)


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
