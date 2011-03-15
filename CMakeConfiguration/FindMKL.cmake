if (NOT MKL_FOUND)

	#Extract MKL
	FIND_PATH (MKL_INCLUDE_DIR_FOUND mkl.h ${PROJECT_SOURCE_DIR}/../Libs/MKL/include)
	IF (NOT MKL_INCLUDE_DIR_FOUND)
		FIND_PATH (MKL_DIR_FOUND mkl-include.tgz ${PROJECT_SOURCE_DIR}/../Libs/MKL)
		IF (MKL_DIR_FOUND)
			MESSAGE (STATUS "Uncompressing MKL...")
			EXECUTE_PROCESS (COMMAND tar xvzf mkl-include.tgz WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../Libs/MKL/ OUTPUT_QUIET)
			IF (HAVE_64_BIT)
				EXECUTE_PROCESS (COMMAND tar xvzf mkl-64.tgz WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../Libs/MKL/ OUTPUT_QUIET)
			ELSE (HAVE_64_BIT)
				EXECUTE_PROCESS (COMMAND tar xvzf mkl-32.tgz WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../Libs/MKL/ OUTPUT_QUIET)
			ENDIF (HAVE_64_BIT)
		ELSE (MKL_DIR_FOUND)
			MESSAGE (STATUS "MKL archives in ../Libs/ not found")
		ENDIF (MKL_DIR_FOUND)
	ENDIF (NOT MKL_INCLUDE_DIR_FOUND)

	include(LibFindMacros)
	
	find_path( MKL_INCLUDE_DIR NAMES mkl.h
		   PATHS ${CMAKE_SOURCE_DIR}/../Libs/MKL/include)

	if ( UNIX )
		# Tell if the unix system is on 64-bit base
		if(CMAKE_SIZEOF_VOID_P MATCHES "8")
		set (MKL_LIB_PATH "${CMAKE_SOURCE_DIR}/../Libs/MKL/64")
			find_library(MKL_SOLVER_LIBRARY
			    NAMES mkl_solver_lp64 PATHS ${MKL_LIB_PATH} )
			find_library(MKL_INTEL_LIBRARY
			    NAMES mkl_intel_lp64 PATHS ${MKL_LIB_PATH} )
			find_library(MKL_CORE_LIBRARY
			    NAMES mkl_core PATHS ${MKL_LIB_PATH} )
			set(MKL_LIBRARIES ${MKL_SOLVER_LIBRARY} ${MKL_INTEL_LIBRARY} ${MKL_CORE_LIBRARY} CACHE STRING "Found MKL Libraries")
		    if(CMAKE_C_COMPILER MATCHES "icc")
		        find_library(MKL_INTEL_THREAD_LIBRARY
				    NAMES mkl_intel_thread PATHS ${MKL_LIB_PATH} )
				find_library(MKL_IOMP5_LIBRARY
				    NAMES iomp5 PATHS ${MKL_LIB_PATH} )
				find_library(MKL_PTHREAD_LIBRARY
				    NAMES pthread PATHS ${MKL_LIB_PATH} )
				set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_IOMP5_LIBRARY} ${MKL_PTHREAD_LIBRARY})
		    else(CMAKE_C_COMPILER MATCHES "icc")
				find_library(MKL_GNU_THREAD_LIBRARY
				    NAMES mkl_gnu_thread PATHS ${MKL_LIB_PATH} )
				set(MKL_LIBRARIES ${MKL_LIBRARIES} ${MKL_GNU_THREAD_LIBRARY})
		    endif(CMAKE_C_COMPILER MATCHES "icc") 
		
		else (CMAKE_SIZEOF_VOID_P MATCHES "8")
			set (MKL_LIB_PATH "${CMAKE_SOURCE_DIR}/../Libs/MKL/32")
			find_library(MKL_SOLVER_LIBRARY
				NAMES mkl_solver PATHS ${MKL_LIB_PATH} )
			find_library(MKL_INTEL_LIBRARY
				NAMES mkl_intel PATHS ${MKL_LIB_PATH} )
			find_library(MKL_GNU_THREAD_LIBRARY
				NAMES mkl_gnu_thread PATHS ${MKL_LIB_PATH} )	
			find_library(MKL_CORE_LIBRARY
				NAMES mkl_core PATHS ${MKL_LIB_PATH} )
			set(MKL_LIBRARIES ${MKL_SOLVER_LIBRARY} ${MKL_INTEL_LIBRARY} ${MKL_GNU_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} CACHE STRING "Found MKL Libraries")
		endif (CMAKE_SIZEOF_VOID_P MATCHES "8")	
	endif ( UNIX )

	set(MKL_LIBRARIES ${MKL_SOLVER_LIBRARY} ${MKL_INTEL_LIBRARY} ${MKL_GNU_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} CACHE STRING "Found MKL Libraries")

	# Set the include dir variables and the libraries and let libfind_process do the rest.
	# NOTE: Singular variables for this library, plural for libraries this lib depends on.
	set(MKL_PROCESS_INCLUDES MKL_INCLUDE_DIR)
	set(MKL_PROCESS_LIBS MKL_LIBRARIES)
	libfind_process(MKL)
	
endif (NOT MKL_FOUND)
