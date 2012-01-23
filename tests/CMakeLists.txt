cmake_minimum_required(VERSION 2.6)

ADD_SUBDIRECTORY( gtest )
ADD_SUBDIRECTORY( gdiff )

SET ( SOURCES
	testrunner.cpp
	testBase.cpp
	testGeo.cpp
)

IF (NOT MSVC)
	SET ( SOURCES ${SOURCES}
		testOgs.cxx
		testBmSkel.cxx)
ENDIF (NOT MSVC)

include_directories(
	gtest
    ../Base
    ../FileIO
	../GEO
	../MathLib
    ../FEM
    ../MSH
    ../OGSProject
)

add_executable (testrunner ${SOURCES})
SET(TESTRUNNER_LIBS 
  GTest
  Base
  FileIO
  GEO
  MathLib
  FEM
  MSH
  OGSProject
)
target_link_libraries(testrunner ${TESTRUNNER_LIBS})

IF (MSVC)
	SET (testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/testrunner)
	SET (QtTestrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/QtTestrunner)
ELSE ()
	SET (testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/testrunner)
	IF ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
		TARGET_LINK_LIBRARIES( testrunner pthread )
	ENDIF ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
	SET (QtTestrunnerExe ${EXECUTABLE_OUTPUT_PATH}/QtTestrunner)
ENDIF()

ADD_CUSTOM_TARGET(gtests testrunner	DEPENDS ${TESTRUNNER_LIBS})

SET_TARGET_PROPERTIES(testrunner gtests PROPERTIES FOLDER Testing)

macro(ADD_GOOGLE_TESTS executable)
	foreach ( source ${ARGN} )
		file(READ "${source}" contents)
		string(REGEX MATCHALL "TEST_?F?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
		foreach(hit ${found_tests})
			string(REGEX REPLACE ".*\\(([A-Za-z_0-9]+)[, ]*([A-Za-z_0-9]+)\\).*" "\\1.\\2" test_name ${hit})
			add_test(${test_name} ${executable}  --gtest_output=xml --gtest_filter=${test_name} ${MI3CTestingDir})
			# message ("Adding test: ${test_name}")
		endforeach(hit)
	endforeach()
endmacro()

# Either collect all gtests in one ctest run of testrunner executable
# or run each gtest individually
IF (OGS_ONE_BIG_GTEST)
	add_test (GTests ${testrunnerExe}  --gtest_output=xml:./)
ELSE (OGS_ONE_BIG_GTEST)
	ADD_GOOGLE_TESTS ( ${testrunnerExe} ${SOURCES})
ENDIF (OGS_ONE_BIG_GTEST)


# copies the model files to the build dir and adds them as targets so that
# the build files are re-built when the source model files change
macro ( UPDATE_MODEL_FILES dirOUT fileLIST )
	get_filename_component( _tdir ${CMAKE_CURRENT_SOURCE_DIR} NAME )
	#message (STATUS "Copying files to ${dirOUT} from ${fileLIST}.\n")
	foreach ( _file1 ${${fileLIST}} )
		set( _file ${CMAKE_CURRENT_SOURCE_DIR}/${_file1} )
		get_filename_component( _fdest ${_file} NAME )
		set( dest ${dirOUT}/${_fdest} )
		#message( STATUS "Copying ${_file} to ${dest} \n" )
  
		#message (STATUS "Adding targets ${_tdir}.${_fdest}\n").\n")
		add_custom_target( ${_tdir}.${_fdest}
			${CMAKE_COMMAND} -E copy_if_different
			${_file} ${dest})
  
		add_dependencies( testrunner ${_tdir}.${_fdest} )
	endforeach(_file1 ${${fileLIST}})
  
endmacro ( UPDATE_MODEL_FILES dirOUT fileLIST)
  
# add new test models here not all of these models will be needed, depending
# on compile flags but it is safe to copy them all across
IF (NOT MSVC)
	add_subdirectory( data/bmskel )
	add_subdirectory( data/GEO )
ENDIF (NOT MSVC)

IF(OGS_USE_QT AND TESTDATA_DIR_FOUND)
  ADD_SUBDIRECTORY( Qt )
ENDIF()