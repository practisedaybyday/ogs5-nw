# - Enable Code Coverage
#
# 2012-01-31, Lars Bilke
#
# USAGE:
# 1. Copy this file into your cmake modules path
# 2. Add the following line to your CMakeLists.txt:
#      INCLUDE(CodeCoverage)
# 
# 3. Use the function SETUP_TARGET_FOR_COVERAGE to create a custom make target
#    which runs your test executable and produces a lcov code coverage report.
#

# Check prereqs
FIND_PROGRAM( GCOV_PATH gcov )
FIND_PROGRAM( LCOV_PATH lcov )
FIND_PROGRAM( GENHTML_PATH genhtml )

IF(NOT GCOV_PATH)
	MESSAGE(FATAL_ERROR "gcov not found! Aborting...")
ENDIF() # NOT GCOV_PATH

IF(NOT LCOV_PATH)
	MESSAGE(FATAL_ERROR "lcov not found! Aborting...")
ENDIF() # NOT LCOV_PATH

IF(NOT GENHTML_PATH)
	MESSAGE(FATAL_ERROR "genhtml not found! Aborting...")
ENDIF() # NOT GENHTML_PATH

IF(NOT CMAKE_COMPILER_IS_GNUCXX)
	MESSAGE(FATAL_ERROR "Compiler is not GNU gcc! Aborting...")
ENDIF() # NOT CMAKE_COMPILER_IS_GNUCXX

IF ( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
  MESSAGE( WARNING "Code coverage results with an optimised (non-Debug) build may be misleading" )
ENDIF() # NOT CMAKE_BUILD_TYPE STREQUAL "Debug"


# Setup compiler options
ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
LINK_LIBRARIES(gcov)


# Param _targetname     The name of new the custom make target
# Param _testrunnerexe  The executable name which runs the tests
# Param _outputname     lcov output is generated as _outputname.info
#                       HTML report is generated in _outputname/index.html
FUNCTION(SETUP_TARGET_FOR_COVERAGE _targetname _testrunnerexe _outputname)

	# Setup target
	ADD_CUSTOM_TARGET(${_targetname}
		
		# Cleanup lcov
		lcov --directory . --zerocounters
		
		# Run tests
		COMMAND ${_testrunnerexe}
		
		# Capturing lcov counters and generating report
		COMMAND lcov --directory . --capture --output-file ${_outputname}.info
		COMMAND lcov --remove ${_outputname}.info 'tests/*' '/usr/*' --output-file ${_outputname}.info.cleaned
		COMMAND genhtml -o ${_outputname} ${_outputname}.info.cleaned
		COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned
		
		DEPENDS ${_testrunnerexe}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		COMMENT "Resetting code coverage counters to zero.\n
			Processing code coverage counters and generating report."
	)
	
	# Show info where to find the report
	ADD_CUSTOM_COMMAND(TARGET ${_targetname} POST_BUILD
		COMMAND ;
		COMMENT "Open ./${_outputname}/index.html in your browser to view the coverage report."
	)

ENDFUNCTION() # SETUP_TARGET_FOR_COVERAGE
