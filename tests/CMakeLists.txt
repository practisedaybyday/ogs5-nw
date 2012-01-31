ADD_SUBDIRECTORY( gtest )

# Add tests which do not need testdata here:
SET ( SOURCES
	testrunner.cpp
	testBase.cpp
)

# Add tests here if they need testdata
IF(TESTDATA_DIR_FOUND)
	SET(SOURCES ${SOURCES}
		GEO/TestPointInTriangle.cpp
		GEO/TestPointInPolygon.cpp
	)
ENDIF()

IF (NOT MSVC)
	SET ( SOURCES ${SOURCES}
		testOgs.cxx
		testBmSkel.cxx)
ENDIF (NOT MSVC)

INCLUDE_DIRECTORIES(
	gtest
    ../Base
    ../FileIO
	../GEO
	../MathLib
    ../FEM
    ../MSH
    ../OGSProject
)

ADD_EXECUTABLE (testrunner ${SOURCES})
TARGET_LINK_LIBRARIES(testrunner 
  GTest
  Base
  FileIO
  GEO
  MathLib
  FEM
  MSH
  OGSProject
)

IF (MSVC)
	SET(testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/testrunner)
	SET(QtTestrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/QtTestrunner)
ELSE ()
	SET(testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/testrunner)
	IF("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
		TARGET_LINK_LIBRARIES( testrunner pthread )
	ENDIF() # "${CMAKE_SYSTEM_NAME}" MATCHES "Linux"
	SET (QtTestrunnerExe ${EXECUTABLE_OUTPUT_PATH}/QtTestrunner)
ENDIF() # MSVC

IF(OGS_COVERAGE)
	SETUP_TARGET_FOR_COVERAGE(gtests testrunner "gtests_coverage_results")
ELSE()
	ADD_CUSTOM_TARGET(gtests testrunner	DEPENDS testrunner)
ENDIF() # OGS_COVERAGE

SET_TARGET_PROPERTIES(testrunner gtests PROPERTIES FOLDER Testing)

# Either collect all gtests in one ctest run of testrunner executable
# or run each gtest individually
IF (OGS_ONE_BIG_GTEST)
	add_test (GTests ${testrunnerExe}  --gtest_output=xml:./)
ELSE (OGS_ONE_BIG_GTEST)
	ADD_GOOGLE_TESTS ( ${testrunnerExe} ${SOURCES})
ENDIF (OGS_ONE_BIG_GTEST)
  
# add new test models here not all of these models will be needed, depending
# on compile flags but it is safe to copy them all across
IF (NOT MSVC)
	ADD_SUBDIRECTORY( data/bmskel )
ENDIF (NOT MSVC)

IF(OGS_USE_QT AND TESTDATA_DIR_FOUND)
  ADD_SUBDIRECTORY( gdiff ) # depends on Qt
  ADD_SUBDIRECTORY( Qt )
ENDIF()