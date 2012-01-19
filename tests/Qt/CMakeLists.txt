ADD_SUBDIRECTORY( gdiff )

INCLUDE( ${VTK_USE_FILE} )

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/Base
  ${CMAKE_SOURCE_DIR}/GEO
  ${CMAKE_SOURCE_DIR}/FEM
  ${CMAKE_SOURCE_DIR}/FileIO
  ${CMAKE_SOURCE_DIR}/MathLib
  ${CMAKE_SOURCE_DIR}/MSH
  ${CMAKE_SOURCE_DIR}/Qt/QtDataView
  ${CMAKE_SOURCE_DIR}/Qt/VtkVis
  ${CMAKE_SOURCE_DIR}/tests/Qt/gdiff
  ${CMAKE_SOURCE_DIR}/tests/Qt
  ${CMAKE_BINARY_DIR}/tests/Qt
)

SET(MOC_HEADERS
  QtTestBase.h
  
  # Add your test file here:
  Vtk/TestMeshFromRaster.h
)


qt4_wrap_cpp( MOC_SOURCES ${MOC_HEADERS} )

ADD_EXECUTABLE(QtTestrunner ${MOC_HEADERS} ${MOC_SOURCES} QtTestBase.cpp QtTestrunner.cpp)

SET(QTTESTRUNNER_LIBS
  Base
  GEO
  FEM
  FileIO
  MathLib
  MSH
  ${QT_LIBRARIES}
  ${libgeotiff_LIBRARIES}
  ${libtiff_LIBRARIES}
  VtkVis
  QtDataView
  GDiff
)

TARGET_LINK_LIBRARIES(QtTestrunner ${QTTESTRUNNER_LIBS})

# Tests can be run with make qttests
ADD_CUSTOM_TARGET(qtests QtTestrunner DEPENDS  ${QTTESTRUNNER_LIBS})

# Prepares the results directory for test diffs
ADD_CUSTOM_COMMAND(TARGET QtTestrunner PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E remove_directory results
	COMMAND ${CMAKE_COMMAND} -E make_directory results
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests COMMENT "Cleaning diff files")

SET_TARGET_PROPERTIES(QtTestrunner qtests PROPERTIES FOLDER Testing)

ADD_TEST(QTests ${QtTestrunnerExe} -xunitxml -o ./QtTestrunner.xml)
