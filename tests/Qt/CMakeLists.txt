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

TARGET_LINK_LIBRARIES(QtTestrunner
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

# Create directory for result files in build directory
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/results)