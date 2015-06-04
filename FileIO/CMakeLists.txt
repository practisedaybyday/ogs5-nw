# Source files
SET( SOURCES
	Writer.cpp
	OGSIOVer4.cpp
	StationIO.cpp
	readNonBlankLineFromInputStream.cpp
	#FEMIO.cpp
	FEMIO/GeoIO.cpp
	FEMIO/ProcessIO.cpp
	FEMIO/BoundaryConditionIO.cpp
	MeshIO/LegacyVtkInterface.cpp
	MeshIO/OGSMeshIO.cpp
	XmlIO/RapidXMLInterface.cpp
)

# Header files
SET( HEADERS
	Writer.h
	OGSIOVer4.h
	StationIO.h
	readNonBlankLineFromInputStream.h
	#FEMIO.h
	FEMIO/GeoIO.h
	FEMIO/ProcessIO.h
	FEMIO/BoundaryConditionIO.h
	MathIO/CRSIO.h
	MeshIO/LegacyVtkInterface.h
	MeshIO/OGSMeshIO.h
	XmlIO/RapidXMLInterface.h
)

IF(OGS_FEM_GEMS OR OGS_FEM_PETSC_GEMS)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

# Create the library
ADD_LIBRARY( FileIO STATIC
	${SOURCES}
	${HEADERS}
)

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/FileIO/MathIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/FileIO/GeoIO
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/GEM
)

TARGET_LINK_LIBRARIES( FileIO
	Base
	#	MathLib
	GEO
	FEM
)

IF (Shapelib_FOUND)
	TARGET_LINK_LIBRARIES( FileIO ${Shapelib_LIBRARIES} )
ENDIF (Shapelib_FOUND)
