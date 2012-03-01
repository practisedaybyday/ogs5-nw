# Source files
SET( SOURCES
	GMSInterface.cpp
	GocadInterface.cpp
	OGSIOVer4.cpp
	PetrelInterface.cpp	
	StationIO.cpp
	readNonBlankLineFromInputStream.cpp
	#FEMIO.cpp
	FEMIO/GeoIO.cpp
	FEMIO/ProcessIO.cpp
	FEMIO/BoundaryConditionIO.cpp
	MeshIO/GMSHInterface.cpp
	MeshIO/LegacyVtkInterface.cpp
	MeshIO/OGSMeshIO.cpp
	MeshIO/TetGenInterface.cpp
	GeoIO/Gmsh2GeoIO.cpp
)

# Header files
SET( HEADERS
	GMSInterface.h
	GocadInterface.h
	OGSIOVer4.h
	PetrelInterface.h
	StationIO.h
	readNonBlankLineFromInputStream.h
	#FEMIO.h
	FEMIO/GeoIO.h
	FEMIO/ProcessIO.h
	FEMIO/BoundaryConditionIO.h
	MathIO/CRSIO.h
	MeshIO/GMSHInterface.h
	MeshIO/LegacyVtkInterface.h
	MeshIO/OGSMeshIO.h
	MeshIO/TetGenInterface.h
	GeoIO/Gmsh2GeoIO.h
)

IF (OGS_USE_NETCDF)
	INCLUDE( ${VTK_USE_FILE} )
	SET (HEADERS NetCDFInterface.h ${HEADERS})
	SET (SOURCES NetCDFInterface.cpp ${SOURCES})
ENDIF (OGS_USE_NETCDF)

IF (Shapelib_FOUND)
	SET( HEADERS ${HEADERS} SHPInterface.h )
	SET (SOURCES ${SOURCES} SHPInterface.cpp )
	INCLUDE_DIRECTORIES (${Shapelib_INCLUDE_DIR})
ENDIF (Shapelib_FOUND)

IF (QT4_FOUND)
	SET( HEADERS ${HEADERS} XmlIO/XMLInterface.h )
	SET (SOURCES ${SOURCES} XmlIO/XMLInterface.cpp )
	SET( HEADERS ${HEADERS} XmlIO/XmlGmlInterface.h )
	SET (SOURCES ${SOURCES} XmlIO/XmlGmlInterface.cpp )
	SET( HEADERS ${HEADERS} XmlIO/XmlStnInterface.h )
	SET (SOURCES ${SOURCES} XmlIO/XmlStnInterface.cpp )
	SET( HEADERS ${HEADERS} XmlIO/XmlCndInterface.h )
	SET (SOURCES ${SOURCES} XmlIO/XmlCndInterface.cpp )
	SET( HEADERS ${HEADERS} XmlIO/XmlGspInterface.h )
	SET (SOURCES ${SOURCES} XmlIO/XmlGspInterface.cpp )
	SET( HEADERS ${HEADERS} FEFLOWInterface.h )
	SET (SOURCES ${SOURCES} FEFLOWInterface.cpp )
ENDIF (QT4_FOUND)

IF(OGS_FEM_GEMS)
	INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF()

# Create the library
ADD_LIBRARY( FileIO STATIC
	${SOURCES}
	${HEADERS}
)

INCLUDE_DIRECTORIES(
	.
	./FEMIO
	./MathIO
	./MeshIO
	./GeoIO
	../Base
	../FEM
	../GEO
	../MathLib
	../MSH
	../GEM
)

TARGET_LINK_LIBRARIES( FileIO
#	Base
#	MathLib
	GEO
	FEM
)

IF (Shapelib_FOUND)
	TARGET_LINK_LIBRARIES( FileIO ${Shapelib_LIBRARIES} )
ENDIF (Shapelib_FOUND)

IF(OGS_PACKAGING)
	# Install xml schema files
	FILE (GLOB XML_SCHEMA_FILES . "*.xs?")
	INSTALL (FILES ${XML_SCHEMA_FILES} DESTINATION bin)
ENDIF(OGS_PACKAGING)
