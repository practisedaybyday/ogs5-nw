# Source files
SET( SOURCES
	GMSHInterface.cpp
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
	MeshIO/OGSMeshIO.cpp
	MeshIO/LegacyVtkInterface.cpp
	GeoIO/Gmsh2GeoIO.cpp
)

# Header files
SET( HEADERS
	GMSHInterface.h
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
	MeshIO/OGSMeshIO.h
	MeshIO/LegacyVtkInterface.h
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
	SET( HEADERS ${HEADERS} XMLInterface.h )
	SET (SOURCES ${SOURCES} XMLInterface.cpp )
	SET( HEADERS ${HEADERS} FEFLOWInterface.h )
	SET (SOURCES ${SOURCES} FEFLOWInterface.cpp )
ENDIF (QT4_FOUND)

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
