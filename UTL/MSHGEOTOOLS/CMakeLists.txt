
INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/UTL/MSHGEOTOOLS
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSHGEOTOOLS
)

# Create executables
ADD_EXECUTABLE( ExtractMeshNodes mainExtractMeshNodes.cpp )
ADD_EXECUTABLE( ModifyMeshProperties mainModifyMeshProperties.cpp ModifyMeshProperties.h ModifyMeshProperties.cpp )
ADD_EXECUTABLE( filterMeshNodes filterMeshNodes.cpp )
ADD_EXECUTABLE( extractSurface mainExtractSurface.cpp ExtractSurface.h ExtractSurface.cpp )
ADD_EXECUTABLE( createMeshFromPolygon createMeshFromPolygon.cpp )
SET_TARGET_PROPERTIES(ExtractMeshNodes ModifyMeshProperties filterMeshNodes
  extractSurface createMeshFromPolygon PROPERTIES FOLDER Utilities)

TARGET_LINK_LIBRARIES( ExtractMeshNodes
	GEO
	MSH
	MSHGEOTOOLS
	FileIO
	FEM
	OGSProject
)

TARGET_LINK_LIBRARIES( ModifyMeshProperties
	GEO
	MSH
	MathLib
	FileIO
	FEM
	OGSProject
)

TARGET_LINK_LIBRARIES( filterMeshNodes
	Base
	MSH
)

TARGET_LINK_LIBRARIES( extractSurface
	GEO
	MSH
	FileIO
	FEM
	OGSProject
	${QT_LIBRARIES}
)

TARGET_LINK_LIBRARIES( createMeshFromPolygon
	GEO
	MSH
	MSHGEOTOOLS
	FileIO
	${QT_LIBRARIES}
)

IF (OGS_PACKAGING)
	INSTALL ( TARGETS ExtractMeshNodes RUNTIME DESTINATION bin COMPONENT
ExtractMeshNodes )
ENDIF (OGS_PACKAGING)
