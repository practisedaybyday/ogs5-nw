
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
	${CMAKE_SOURCE_DIR}/Qt/VtkVis
)

# Create executables
ADD_EXECUTABLE( ExtractMeshNodes mainExtractMeshNodes.cpp )
ADD_EXECUTABLE( ModifyMeshProperties mainModifyMeshProperties.cpp ModifyMeshProperties.h ModifyMeshProperties.cpp )
ADD_EXECUTABLE( filterMeshNodes filterMeshNodes.cpp )
ADD_EXECUTABLE( rotateMesh rotateMesh.cpp )

SET_TARGET_PROPERTIES(ExtractMeshNodes ModifyMeshProperties filterMeshNodes rotateMesh
  PROPERTIES FOLDER Utilities)

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


IF(QT4_FOUND)

	ADD_EXECUTABLE( extractSurface mainExtractSurface.cpp ExtractSurface.h ExtractSurface.cpp )
	ADD_EXECUTABLE( createMeshFromPolygon createMeshFromPolygon.cpp )

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

	ADD_EXECUTABLE( mapGeoToDEM mapGeoToDEM.cpp )
	TARGET_LINK_LIBRARIES( mapGeoToDEM
		GEO
		MSH
		FEM
		FileIO
		VtkVis
	)
	IF(libgeotiff_FOUND)
		TARGET_LINK_LIBRARIES( mapGeoToDEM ${libgeotiff_LIBRARIES})
	ENDIF() # libgeotiff_FOUND
	SET_TARGET_PROPERTIES(mapGeoToDEM PROPERTIES FOLDER Utilities)

	TARGET_LINK_LIBRARIES( rotateMesh
		GEO
		MSH
		MSHGEOTOOLS
		FileIO
		${QT_LIBRARIES}
	)	

	ADD_EXECUTABLE( testMeshSearchAlgorithms testMeshSearchAlgorithms.cpp )
	TARGET_LINK_LIBRARIES( testMeshSearchAlgorithms
		GEO
		MSH
		MSHGEOTOOLS
		FileIO
		${QT_LIBRARIES}
	)	
	
	ADD_EXECUTABLE( removeMeshNodes removeMeshNodes.cpp )
	TARGET_LINK_LIBRARIES( removeMeshNodes
		Base
		FileIO
		MSH
		${QT_LIBRARIES}
	)

	SET_TARGET_PROPERTIES(extractSurface createMeshFromPolygon removeMeshNodes
		PROPERTIES FOLDER Utilities)

ENDIF() # QT4_FOUND

IF (OGS_PACKAGING)
	INSTALL ( TARGETS ExtractMeshNodes RUNTIME DESTINATION bin COMPONENT
ExtractMeshNodes )
ENDIF (OGS_PACKAGING)
