
INCLUDE_DIRECTORIES(
	./
	../../Base
	../../GEO
	../../FileIO
	../../MathLib
	../../FEM
	../../MSH
)

# Create executables
ADD_EXECUTABLE( ExtractMeshNodes mainExtractMeshNodes.cpp ExtractMeshNodes.h ExtractMeshNodes.cpp )
ADD_EXECUTABLE( ModifyMeshProperties mainModifyMeshProperties.cpp ModifyMeshProperties.h ModifyMeshProperties.cpp )
ADD_EXECUTABLE( filterMeshNodes filterMeshNodes.cpp )

TARGET_LINK_LIBRARIES( ExtractMeshNodes
	GEO
	MSH
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
)

IF (OGS_PACKAGING)
	INSTALL ( TARGETS ExtractMeshNodes RUNTIME DESTINATION bin COMPONENT
ExtractMeshNodes )
ENDIF (OGS_PACKAGING)
