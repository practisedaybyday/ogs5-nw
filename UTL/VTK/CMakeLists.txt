INCLUDE( ${VTK_USE_FILE} )

INCLUDE_DIRECTORIES(
	./
	../../Base
	../../GEO
	../../FileIO
	../../MathLib
	../../FEM
	../../MSH
	../../Qt/OpenSG
	../../Qt/VtkVis
)

# Create executables
ADD_EXECUTABLE( testMeshFromRaster testMeshFromRaster.cpp )

TARGET_LINK_LIBRARIES( testMeshFromRaster
	FEM
        MSH
       	GEO
        MathLib
       	FileIO
        Base
       	OGSProject
       	QtDataView
       	VtkVis
	${QT_LIBRARIES}
)
