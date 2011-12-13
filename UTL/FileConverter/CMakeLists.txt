
INCLUDE_DIRECTORIES(
	./
	../../Base
	../../GEO
	../../FileIO
	../../MathLib
	../../FEM
	../../MSH
	../../Qt/OpenSG
)

## Shape to GLI file converter ##
IF (Shapelib_FOUND)
	
	INCLUDE_DIRECTORIES (${Shapelib_INCLUDE_DIR})
	
	# Create executables
	ADD_EXECUTABLE( ConvertSHPToGLI ConvertSHPToGLI.cpp )

	TARGET_LINK_LIBRARIES( ConvertSHPToGLI
		FEM
	        MSH
	        MSHGEOTOOLS
        	GEO
	        MathLib
        	FileIO
	        Base
        	OGSProject
		${Shapelib_LIBRARIES}
		${QT_LIBRARIES}
	)

	IF (OGS_PACKAGING)
		INSTALL ( TARGETS ConvertSHPToGLI RUNTIME DESTINATION bin COMPONENT ConvertSHPToGLI )
	ENDIF (OGS_PACKAGING)
	
ENDIF (Shapelib_FOUND)

## VTK to OpenSG file converter ##
IF (OGS_USE_OPENSG AND VTK_FOUND AND Boost_FOUND)
	
	INCLUDE( ${VTK_USE_FILE} )
	
	ADD_EXECUTABLE (ConvertVtkToOsg ConvertVtkToOsg.cpp)
	
	USE_OPENSG ( ConvertVtkToOsg )
	INCLUDE_DIRECTORIES ( ../OpenSG ${Boost_INCLUDE_DIRS} )
	TARGET_LINK_LIBRARIES ( ConvertVtkToOsg OgsOpenSG vtkRendering ${Boost_LIBRARIES} )
	
	IF (OGS_PACKAGING)
		INSTALL ( TARGETS ConvertVtkToOsg RUNTIME DESTINATION bin COMPONENT ConvertVtkToOsg )
	ENDIF (OGS_PACKAGING)
	
ENDIF () # OGS_USE_OPENSG AND VTK_FOUND AND Boost_FOUND

## Gli to VTK file converter ##
IF(VTK_FOUND AND Boost_FOUND)
	
	INCLUDE( ${VTK_USE_FILE} )
	
	ADD_EXECUTABLE (ConvertGliToVtk ConvertGliToVtk.cpp)
	
	INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} ../../Qt/VtkVis )
	TARGET_LINK_LIBRARIES ( ConvertGliToVtk FileIO VtkVis vtkRendering ${Boost_LIBRARIES} )
	
	IF (OGS_PACKAGING)
		INSTALL ( TARGETS ConvertGliToVtk RUNTIME DESTINATION bin COMPONENT ConvertGliToVtk )
	ENDIF (OGS_PACKAGING)
	
ENDIF() # VTK_FOUND AND Boost_FOUND

IF (QT4_FOUND)
	ADD_EXECUTABLE (generateBCandGLI generateBCandGLI.cpp )
	TARGET_LINK_LIBRARIES (generateBCandGLI 
		FEM
		MSH
		MSHGEOTOOLS
		GEO
		MathLib
		FileIO
		Base
		OGSProject
		${QT_LIBRARIES}
	)

	ADD_EXECUTABLE (generateBCFromPolyline generateBCFromPolyline.cpp )
	TARGET_LINK_LIBRARIES (generateBCFromPolyline
		FEM
		MSH
		MSHGEOTOOLS
		GEO
		MathLib
		FileIO
		Base
		OGSProject
		${QT_LIBRARIES}
	)
ENDIF () # QT4_FOUND

ADD_EXECUTABLE (GMSH2OGS GMSH2OGS.cpp)
TARGET_LINK_LIBRARIES (GMSH2OGS
	MSH
	FileIO
)
