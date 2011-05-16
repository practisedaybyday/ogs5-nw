
SET ( SOURCES
	ExtractMeshNodes.h
	ExtractMeshNodes.cpp
)

# include directories
INCLUDE_DIRECTORIES(
	./
	../Base
	../GEO
	../MSH
	../MathLib
	../FileIO
	../FEM
)

# create library
ADD_LIBRARY( MSHGEOTOOLS STATIC ${SOURCES} )

TARGET_LINK_LIBRARIES (
        MSHGEOTOOLS
        MSH
        GEO
        ${CMAKE_THREAD_LIBS_INIT}
)

