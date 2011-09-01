# Source files
set( SOURCES
	AxisAlignedBoundingBox.cpp
	BruteForceClosestPair.cpp
	Color.cpp
#	DivideAndConquerClosestPair.cpp
	geo_dom.cpp
	geo_lib.cpp
	geo_lin.cpp
	geo_mathlib.cpp
	geo_ply.cpp
	geo_pnt.cpp
	geo_sfc.cpp
	geo_vol.cpp
	GEOObjects.cpp
	GeoType.cpp
	Point.cpp
	PointVec.cpp
	Polygon.cpp
	Polyline.cpp
	PolylineSmoother.cpp
	ProjectData.cpp
	SimplePolygonTree.cpp
	Station.cpp
	Surface.cpp
	Triangle.cpp
)

# Header files
set( HEADERS
	AxisAlignedBoundingBox.h
	BruteForceClosestPair.h
	ClosestPair.h
	Color.h
#	DivideAndConquerClosestPair.h
	geo_dom.h
	geo_lib.h
	geo_lin.h
	geo_mathlib.h
	geo_ply.h
	geo_pnt.h
	geo_sfc.h	
	geo_vol.h
	GeoObject.h
	GEOObjects.h
	GeoType.h
	Point.h
	PointVec.h
	PointWithID.h
	Polygon.h
	Polyline.h
	PolylineSmoother.h
	PolylineVec.h
	ProjectData.h
	PropertyBounds.h
	QuadTree.h
	SimplePolygonTree.h
	Station.h
	Surface.h
	SurfaceVec.h
	TemplatePoint.h
	TemplateVec.h
	Triangle.h
)

# Create the library
add_library( GEO STATIC
	${SOURCES}
	${HEADERS}
)


include_directories(
	.
	../Base
	../FEM
	../FileIO
	../MathLib
	../MSH
)


target_link_libraries (
	GEO
	Base
	FEM
	MathLib
)

