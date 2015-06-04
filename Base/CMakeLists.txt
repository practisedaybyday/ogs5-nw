# Source files
set( HEADERS
	binarySearch.h
	Configure.h.in
	DateTools.h
	FileFinder.h
	FileTools.h
	MemWatch.h
	quicksort.h
	printList.h
	uniqueListInsert.h
	StringTools.h
	swap.h
	wait.h	
	display.h
	makros.h
	timer.h
	memory.h
	break.h
)

set( SOURCES
	DateTools.cpp
	FileTools.cpp
	MemWatch.cpp
	StringTools.cpp
	binarySearch.cpp
	display.cpp
	timer.cpp
	memory.cpp
	break.cpp
	makros.cpp
)

# Create the library
add_library( Base STATIC ${HEADERS} ${SOURCES} )

set_target_properties(Base PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Base
)

