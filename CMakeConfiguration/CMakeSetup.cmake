
# Set additional CMake modules path
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeConfiguration")

# Adds useful macros and variables
INCLUDE( CMakeConfiguration/Macros.cmake )

# Suppress warning on setting policies
CMAKE_POLICY(SET CMP0011 OLD)

# Suppress warning on add_subdirectory(dir) where dir contains no CMakeLists.txt
IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.7)
	CMAKE_POLICY(SET CMP0014 OLD)
ENDIF ()