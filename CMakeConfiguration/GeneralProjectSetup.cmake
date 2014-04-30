# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
IF (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/release/ogs)
ELSE (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/ogs)
ENDIF (MSVC)

# Collect build information such as revision/commit and timestamp
IF (OGS_BUILD_INFO)
	IF(Git_FOUND)
		# Get git commit
		EXECUTE_PROCESS(
			COMMAND ${GIT_EXECUTABLE} "log" "--name-status" "HEAD^..HEAD"
			COMMAND ${GREP_TOOL_PATH} "-m" "1" "commit"
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE GIT_COMMIT_INFO
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
		MESSAGE(STATUS "Git commit: ${GIT_COMMIT_INFO}")
	ENDIF() # GIT_FOUND

	FIND_PATH(HIDDEN_SVN_DIR entries ${CMAKE_SOURCE_DIR}/.svn)
	IF(Subversion_FOUND AND HIDDEN_SVN_DIR)
		Subversion_WC_INFO(${PROJECT_SOURCE_DIR} Project)
		SET(SVN_REVISION ${Project_WC_REVISION})
	ENDIF() # Subversion_FOUND AND HIDDEN_SVN_DIR
	UNSET(HIDDEN_SVN_DIR)

	EXECUTE_PROCESS(
		COMMAND ${DATE_TOOL_PATH} "+%Y-%m-%d" # %H:%M:%S"
		OUTPUT_VARIABLE BUILD_TIMESTAMP
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

ENDIF () # OGS_BUILD_INFO

# This is for Configure.h which is generated later
INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR}/Base )