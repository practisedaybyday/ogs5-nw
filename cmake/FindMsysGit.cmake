# The module defines the following variables:
#   MSYSGIT_BIN_DIR - path to the tool binaries
#   MSYSGIT_FOUND - true if the command line client was found
# Example usage:
#   FIND_PACKAGE(MsysGit)
#   IF(MSYSGIT_FOUND)
#     MESSAGE("msysGit tools found in: ${MSYSGIT_BIN_DIR}")
#   ENDIF()


FIND_PATH(MSYSGIT_BIN_DIR
	NAMES bash.exe git.exe PATH_SUFFIXES Git/bin)

# Handle the QUIETLY and REQUIRED arguments and set MSYSGIT_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MsysGit
                                  REQUIRED_VARS MSYSGIT_BIN_DIR)