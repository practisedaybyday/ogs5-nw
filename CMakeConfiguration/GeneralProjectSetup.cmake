# Set build directories
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
IF (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/release/ogs)
ELSE (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/ogs)
ENDIF (MSVC)

IF (OGS_BUILD_INFO)
	find_package( Subversion )
	if(Subversion_FOUND)
	  #Subversion_WC_INFO(${CMAKE_SOURCE_DIR} "svn") , not a working copy
	  execute_process(
	    COMMAND "svnversion" "--no-newline" ${CMAKE_SOURCE_DIR}
	    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	    OUTPUT_VARIABLE svn_revision
	    ERROR_VARIABLE svn_error
	    )
	  set( SVN_REVISION "${svn_revision}" )
	endif( Subversion_FOUND )

	if( NOT MSVC )
	  # see if there is a file that will work with find_package
	  find_file( GIT FindGit.cmake CMakeConfiguration/ )
	  if( GIT STREQUAL GIT-NOTFOUND )
	    message( WARNING "there is no FindGit.cmake file, so I am downloading one that is licensed under the GPL, copyright Szilard Pall (pszilard@cbr.su.se).  Do not commit it to the repository." )
	    file( DOWNLOAD http://repo.or.cz/w/gromacs.git/blob_plain/f81efa276ee7e1c41197d9d8b56b91c7764e0641:/cmake/FindGit.cmake ${PROJECT_SOURCE_DIR}/CMakeConfiguration/FindGit.cmake TIMEOUT 10 )
	  endif( GIT STREQUAL GIT-NOTFOUND )
	    execute_process(
	      COMMAND "date" "+%Y-%m-%d %H:%M:%S"
	      OUTPUT_VARIABLE build_timestamp
	      OUTPUT_STRIP_TRAILING_WHITESPACE
	      )
	    if( build_timestamp )
	       set( BUILD_TIMESTAMP "${build_timestamp}" )
	       message( STATUS "Build timestamp: ${BUILD_TIMESTAMP}" )
	    endif()
	  find_package( Git )
	  if(Git_FOUND)
	    execute_process(
	      COMMAND "git" "describe" "--always"
	      COMMAND "xargs" "-Ixx" "git" "name-rev" "xx"
	      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	      OUTPUT_VARIABLE git_commit_info
	      OUTPUT_STRIP_TRAILING_WHITESPACE
	      ERROR_VARIABLE not_git_repo
	      )
	    if( git_commit_info )
	      #message( STATUS "Git executable version: " ${Git_VERSION} )
	      message( STATUS "Git commit: " ${git_commit_info} )
	      set( GIT_COMMIT_INFO "${git_commit_info}" )
	    else( git_commit_info )
	      if( not_git_repo )
		message( STATUS "Git not versioning the source" )
	      endif( not_git_repo )
	    endif( git_commit_info )
	  endif( Git_FOUND )
	endif( NOT MSVC )
ENDIF (OGS_BUILD_INFO)

# This is for Configure.h which is generated later
INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR}/Base )