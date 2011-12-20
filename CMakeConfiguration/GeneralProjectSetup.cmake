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
	    execute_process(
			COMMAND "date" "+%Y-%m-%d %H:%M:%S"
			OUTPUT_VARIABLE build_timestamp
			OUTPUT_STRIP_TRAILING_WHITESPACE
		)
	    if( build_timestamp )
	       set( BUILD_TIMESTAMP "${build_timestamp}" )
	       # message( STATUS "Build timestamp: ${BUILD_TIMESTAMP}" )
	    endif()
	  find_package( Git )
	  if(GIT_FOUND)
	    # Get git commit
	    execute_process(
			COMMAND "git" "log" "--name-status" "HEAD^..HEAD"
			COMMAND "grep" "-m" "1" "commit"
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			OUTPUT_VARIABLE git_commit_info
			OUTPUT_STRIP_TRAILING_WHITESPACE
			ERROR_VARIABLE not_git_repo
		)
	    if( git_commit_info )
	      # message( STATUS "Git commit: " ${git_commit_info} )
	      set( GIT_COMMIT_INFO "${git_commit_info}" )
	    else( git_commit_info )
	      if( not_git_repo )
		message( STATUS "Git not versioning the source" )
	      endif( not_git_repo )
	    endif( git_commit_info )

	    # Get git branch
			execute_process(
				COMMAND "git" "branch"
				COMMAND "grep" "\\*"
				WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
				OUTPUT_STRIP_TRAILING_WHITESPACE
				OUTPUT_VARIABLE git_branch_info)
			if (git_branch_info)
				set (GIT_BRANCH_INFO "${git_branch_info}")
			else(git_branch_info)
				message(STATUS "Git branch could not be determined!")
			endif(git_branch_info)
	  else( GIT_FOUND )
		MESSAGE(STATUS "Git not found!")
	  endif( GIT_FOUND )
	endif( NOT MSVC )
ENDIF (OGS_BUILD_INFO)

# This is for Configure.h which is generated later
INCLUDE_DIRECTORIES( ${PROJECT_BINARY_DIR}/Base )