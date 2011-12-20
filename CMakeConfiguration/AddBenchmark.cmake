# This script is called from AddBenchmark in Macros.cmake
# It deletes the benchmark output files and then runs the benchmark.
IF (WIN32)

	SEPARATE_ARGUMENTS(FILES_TO_DELETE_VARS WINDOWS_COMMAND ${FILES_TO_DELETE})

	EXECUTE_PROCESS (
		COMMAND del /S /Q ${FILES_TO_DELETE_VARS}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks)

	EXECUTE_PROCESS (
		COMMAND ${EXECUTABLE_OUTPUT_PATH}/Release/ogs ${benchmarkStrippedName}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
		TIMEOUT ${BENCHMARK_TIMEOUT})

ELSE (WIN32)

	SEPARATE_ARGUMENTS(FILES_TO_DELETE_VARS UNIX_COMMAND ${FILES_TO_DELETE})

	EXECUTE_PROCESS (
		COMMAND rm -f ${FILES_TO_DELETE_VARS}
		WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks)
	IF(OGS_PROFILE)
		MESSAGE(STATUS "Profiling benchmark")
		IF(OGS_OUTPUT_PROFILE)
			MESSAGE(STATUS "Executing gprof2dot.py")
			EXECUTE_PROCESS (
				COMMAND ${EXECUTABLE_OUTPUT_PATH}/ogs ${benchmarkStrippedName}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
				TIMEOUT ${BENCHMARK_TIMEOUT})

			# Run gprof2dot.py
			EXECUTE_PROCESS (
				COMMAND ${GPROF_PATH} ${EXECUTABLE_OUTPUT_PATH}/ogs
				COMMAND ${PROJECT_SOURCE_DIR}/scripts/gprof2dot.py -s -n 5.0 -e 1.0
				COMMAND ${DOT_TOOL_PATH} -Tpng -o ${PROJECT_SOURCE_DIR}/../benchmarks/results/${benchmarkStrippedName}.png
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir})

		ELSE()
			EXECUTE_PROCESS (
				COMMAND ${EXECUTABLE_OUTPUT_PATH}/ogs ${benchmarkStrippedName}
				WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
				TIMEOUT ${BENCHMARK_TIMEOUT})
		ENDIF()
	ELSE()
		EXECUTE_PROCESS (
			COMMAND ${EXECUTABLE_OUTPUT_PATH}/ogs ${benchmarkStrippedName}
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/../benchmarks/${benchmarkDir}
			TIMEOUT ${BENCHMARK_TIMEOUT})
	ENDIF()

ENDIF (WIN32)
