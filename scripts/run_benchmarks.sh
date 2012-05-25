#!/usr/bin/env bash

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/.."
source $SOURCE_LOCATION/scripts/base/configure_compiler.sh

# Parse options
while getopts "a:ed:" opt; do
	case $opt in
		e)
			RUN_EXCEEDING=true
			;;
		d)
			BUILD_LOCATION="$SOURCE_LOCATION/$OPTARG"
			;;
		\?)
			echo "Invalid option: -$OPTARG"
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument."
			exit 1
			;;
	esac
done

# Clean results directory
rm -rf $SOURCE_LOCATION/../benchmarks/results/*.html

# Goto sources directory
cd $SOURCE_LOCATION >/dev/null

if [ -d ".svn" ]; then
	# Get svn information
	svn info > $BUILD_LOCATION/svnInfo.txt
elif [ -d ".git" ]; then
	# Get git information
	git log HEAD~1..HEAD > $BUILD_LOCATION/svnInfo.txt
else
	echo "Aborting: Version information not found."
	exit 1
fi

# Don´t abort on errors
set +e >/dev/null

# Run FEM benchmarks
cd $BUILD_LOCATION/build_fem
if  [ $RUN_EXCEEDING ]; then
	ctest -R 'EXCEED' -E 'JOD|Tests|FILE' -j $NUM_PROCESSORS > ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' -E 'JOD' >> ../benchOut.txt
	
	cd $BUILD_LOCATION/build_brns
	ctest -R 'EXCEED' -E 'Tests|FILE' -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
	
	cd $BUILD_LOCATION/build_pqc
	ctest -R 'EXCEED' -E 'Tests|FILE' -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
	
	cd $BUILD_LOCATION/build_gems
	ctest -R 'EXCEED' -E 'Tests|FILE' -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
else
	ctest -E 'Tests|FILE|EXCEED' -j $NUM_PROCESSORS > ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
	
	cd $BUILD_LOCATION/build_brns
	ctest -E 'FILE|EXCEED|Tests' -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt

	cd $BUILD_LOCATION/build_pqc
	ctest -E 'FILE|EXCEED|Tests' -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
	
	cd $BUILD_LOCATION/build_gems
	ctest -E 'FILE|EXCEED|Tests' -j $NUM_PROCESSORS >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
fi
cd $BUILD_LOCATION >/dev/null

# Print results
cat benchOut.txt

cd $SOURCE_LOCATION/scripts
# Send emails on errors
FILESIZE=$(stat -c %s $BUILD_LOCATION/svnInfo.txt)
if [ "$FILESIZE" > "0" ] ; then
  echo "Running process_benchmark_job.rb"
  cd process_benchmark_job
  ruby process_benchmark_job.rb $BUILD_LOCATION/svnInfo.txt $BUILD_LOCATION/benchOut.txt $HUDSON_EMAIL $1
  cd ..
fi

set -e >/dev/null
