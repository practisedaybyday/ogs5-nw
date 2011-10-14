#!/usr/bin/env bash

# Clean results directory
rm -rf ../../benchmarks/results/*.html

# Goto sources directory
cd .. >/dev/null

if [ -d ".svn" ]; then
	# Get svn information
	svn info > svnInfo.txt
elif [ -d ".git" ]; then
	# Get git information
	git log HEAD~1..HEAD > svnInfo.txt
else
	echo "Aborting: Version information not found."
	exit 1
fi

# Don´t abort on errors
set +e >/dev/null

# Run FEM benchmarks
cd build_fem
if  [ "$1" = "ex" ]; then
	ctest -R 'EXCEED' -E 'JOD|AllTests|FILE' -j 32 > ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' -E 'JOD' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_brns
	ctest -R 'EXCEED' -E 'AllTests|FILE' -j 32 >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_pqc
	ctest -R 'EXCEED' -E 'AllTests|FILE' -j 32 >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_gems
	ctest -R 'EXCEED' -E 'AllTests|FILE' -j 32 >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
else
	ctest -E 'AllTests|FILE|EXCEED' -j 32 > ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_brns
	ctest -E 'FILE|EXCEED|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
	cd .. >/dev/null

	cd build_pqc
	ctest -E 'FILE|EXCEED|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_gems
	ctest -E 'FILE|EXCEED|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
fi
cd .. >/dev/null

# Print results
cat benchOut.txt

cd scripts
# Send emails on errors
FILESIZE=$(stat -c %s ./../svnInfo.txt)
if [ "$FILESIZE" > "0" ] ; then
  echo "Running process_benchmark_job.rb"
  cd process_benchmark_job
  ruby process_benchmark_job.rb ./../../svnInfo.txt ./../../benchOut.txt $HUDSON_EMAIL $1
  cd ..
fi

set -e >/dev/null
