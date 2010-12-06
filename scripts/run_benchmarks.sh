#!/usr/bin/env bash

# Clean results directory
rm -rf ../../benchmarks/results/*.html

# Goto sources directory
cd .. >/dev/null

# Get svn information
svn info > svnInfo.txt

# Don´t abort on errors
set +e >/dev/null

# Run FEM benchmarks
cd build_fem
if  [ "$1" = "ex" ]
then
	ctest -R 'EXCEED' -E 'JOD|AllTests|FILE' -j 32 >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' -E 'JOD' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_brns
	ctest -R 'EXCEED' -E 'FILE|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_pqc
	ctest -R 'EXCEED' -E 'FILE|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'EXCEEDING_FILECOMPARE' >> ../benchOut.txt
	
else
	ctest -E 'JOD|AllTests|FILE|EXCEED' -j 32 > ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'JOD|EXCEED' >> ../benchOut.txt
	cd .. >/dev/null
	
	cd build_brns
	ctest -E 'FILE|EXCEED|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
	cd .. >/dev/null

	cd build_pqc
	ctest -E 'FILE|EXCEED|AllTests' -j 32 >> ../benchOut.txt
	ctest -R 'FILECOMPARE' -E 'EXCEED' >> ../benchOut.txt
fi
cd .. >/dev/null

# Print results
cat benchOut.txt

# Send emails on errors
python26 ~/bin/sendmail/sendmail.py

set -e >/dev/null
