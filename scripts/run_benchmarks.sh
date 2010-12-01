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
ctest -R 'NB|MCGB|FS|AKS|YW|UJG|MDL|NW|WW|HS|PCH|CB' -E 'FILE|EXCEED|thm_decov|h_us_line|h_us_quad|h_us_tri|w_exp' -j 32 > ../benchOut.txt
ctest -R FILECOMPARE -E 'JOD|EXCEED|thm_decov|h_us_line|h_us_quad|h_us_tri|w_exp' >> ../benchOut.txt
cd .. >/dev/null

# Run FEM_BRNS benchmarks
cd build_brns
ctest -R FC -E "FILE|EXCEED" -j 32 >> ../benchOut.txt
ctest -R FILECOMPARE -E "EXCEED" >> ../benchOut.txt
cd .. >/dev/null

# Run FEM_PQC benchmarks
#cd build_pqc
#ctest -E "FILE|EXCEED" -j 32 >> ../benchOut.txt
#ctest -R FILECOMPARE -E "EXCEED" >> ../benchOut.txt
#cd .. >/dev/null


# Print results
cat benchOut.txt

# Send emails on errors
python26 ~/bin/sendmail/sendmail.py

set -e >/dev/null
