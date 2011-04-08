:: Setup Visual Studio environment
call "%VS80COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
call "%VS90COMNTOOLS%\..\..\VC\bin\vcvars32.bat"

:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q build_fem

:: Build
mkdir build_fem
cd build_fem
cmake -DOGS_FEM=ON -DOGS_PACKAGING=ON -DOGS_PACKAGING_ZIP=ON ..
cmake ..
devenv OGS-FEM-5.sln /Build Release

:: Package
cmake ..
devenv OGS-FEM-5.sln /Build Release /Project PACKAGE
