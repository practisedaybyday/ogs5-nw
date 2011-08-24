:: Setup Visual Studio environment
if %1!==! goto msvc2005

if "%1" == "msvc2005" (
  :msvc2005
  call "%VS80COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
  set generator="Visual Studio 8 2005"
)

if "%1" == "msvc2008" (
  call "%VS90COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
  set generator="Visual Studio 9 2008"
)

:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q build_gui

:: Build
mkdir build_gui
cd build_gui
cmake -G %generator% -DOGS_USE_QT=ON -DOGS_PACKAGING=ON ..
cmake ..
devenv OGS-5-GUI.sln /Build Release

:: Package
cmake ..
devenv OGS-5-GUI.sln /Build Release /Project PACKAGE
del CMakeCache.txt
cmake -G %generator% -DOGS_USE_QT=ON -DOGS_PACKAGING=ON -DOGS_PACKAGING_ZIP=ON ..
cmake ..
devenv OGS-5-GUI.sln /Build Release /Project PACKAGE

cd ..
cd scripts
