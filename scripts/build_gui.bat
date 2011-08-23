:: Setup Visual Studio environment
call "%VS80COMNTOOLS%\..\..\VC\bin\vcvars32.bat"
call "%VS90COMNTOOLS%\..\..\VC\bin\vcvars32.bat"

:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q build_gui

:: Build
mkdir build_gui
cd build_gui
rem cmake -G "Visual Studio 8 2005" -DOGS_USE_QT=ON -DOGS_PACKAGING=ON ..
cmake -DOGS_USE_QT=ON -DOGS_PACKAGING=ON ..
cmake ..
devenv OGS-5-GUI.sln /Build Release

:: Package
cmake ..
devenv OGS-5-GUI.sln /Build Release /Project PACKAGE
del CMakeCache.txt
rem cmake -G "Visual Studio 8 2005" -DOGS_USE_QT=ON -DOGS_PACKAGING=ON -DOGS_PACKAGING_ZIP=ON ..
cmake -DOGS_USE_QT=ON -DOGS_PACKAGING=ON -DOGS_PACKAGING_ZIP=ON ..
cmake ..
devenv OGS-5-GUI.sln /Build Release /Project PACKAGE
