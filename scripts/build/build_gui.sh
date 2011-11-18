#!/usr/bin/env bash

SOURCE_LOCATION=`pwd`
SOURCE_LOCATION="$SOURCE_LOCATION/../.."

# Parse options
if [ "$OSTYPE" == 'msys' ]; then
	while getopts "a:d:" opt; do
		case $opt in
			a)
				if [ "$OPTARG" == "x32" ]; then
					ARCHITECTURE="x32"
					WIN_ARCHITECTURE="x86"
				elif [ "$OPTARG" == "x64" ]; then
					ARCHITECTURE="x64"
					WIN_ARCHITECTURE="x64"
				else
					echo "$OPTARG is not a valid argument. Specify x32 or x64."
					exit 1
				fi
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
fi

# Cleanup
rm -rf $BUILD_LOCATION
mkdir -p $BUILD_LOCATION && cd $BUILD_LOCATION

# Configure compiler
source $SOURCE_LOCATION/scripts/base/configure_compiler.sh

# CMake
cmake -DOGS_USE_QT=ON -DOGS_PACKAGING=ON -DDOCS_GENERATE_DIAGRAMS=ON -DDOCS_GENERATE_COLLABORATION_GRAPHS=ON -DCMAKE_BUILD_TYPE=Release -G "$CMAKE_GENERATOR" $SOURCE_LOCATION
cmake $SOURCE_LOCATION

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then
	# Installer
	C:/Windows/system32/cmd.exe \/c "devenv OGS.sln /Build Release /Project PACKAGE"
	rm CMakeCache.txt
	# Zip
	cmake -DOGS_USE_QT=ON -DOGS_PACKAGING=ON -DOGS_PACKAGING_ZIP=ON -G "$CMAKE_GENERATOR" $SOURCE_LOCATION
	cmake $SOURCE_LOCATION
	C:/Windows/system32/cmd.exe \/c "devenv OGS.sln /Build Release /Project PACKAGE"
exit
else
	make -j
	cmake $SOURCE_LOCATION
	make -j
	make package
fi

cd "$SOURCE_LOCATION/scripts/build"