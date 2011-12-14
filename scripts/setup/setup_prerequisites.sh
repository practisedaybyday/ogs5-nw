#!/usr/bin/env bash

# Check for prerequisites
WGET_LOCATION=`which wget`
if [ -z "$WGET_LOCATION" ]; then
	echo "wget not found! Aborting..."
	exit 1
fi
CMAKE_LOCATION=`which cmake`
if [ -z "$CMAKE_LOCATION" ]; then
	echo "CMake not found! Aborting..."
	exit 1
fi

source $SOURCE_LOCATION/scripts/base/configure_compiler.sh

## Windows specific
if [ "$OSTYPE" == 'msys' ]; then
	
	# 7-zip
	SEVENZIP_LOCATION=`which 7za`
	if [ ! -z "$SEVENZIP_LOCATION" ]; then
		echo "7-zip found."
	else
		cd ~/bin
		wget --no-check-certificate https://github.com/downloads/ufz/devguide/7za.exe
	fi

	# jom
	JOM_LOCATION=`which jom`
	if [ ! -z "$JOM_LOCATION" ]; then
		echo "jom found."
	else
		cd ~/bin
		wget --no-check-certificate https://github.com/downloads/ufz/devguide/jom.exe
	fi

fi


cd "$SOURCE_LOCATION/scripts/setup"