if [ "$OSTYPE" == 'msys' ]; then
	source $SOURCE_LOCATION/scripts/base/configure_win_vs.sh
else
	CMAKE_GENERATOR="Unix Makefiles"
fi