if [ "$OSTYPE" == 'msys' ]; then
	if [ "$ARCHITECTURE" != "x32" -a "$ARCHITECTURE" != "x64" ]; then
		echo "You did not pass the cpu architecture as a script parameter. Aborting..."
		echo "Usage: $0 -x [x32|x64]"
		exit 1
	fi
fi