/**
 * \file FileTools.h
 * 26/4/2010 LB Initial implementation
 *
 */

#include "FileTools.h"

#include <iostream>
#include <fstream>
#include <cstdio>

#include <sys/stat.h>

#include "display.h"

/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(std::string const& strFilename)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);

	if(intStat == 0)
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	else
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;

	return blnReturn;
}

bool HasCRInLineEnding(std::string const& strFilename)
{
	std::ifstream is(strFilename.c_str());
	if (!is) {
		ScreenMessage2("*** error: could not open %s\n", strFilename.data());
		return false;
	}

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	bool foundCR = false;
	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\r':
			foundCR = true;
			break;
		case '\n':
		case EOF:
			break;
		}
	}
	is.close();

	return foundCR;
}

