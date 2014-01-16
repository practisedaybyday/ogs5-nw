/**
 * \file FileTools.h
 * 26/4/2010 LB Initial implementation
 *
 */

#ifndef FILETOOLS_H
#define FILETOOLS_H

// ** INCLUDES **
#include <string>
#include <sys/stat.h>

/**
 * Returns true if given file exists. From http://www.techbytes.ca/techbyte103.html
 */
bool IsFileExisting(std::string const& strFilename);

#endif // FILETOOLS_H
