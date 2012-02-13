/**
 * \file Writer.cpp
 * 13/02/2012 LB Initial implementation
 * 
 * Implementation of the Writer class
 */

// ** INCLUDES **
#include "Writer.h"

#include <fstream>

namespace FileIO
{

Writer::Writer()
{
}

std::string Writer::writeToString()
{
	// Empty stream and clear error states.
	_stream.str("");
	_stream.clear();
	
	this->write(_stream);
	return _stream.str();
}

void Writer::writeToFile(std::string filename)
{
	std::ofstream fileStream;
	fileStream.open (filename.c_str());
	
	// check file stream
	if (!fileStream)
	{
		std::cerr << "Could not open file " << filename << " !" << std::endl;
		return;
	}

	fileStream << this->writeToString();

	fileStream.close();
}

void Writer::setPrecision(unsigned int precision)
{
	_stream.precision(precision);
}

void Writer::setFormat(std::ios_base::fmtflags flags)
{
	_stream.setf(flags);
}

} // namespace FileIO
