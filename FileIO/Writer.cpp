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
	
	if (this->write(_stream))
		return _stream.str();
	else 
		return std::string("");
}

int Writer::writeToFile(std::string filename)
{
	std::string file_content = this->writeToString();
	if (!file_content.empty())
	{
		std::ofstream fileStream;
		fileStream.open (filename.c_str());
	
		// check file stream
		if (!fileStream)
		{
			std::cerr << "Could not open file " << filename << " !" << std::endl;
			return 0;
		}

		fileStream << file_content;

		fileStream.close();
		return 1;
	}
	return 0;
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
