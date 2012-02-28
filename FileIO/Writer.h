/**
 * \file Writer.h
 * 13/02/2012 LB Initial implementation
 */

#ifndef WRITER_H
#define WRITER_H

#include <string>
#include <iostream>
#include <sstream>

namespace FileIO
{

/// @brief Base class which enables writing an object to string, stringstream
/// or file. Also formatting (precision, scientific notation of decimal values)
/// can be set.
///
/// When subclassing you only need to implement void write(std::ostream& stream).
class Writer
{
public:
	Writer();
	virtual ~Writer() {};

	/// @brief Writes the object to a string.
	std::string writeToString();

	/// @brief Writes the object to the given file.
	void writeToFile(std::string filename);

	/// @brief Sets the decimal precision.
	void setPrecision(unsigned int precision);

	/// @brief Sets the format (either ios::scientific or ios::fixed);
	void setFormat(std::ios_base::fmtflags flags);

protected:
	/// @brief Writes the object to the given stream.
	/// This method must be implemented by a subclass.
	virtual void write(std::ostream& stream) = 0;
	
	/// @brief The stream to write to.
	std::stringstream _stream;

private:

};

} // namespace FileIO

#endif // WRITER_H