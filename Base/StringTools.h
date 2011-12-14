#ifndef STRINGTOOLS_H
#define STRINGTOOLS_H

#include <ctype.h>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>

/**
 *   Splits a string into a list of strings.
 *  \param str String to be splitted
 *  \param delim Character indicating that the string should be splitted
 *  \return
 */
std::list<std::string> splitString(const std::string &str, char delim);

/**
 *   Replaces a substring with another in a string
 *  \param searchString Search for this string
 *  \param replaceString Replace with this string
 *  \param stringToReplace Search and replace in this string
 *  \return The modified string
 */
std::string replaceString(const std::string &searchString,
                          const std::string &replaceString,
                          std::string stringToReplace);

/**
 *   Converts a number (double, float, int, ...) into a string
 *  \param d The number to be converted
 *  \return The number as string
 */
template<typename T> std::string number2str(T d)
{
	std::stringstream out;
	out << d;
	return out.str();
}

/**
 *   Converts a string into a number (double, float, int, ...)
 *  Example: size_t number (str2number<size_t> (str));
 *  \param str string to be converted
 *  \return the number
 */
template<typename T> T str2number (const std::string &str)
{
	std::stringstream strs (str, std::stringstream::in | std::stringstream::out);
	T v;
	strs >> v;
	return v;
}

/**
 * Strip whitespace (or other characters) from the beginning and end of a string.
 */
void trim(std::string &str, char ch = ' ');

/**
 * Extract the filename from a path
 */
std::string getFileNameFromPath(const std::string &str, bool with_extension = false);



#ifdef MSVC
void correctScientificNotation(std::string filename, size_t precision = 0);
#endif

#endif //STRINGTOOLS_H
