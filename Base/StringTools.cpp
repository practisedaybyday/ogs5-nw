/*
 * StringTools.cpp
 *
 *  Created on: Jun 16, 2010
 *      Author: TF
 */

#include "StringTools.h"

std::list<std::string> splitString(const std::string &str, char delim)
{
	std::list<std::string> strList;
	std::stringstream ss(str);
	std::string item;
	while(getline(ss, item, delim))
		strList.push_back(item);
	return strList;
}

std::string replaceString(const std::string &searchString,
                          const std::string &replaceString,
                          std::string stringToReplace)
{
	std::string::size_type pos = stringToReplace.find(searchString, 0);
	size_t intLengthSearch = searchString.length();

	while (std::string::npos != pos)
	{
		stringToReplace.replace(pos, intLengthSearch, replaceString);
		pos = stringToReplace.find(searchString, 0);
	}
	return stringToReplace;
}

void trim(std::string &str, char ch)
{
	std::string::size_type pos = str.find_last_not_of(ch);
	if(pos != std::string::npos)
	{
		str.erase(pos + 1);
		pos = str.find_first_not_of(ch);
		if(pos != std::string::npos)
			str.erase(0, pos);
	}
	else
		str.erase(str.begin(), str.end());
}

std::string getFileNameFromPath(const std::string &str)
{
	std::string::size_type beg1 = str.find_last_of('/');
	std::string::size_type beg2 = str.find_last_of('\\');
	std::string::size_type beg;
	if (beg1 == std::string::npos && beg2 == std::string::npos) beg = 0;
	else if (beg1 == std::string::npos) beg = beg2;
	else if (beg2 == std::string::npos) beg = beg1;
	else beg = (beg1<beg2) ? beg2 : beg1;
	std::string file ( str.substr(beg+1) );
	std::string::size_type end  = file.find_last_of('.');
	return file.substr(0,end);
}

#ifdef MSVC
void correctScientificNotation(std::string filename, size_t precision)
{
	std::ifstream stream;
	std::ofstream outputStream;

	stream.open(filename.c_str());
	std::string tmpFilename = filename + ".tmp";
	outputStream.open(tmpFilename.c_str());

	if (!stream)
	{
		std::cout << "correctScientificNotation: fstream is not open" << std::endl;
		return;
	}

	std::string line;

	// Iterate over lines in stream
	while (getline(stream, line))
	{
		std::string word;
		std::istringstream iss(line);
		// Iterate over all words in line
		while (iss >> word)
		{
			// Search for e+0
			std::size_t exponentPosition = word.find("e+0", precision);
			if (exponentPosition == std::string::npos)
				// If not found search for e-0
				exponentPosition = word.find("e-0", precision);
			if (exponentPosition != std::string::npos)
			{
				std::size_t wordSize = word.size();
				std::size_t exponentSize = wordSize - exponentPosition;

				if(exponentSize > 4)
				{
					// Erase the leading zero considering trailing characters
					size_t i = wordSize - 1;
					while (!isdigit(word[i]))
						--i;

					size_t erasePos = wordSize - 3 - (wordSize - 1 - i);
					std::string eraseString = word.substr(erasePos, 1);
					if (eraseString.find("0") != std::string::npos)
						word.erase(erasePos, 1);
				}
			}

			outputStream << word << " ";
		}
		outputStream << std::endl;
	}

	stream.close();
	outputStream.close();

	remove(filename.c_str());
	rename(tmpFilename.c_str(), filename.c_str());
}
#endif

