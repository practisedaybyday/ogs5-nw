/**
 * \file XmlCndInterface.h
 * 2011/11/23 KR as derived class from XMLInterface
 */

#ifndef XMLCNDINTERFACE_H
#define XMLCNDINTERFACE_H

#include "XMLInterface.h"

class FEMCondition;

/**
 * \brief Reads and writes FEM Conditions to and from XML files.
 */
class XmlCndInterface : public XMLInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XmlCndInterface(ProjectData* project, const std::string &schemaFile);

	/// Dummy function so class hierarchy works. This needs to be implemented later.
	int readFile(const QString &fileName)
	{
		std::cout << "There is currently no implementation for XmlCndInterface::readFile(const QString&)." << std::endl;
		return 0;
	}

	/// Reads an xml-file containing containing FEM Conditions such as Boundary- or Initial Conditions
	int readFile(std::vector<FEMCondition*> &conditions, const QString &fileName);

	/// Writes an xml-file containing containing FEM Conditions such as Boundary- or Initial Conditions
	int writeFile(const QString &fileName, const QString &geoName) const;

private:
	/// Read the details of various FEM Conditions from an xml-file
	void readConditions( const QDomNode &condRoot,
	                     std::vector<FEMCondition*> &conditions,
	                     FEMCondition::CondType type);

};

#endif // XMLCNDINTERFACE_H
