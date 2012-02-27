/**
 * \file XmlStnInterface.h
 * 2011/11/23 KR as derived class from XMLInterface
 */

#ifndef XMLSTNINTERFACE_H
#define XMLSTNINTERFACE_H

#include "XMLInterface.h"

/**
 * \brief Reads and writes Observation Sites to and from XML files.
 */
class XmlStnInterface : public XMLInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XmlStnInterface(ProjectData* project, const std::string &schemaFile);

	/// Reads an xml-file containing station object definitions into the GEOObjects used in the contructor
	int readFile(const QString &fileName);

protected:
	int write(std::ostream& stream);

private:
	/// Reads GEOLIB::Station- or StationBorehole-objects from an xml-file
	void readStations  ( const QDomNode &stationsRoot, std::vector<GEOLIB::Point*>* stations );

	/// Writes borehole-specific data to a station-xml-file.
	void writeBoreholeData(QDomDocument &doc,
	                       QDomElement &boreholeTag,
	                       GEOLIB::StationBorehole* borehole) const;

	/// Reads the stratigraphy of a borehole from an xml-file
	void readStratigraphy( const QDomNode &stratRoot, GEOLIB::StationBorehole* borehole );
};

#endif // XMLSTNINTERFACE_H
