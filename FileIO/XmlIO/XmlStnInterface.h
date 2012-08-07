/**
 * \file XmlStnInterface.h
 * 2011/11/23 KR as derived class from XMLInterface
 */

#ifndef XMLSTNINTERFACE_H
#define XMLSTNINTERFACE_H

#include "XMLInterface.h"

#include "RapidXML/rapidxml.hpp"

namespace FileIO
{

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

	/// Reads an xml-file containing station object definitions into the GEOObjects used in the contructor (requires Qt)
	int readFile(const QString &fileName);

	/// Reads an xml-file using the RapidXML parser integrated in the source code (i.e. this function is usable without Qt)
	int rapidReadFile(const std::string &fileName);

protected:
	int write(std::ostream& stream);

private:
	/// Reads GEOLIB::Station- or StationBorehole-objects from an xml-file
	void readStations  ( const QDomNode &stationsRoot, std::vector<GEOLIB::Point*>* stations, const std::string &filename);

	/// Writes borehole-specific data to a station-xml-file.
	void writeBoreholeData(QDomDocument &doc,
	                       QDomElement &boreholeTag,
	                       GEOLIB::StationBorehole* borehole) const;

	/// Reads the stratigraphy of a borehole from an xml-file
	void readStratigraphy( const QDomNode &stratRoot, GEOLIB::StationBorehole*  borehole );

	/// Reads GEOLIB::Station- or StationBorehole-objects from an xml-file using the RapidXML parser
	void rapidReadStations(const rapidxml::xml_node<>* station_root, std::vector<GEOLIB::Point*> *stations, const std::string &file_name);
	
	/// Reads the stratigraphy of a borehole from an xml-file using the RapidXML parser
	void rapidReadStratigraphy(const rapidxml::xml_node<>* strat_root, GEOLIB::StationBorehole* borehole);
};

}

#endif // XMLSTNINTERFACE_H
