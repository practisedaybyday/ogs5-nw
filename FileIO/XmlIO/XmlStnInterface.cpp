/**
 * \file XmlStnInterface.cpp
 * 2011/11/23 KR as derived class from XMLInterface
 */

#include "XmlStnInterface.h"
#include "DateTools.h"

#include <QFile>
#include <QtXml/QDomDocument>


XmlStnInterface::XmlStnInterface(ProjectData* project, const std::string &schemaFile)
: XMLInterface(project, schemaFile)
{
}

int XmlStnInterface::readFile(const QString &fileName)
{
	GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();
	QFile* file = new QFile(fileName);
	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XmlStnInterface::readFile() - Can't open xml-file." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName))
	{
		delete file;
		return 0;
	}

	QDomDocument doc("OGS-STN-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //root element, used for identifying file-type
	if (docElement.nodeName().compare("OpenGeoSysSTN"))
	{
		std::cout << "XmlStnInterface::readFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	QDomNodeList lists = docElement.childNodes();
	for (int i = 0; i < lists.count(); i++)
	{
		// read all the station lists
		QDomNodeList stationList = lists.at(i).childNodes();
		std::vector<GEOLIB::Point*>* stations = new std::vector<GEOLIB::Point*>;
		std::string stnName("[NN]");

		for (int j = 0; j < stationList.count(); j++)
		{
			if (stationList.at(j).nodeName().compare("name") == 0)
				stnName = stationList.at(j).toElement().text().toStdString();
			else if (stationList.at(j).nodeName().compare("stations") == 0)
				readStations(stationList.at(j), stations);
			else if (stationList.at(j).nodeName().compare("boreholes") == 0)
				readStations(stationList.at(j), stations);
		}

		GEOLIB::Color* color = GEOLIB::getRandomColor();
		if (!stations->empty())
			geoObjects->addStationVec(stations, stnName, color);
		else
			delete stations;
	}

	delete file;

	return 1;
}

void XmlStnInterface::readStations( const QDomNode &stationsRoot,
                                 std::vector<GEOLIB::Point*>* stations )
{
	QDomElement station = stationsRoot.firstChildElement();
	while (!station.isNull())
	{
		if (station.hasAttribute("id") && station.hasAttribute("x") &&
		    station.hasAttribute("y"))
		{
			std::string stationName("[NN]");
			std::string boreholeDate("0000-00-00");
			double boreholeDepth(0.0);

			QDomNodeList stationFeatures = station.childNodes();
			for(int i = 0; i < stationFeatures.count(); i++)
			{
				if (stationFeatures.at(i).nodeName().compare("name") == 0)
					stationName =
					        stationFeatures.at(i).toElement().text().
					        toStdString();
				/* add other station features here */

				else if (stationFeatures.at(i).nodeName().compare("bdepth") == 0)
					boreholeDepth = strtod(stationFeatures.at(
					                               i).toElement().text().
					                       toStdString().c_str(), 0);
				else if (stationFeatures.at(i).nodeName().compare("bdate") == 0)
					boreholeDate  =
					        stationFeatures.at(i).toElement().text().
					        toStdString();
				/* add other borehole features here */
			}

			double zVal = (station.hasAttribute("z")) ? strtod((station.attribute(
			                                                            "z")).
			                                                   toStdString().c_str(),
			                                                   0) : 0.0;

			if (station.nodeName().compare("station") == 0)
			{
				GEOLIB::Station* s =
				        new GEOLIB::Station(strtod((station.attribute("x")).
				                                   toStdString().
				                                   c_str(), 0),
				                            strtod((station.attribute(
				                                            "y")).toStdString().
				                                   c_str(), 0),
				                            zVal, stationName);
				stations->push_back(s);
			}
			else if (station.nodeName().compare("borehole") == 0)
			{
				GEOLIB::StationBorehole* s = GEOLIB::StationBorehole::createStation(
				        stationName,
				        strtod((
				                       station.attribute(
				                               "x")).toStdString().c_str(), 0),
				        strtod((
				                       station.attribute(
				                               "y")).toStdString().c_str(), 0),
				        zVal,
				        boreholeDepth,
				        boreholeDate);
				/* add stratigraphy to the borehole */
				for(int i = 0; i < stationFeatures.count(); i++)
					if (stationFeatures.at(i).nodeName().compare("strat") == 0)
						this->readStratigraphy(stationFeatures.at(i), s);

				stations->push_back(s);
			}
		}
		else
			std::cout <<
			"XmlStnInterface::readStations() - Attribute missing in <station> tag ..." <<
			std::endl;
		station = station.nextSiblingElement();
	}
}

void XmlStnInterface::readStratigraphy( const QDomNode &stratRoot, GEOLIB::StationBorehole* borehole )
{
	//borehole->addSoilLayer((*borehole)[0], (*borehole)[1], (*borehole)[2], "");
	QDomElement horizon = stratRoot.firstChildElement();
	while (!horizon.isNull())
	{
		if (horizon.hasAttribute("id") && horizon.hasAttribute("x") &&
		    horizon.hasAttribute("y") && horizon.hasAttribute("z"))
		{
			std::string horizonName("[NN]");

			QDomNodeList horizonFeatures = horizon.childNodes();
			for(int i = 0; i < horizonFeatures.count(); i++)
				if (horizonFeatures.at(i).nodeName().compare("name") == 0)
					horizonName =
					        horizonFeatures.at(i).toElement().text().
					        toStdString();
				/* add other horizon features here */
			borehole->addSoilLayer(strtod((horizon.attribute("x")).toStdString().c_str(),
			                              0),
			                       strtod((horizon.attribute("y")).toStdString().c_str(),
			                              0),
			                       strtod((horizon.attribute("z")).toStdString().c_str(),
			                              0),
			                       horizonName);
		}
		else
			std::cout <<
			"XmlStnInterface::readStratigraphy() - Attribute missing in <horizon> tag ..."
			          << std::endl;
		horizon = horizon.nextSiblingElement();
	}
}

int XmlStnInterface::writeFile(const QString &filename, const QString &stnName) const
{
	GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();
	std::fstream stream(filename.toStdString().c_str(), std::ios::out);
	if (!stream.is_open())
	{
		std::cout << "XmlStnInterface::writeFile() - Could not open file...\n";
		return 0;
	}
	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysSTN.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-STN-DOM");
	QDomElement root = doc.createElement("OpenGeoSysSTN");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.net" );
	root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
	root.setAttribute( "xsi:noNamespaceSchemaLocation", "http://141.65.34.25/OpenGeoSysSTN.xsd" );

	const std::vector<GEOLIB::Point*>* stations (geoObjects->getStationVec(stnName.toStdString()));
	bool isBorehole =
	        (static_cast<GEOLIB::Station*>((*stations)[0])->type() ==
	         GEOLIB::Station::BOREHOLE) ? true : false;

	doc.appendChild(root);
	QDomElement stationListTag = doc.createElement("stationlist");
	root.appendChild(stationListTag);

	QDomElement listNameTag = doc.createElement("name");
	stationListTag.appendChild(listNameTag);
	QDomText stationListNameText = doc.createTextNode(stnName);
	listNameTag.appendChild(stationListNameText);
	QString listType = (isBorehole) ? "boreholes" : "stations";
	QDomElement stationsTag = doc.createElement(listType);
	stationListTag.appendChild(stationsTag);

	size_t nStations(stations->size());
	for (size_t i = 0; i < nStations; i++)
	{
		QString stationType =  (isBorehole) ? "borehole" : "station";
		QDomElement stationTag = doc.createElement(stationType);
		stationTag.setAttribute( "id", QString::number(i) );
		stationTag.setAttribute( "x",  QString::number((*(*stations)[i])[0], 'f') );
		stationTag.setAttribute( "y",  QString::number((*(*stations)[i])[1], 'f') );
		stationTag.setAttribute( "z",  QString::number((*(*stations)[i])[2], 'f') );
		stationsTag.appendChild(stationTag);

		QDomElement stationNameTag = doc.createElement("name");
		stationTag.appendChild(stationNameTag);
		QDomText stationNameText =
		        doc.createTextNode(QString::fromStdString(static_cast<GEOLIB::Station*>((*
		                                                                                 stations)
		                                                                                [
		                                                                                        i
		                                                                                ])
		                                                  ->getName()));
		stationNameTag.appendChild(stationNameText);

		if (isBorehole)
			writeBoreholeData(doc, stationTag,
			                  static_cast<GEOLIB::StationBorehole*>((*stations)[i]));
	}

	std::string xml = doc.toString().toStdString();
	stream << xml;
	stream.close();
	return 1;
}

void XmlStnInterface::writeBoreholeData(QDomDocument &doc,
                                     QDomElement &boreholeTag,
                                     GEOLIB::StationBorehole* borehole) const
{
	QDomElement stationDepthTag = doc.createElement("bdepth");
	boreholeTag.appendChild(stationDepthTag);
	QDomText stationDepthText = doc.createTextNode(QString::number(borehole->getDepth(), 'f'));
	stationDepthTag.appendChild(stationDepthText);
	if (borehole->getDate() != 0)
	{
		QDomElement stationDateTag = doc.createElement("bdate");
		boreholeTag.appendChild(stationDateTag);
		QDomText stationDateText =
		        doc.createTextNode(QString::fromStdString(date2string(borehole->getDate())));
		stationDateTag.appendChild(stationDateText);
	}

	std::vector<GEOLIB::Point*> profile = borehole->getProfile();
	std::vector<std::string> soilNames = borehole->getSoilNames();
	size_t nHorizons(profile.size());

	if (nHorizons > 1)
	{
		QDomElement stratTag = doc.createElement("strat");
		boreholeTag.appendChild(stratTag);

		for (size_t j = 1; j < nHorizons; j++) /// the first entry in the profile vector is just the position of the borehole
		{
			QDomElement horizonTag = doc.createElement("horizon");
			horizonTag.setAttribute( "id", QString::number(j) );
			horizonTag.setAttribute( "x",  QString::number((*profile[j])[0], 'f') );
			horizonTag.setAttribute( "y",  QString::number((*profile[j])[1], 'f') );
			horizonTag.setAttribute( "z",  QString::number((*profile[j])[2], 'f') );
			stratTag.appendChild(horizonTag);
			QDomElement horizonNameTag = doc.createElement("name");
			horizonTag.appendChild(horizonNameTag);
			QDomText horizonNameText =
			        doc.createTextNode(QString::fromStdString(soilNames[j]));
			horizonNameTag.appendChild(horizonNameText);
		}
	}
}
