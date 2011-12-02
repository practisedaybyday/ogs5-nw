/**
 * \file XmlCndInterface.cpp
 * 2011/11/23 KR as derived class from XMLInterface
 */

#include "XmlCndInterface.h"
#include "FEMCondition.h"

#include <QFile>
#include <QTextCodec>
#include <QtXml/QDomDocument>


XmlCndInterface::XmlCndInterface(ProjectData* project, const std::string &schemaFile)
: XMLInterface(project, schemaFile)
{
}

int XmlCndInterface::readFile(std::vector<FEMCondition*> &conditions, const QString &fileName)
{
	QFile* file = new QFile(fileName);
	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XMLInterface::readFEMCondFile() - Can't open xml-file." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName))
	{
		delete file;
		return 0;
	}

	QDomDocument doc("OGS-Cond-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //root element, used for identifying file-type
	if (docElement.nodeName().compare("OpenGeoSysCond"))
	{
		std::cout << "XMLInterface::readFEMCondFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	//std::vector<FEMCondition*> conditions;
	QDomNodeList lists = docElement.childNodes();
	for (int i = 0; i < lists.count(); i++)
	{
		if (lists.at(i).nodeName().compare("BoundaryConditions") == 0)
			readConditions(lists.at(i), conditions, FEMCondition::BOUNDARY_CONDITION);
		else if (lists.at(i).nodeName().compare("InitialConditions") == 0)
			readConditions(lists.at(i), conditions, FEMCondition::INITIAL_CONDITION);
		else if (lists.at(i).nodeName().compare("SourceTerms") == 0)
			readConditions(lists.at(i), conditions, FEMCondition::SOURCE_TERM);
	}
	if (!conditions.empty())
		return 1;             //do something like _geoObjects->addStationVec(stations, stnName, color);
	else
	{
		std::cout << "XMLInterface::readFEMCondFile() - No FEM Conditions found..." <<
		std::endl;
		return 0;
	}

	delete file;

	return 1;
}

void XmlCndInterface::readConditions( const QDomNode &listRoot,
                                   std::vector<FEMCondition*> &conditions,
                                   FEMCondition::CondType type)
{
	QDomElement cond = listRoot.firstChildElement();
	while (!cond.isNull())
	{
		std::string geometry_name ( cond.attribute("geometry").toStdString() );
		if (this->_project->getGEOObjects()->exists(geometry_name) >= 0)
		{

			FEMCondition* c ( new FEMCondition(geometry_name, type) );

			QDomNodeList condProperties = cond.childNodes();
			for (int i = 0; i < condProperties.count(); i++)
			{
				if (condProperties.at(i).nodeName().compare("Process") == 0)
				{
					QDomNodeList processProps = condProperties.at(i).childNodes();
					for (int j = 0; j < processProps.count(); j++)
					{
						if (processProps.at(j).nodeName().compare("Type") == 0)
							c->setProcessType(FiniteElement::convertProcessType(processProps.at(j).toElement().text().toStdString()));
						else if (processProps.at(j).nodeName().compare("Variable") == 0)
							c->setProcessPrimaryVariable(FiniteElement::convertPrimaryVariable(processProps.at(j).toElement().text().toStdString()));
					}
				}
				else if (condProperties.at(i).nodeName().compare("Geometry") == 0)
				{
					QDomNodeList geoProps = condProperties.at(i).childNodes();
					for (int j = 0; j < geoProps.count(); j++)
					{
						if (geoProps.at(j).nodeName().compare("Type") == 0)
							c->setGeoType(GEOLIB::convertGeoType(geoProps.at(j).toElement().text().toStdString()));
						else if (geoProps.at(j).nodeName().compare("Name") == 0)
							c->setGeoName(geoProps.at(j).toElement().text().toStdString());
					}
				}
				else if (condProperties.at(i).nodeName().compare("Distribution") == 0)
				{
					QDomNodeList distProps = condProperties.at(i).childNodes();
					for (int j = 0; j < distProps.count(); j++)
					{
						if (distProps.at(j).nodeName().compare("Type") == 0)
							c->setProcessDistributionType(FiniteElement::convertDisType(distProps.at(j).toElement().text().toStdString()));
						else if (distProps.at(j).nodeName().compare("Value") == 0)
							c->setDisValue(strtod(distProps.at(j).toElement().text().toStdString().c_str(), 0));
					}
				}
			}
			conditions.push_back(c);
		}
		else
		{
			std::cout << "Error loading FEM Conditions: No geometry \"" << geometry_name << "\" found." << std::endl;
		}
		cond = cond.nextSiblingElement();
	}
}

int XmlCndInterface::writeFile(const QString &fileName, const QString &geoName, FEMCondition::CondType type) const
{
	std::fstream stream(fileName.toStdString().c_str(), std::ios::out);
	if (!stream.is_open())
	{
		std::cout << "XmlCndInterface::writeFile() - Could not open file...\n";
		return 0;
	}
	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysCND.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-CND-DOM");
	QDomElement root = doc.createElement("OpenGeoSysCond");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.net" );
	root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
	root.setAttribute( "xsi:noNamespaceSchemaLocation", "http://141.65.34.25/OpenGeoSysCND.xsd" );

	const std::vector<FEMCondition*> conditions (_project->getConditions(FiniteElement::INVALID_PROCESS, geoName.toStdString(), type) );

	if (conditions.empty()) return 1;

	doc.appendChild(root);


	size_t nConditions (conditions.size());
	for (size_t i=0; i<nConditions; i++)
	{
		FEMCondition::CondType current_type = conditions[i]->getCondType();
		if (current_type == type || type == FEMCondition::UNSPECIFIED)
		{
			QDomElement listTag;
			QString condText;

			if (current_type == FEMCondition::BOUNDARY_CONDITION)
			{
				listTag = this->getCondListElement(doc, root, "BoundaryConditions");
				condText = "BC";
			}
			else if (current_type == FEMCondition::INITIAL_CONDITION)
			{
				listTag = this->getCondListElement(doc, root, "InitialConditions");
				condText = "IC";
			}
			else if (current_type == FEMCondition::SOURCE_TERM)
			{
				listTag = this->getCondListElement(doc, root, "SourceTerms");
				condText = "ST";
			}
			else 
			{
				std::cout << "Error in XmlCndInterface::writeFile() - Unspecified FEMConditions found ... Abort writing." << std::endl;
				return 0;
			}
			this->writeCondition(doc, listTag, conditions[i], condText, geoName);
		}
	}
	std::string xml = doc.toString().toStdString();
	stream << xml;
	stream.close();

	return 1;
}

void XmlCndInterface::writeCondition( QDomDocument doc, QDomElement &listTag, const FEMCondition* cond, const QString &condText, const QString &geometryName) const
{
	QString geoName (QString::fromStdString(cond->getAssociatedGeometryName()));

	if (geoName.compare(geometryName) != 0)
	{
		std::cout << "Geometry name not matching, skipping condition \"" << cond->getGeoName() << "\"..." << std::endl;
	}

	QDomElement condTag ( doc.createElement(condText) );
	condTag.setAttribute("geometry", geoName);
	listTag.appendChild(condTag);

	QDomElement processTag ( doc.createElement("Process") );
	condTag.appendChild(processTag);
	QDomElement processTypeTag ( doc.createElement("Type") );
	processTag.appendChild(processTypeTag);
	QDomText processTypeText ( doc.createTextNode(
		QString::fromStdString(FiniteElement::convertProcessTypeToString(cond->getProcessType()))) );
	processTypeTag.appendChild(processTypeText);
	QDomElement processPVTag ( doc.createElement("Variable") );
	processTag.appendChild(processPVTag);
	QDomText processPVText ( doc.createTextNode(
		QString::fromStdString(FiniteElement::convertPrimaryVariableToString(cond->getProcessPrimaryVariable()))) );
	processPVTag.appendChild(processPVText);

	QDomElement geoTag ( doc.createElement("Geometry") );
	condTag.appendChild(geoTag);
	QDomElement geoTypeTag ( doc.createElement("Type") );
	geoTag.appendChild(geoTypeTag);
	QDomText geoTypeText ( doc.createTextNode(
		QString::fromStdString(GEOLIB::convertGeoTypeToString(cond->getGeoType()))) );
	geoTypeTag.appendChild(geoTypeText);
	QDomElement geoNameTag ( doc.createElement("Name") );
	geoTag.appendChild(geoNameTag);
	QString geo_obj_name ( QString::fromStdString(cond->getGeoName()) );
	QDomText geoNameText ( doc.createTextNode(geo_obj_name) );
	geoNameTag.appendChild(geoNameText);

	QDomElement disTag ( doc.createElement("Distribution") );
	condTag.appendChild(disTag);
	QDomElement disTypeTag ( doc.createElement("Type") );
	disTag.appendChild(disTypeTag);
	QDomText disTypeText ( doc.createTextNode(
		QString::fromStdString(FiniteElement::convertDisTypeToString(cond->getProcessDistributionType()))) );
	disTypeTag.appendChild(disTypeText);
	QDomElement disValueTag ( doc.createElement("Value") );
	disTag.appendChild(disValueTag);
	double dis_value (cond->getDisValue()[0]); //TODO: do this correctly!
	QDomText disValueText ( doc.createTextNode(QString::number(dis_value)) );
	disValueTag.appendChild(disValueText);
}

QDomElement XmlCndInterface::getCondListElement( QDomDocument doc, QDomElement &root, const QString &text ) const
{
	QDomNodeList list = root.elementsByTagName(text);
	if (list.isEmpty()) 
	{ 
		QDomElement newListTag ( doc.createElement(text) );
		root.appendChild(newListTag);
		return newListTag;
	}
	return list.at(0).toElement();
}