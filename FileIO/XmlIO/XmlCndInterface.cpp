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
			//bool cancelLoading = OGSError();
			std::cout << "Error loading FEM Conditions: No geometry \"" << geometry_name << "\" found." << std::endl;
		}
		cond = cond.nextSiblingElement();
	}
}

int XmlCndInterface::writeFile(const QString &fileName, const QString &geoName) const
{
	Q_UNUSED(fileName)
	Q_UNUSED(geoName)
	return 0;
}
