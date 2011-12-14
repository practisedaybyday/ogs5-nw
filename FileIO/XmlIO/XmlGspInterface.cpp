/**
 * \file XmlGspInterface.cpp
 * 2011/11/23 KR as derived class from XMLInterface
 */

#include "XmlGspInterface.h"

#include "XmlGmlInterface.h"
#include "XmlStnInterface.h"
#include "XmlCndInterface.h"

#include <QFileInfo>
#include <QFile>
#include <QtXml/QDomDocument>

XmlGspInterface::XmlGspInterface(ProjectData* project, const std::string &schemaFile)
: XMLInterface(project, schemaFile)
{
}

int XmlGspInterface::readFile(const QString &fileName)
{
	QFile* file = new QFile(fileName);
	QFileInfo fi(fileName);
	QString path = (fi.path().length() > 3) ? QString(fi.path() + "/") : fi.path();

	QFileInfo si(QString::fromStdString(_schemaName));
	QString schemaPath(si.absolutePath() + "/");

	if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "XmlGspInterface::readFile() - Can't open xml-file " <<
		fileName.toStdString() << "." << std::endl;
		delete file;
		return 0;
	}
	if (!checkHash(fileName))
	{
		delete file;
		return 0;
	}

	QDomDocument doc("OGS-PROJECT-DOM");
	doc.setContent(file);
	QDomElement docElement = doc.documentElement(); //OpenGeoSysProject
	if (docElement.nodeName().compare("OpenGeoSysProject"))
	{
		std::cout << "XmlGspInterface::readFile() - Unexpected XML root." << std::endl;
		delete file;
		return 0;
	}

	QDomNodeList fileList = docElement.childNodes();

	for(int i = 0; i < fileList.count(); i++)
	{
		if (fileList.at(i).nodeName().compare("geo") == 0)
		{
			XmlGmlInterface gml(_project, schemaPath.toStdString() + "OpenGeoSysGLI.xsd");
			gml.readFile(QString(path + fileList.at(i).toElement().text()));
		}
		else if (fileList.at(i).nodeName().compare("stn") == 0)
		{
			XmlStnInterface stn(_project, schemaPath.toStdString() + "OpenGeoSysSTN.xsd");
			QDomNodeList childList = fileList.at(i).childNodes();
			for(int j = 0; j < childList.count(); j++)
				if (childList.at(j).nodeName().compare("file") == 0)
					stn.readFile(QString(path + childList.at(j).toElement().text()));
		}
		else if (fileList.at(i).nodeName().compare("msh") == 0)
		{
			std::string msh_name = path.toStdString() +
			                       fileList.at(i).toElement().text().toStdString();
			MeshLib::CFEMesh* msh = FileIO::OGSMeshIO::loadMeshFromFile(msh_name);
			QFileInfo fi(QString::fromStdString(msh_name));
			std::string name = fi.fileName().toStdString();
			_project->addMesh(msh, name);
			//GridAdapter msh(fileList.at(i).toElement().text().toStdString());
			// TODO gridadapter to mesh-models
		}
	}

	return 1;
}

int XmlGspInterface::writeFile(const QString &fileName, const QString &tmp) const
{
	Q_UNUSED(tmp)
	GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();
	std::fstream stream(fileName.toStdString().c_str(), std::ios::out);
	QFileInfo fi(fileName);
	QString path(fi.absolutePath() + "/");
	if (!stream.is_open())
	{
		std::cout << "XmlGspInterface::writeFile() - Could not open file...\n";
		return 0;
	}

	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysProject.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-PROJECT-DOM");
	QDomElement root = doc.createElement("OpenGeoSysProject");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.net" );
	root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
	root.setAttribute( "xsi:noNamespaceSchemaLocation",
	                   "http://141.65.34.25/OpenGeoSysProject.xsd" );

	doc.appendChild(root);

	// GLI
	std::vector<std::string> geoNames;
	geoObjects->getGeometryNames(geoNames);
	for (std::vector<std::string>::const_iterator it(geoNames.begin()); it != geoNames.end();
	     ++it)
	{
		// write GLI file
		XmlGmlInterface gml(_project, path.toStdString() + "OpenGeoSysGLI.xsd");
		QString name(QString::fromStdString(*it));
		gml.writeFile(QString(path + name + ".gml"), name);

		// write entry in project file
		QDomElement geoTag = doc.createElement("geo");
		root.appendChild(geoTag);
		QDomElement fileNameTag = doc.createElement("file");
		geoTag.appendChild(fileNameTag);
		QDomText fileNameText = doc.createTextNode(QString(name + ".gml"));
		fileNameTag.appendChild(fileNameText);
	}

	// MSH
	const std::map<std::string, MeshLib::CFEMesh*> msh_vec = _project->getMeshObjects();
	for (std::map<std::string, MeshLib::CFEMesh*>::const_iterator it(msh_vec.begin());
	     it != msh_vec.end(); ++it)
	{
		// write mesh file
		QString fileName(path + QString::fromStdString(it->first));
		std::ofstream out (fileName.toStdString().c_str(), std::fstream::out);
		if (out.is_open())
		{
			FileIO::OGSMeshIO::write (it->second, out);
			out.close();
		}
		else
			std::cout << "MshTabWidget::saveMeshFile() - Could not create file..." <<
			std::endl;

		// write entry in project file
		QDomElement mshTag = doc.createElement("msh");
		root.appendChild(mshTag);
		QDomElement fileNameTag = doc.createElement("file");
		mshTag.appendChild(fileNameTag);
		QDomText fileNameText = doc.createTextNode(QString::fromStdString(it->first));
		fileNameTag.appendChild(fileNameText);
	}

	// STN
	std::vector<std::string> stnNames;
	geoObjects->getStationVectorNames(stnNames);
	for (std::vector<std::string>::const_iterator it(stnNames.begin()); it != stnNames.end();
	     ++it)
	{
		// write STN file
		XmlStnInterface stn(_project, path.toStdString() + "OpenGeoSysSTN.xsd");
		QString name(QString::fromStdString(*it));

		if (stn.writeFile(QString(path + name + ".stn"), name))
		{
			// write entry in project file
			QDomElement geoTag = doc.createElement("stn");
			root.appendChild(geoTag);
			QDomElement fileNameTag = doc.createElement("file");
			geoTag.appendChild(fileNameTag);
			QDomText fileNameText = doc.createTextNode(QString(name + ".stn"));
			fileNameTag.appendChild(fileNameText);
		}
		else
			std::cout << "XmlGspInterface::writeFile() -  Error writing file: " <<
			name.toStdString() << std::endl;
	}

	std::string xml = doc.toString().toStdString();
	stream << xml;
	stream.close();
	return 1;
}
