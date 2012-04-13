/**
 * \file OGSFileConverter.cpp
 * 2012/04/04 KR Initial implementation
 */

#include "OGSFileConverter.h"
#include "FileListDialog.h"
#include "ConversionTools.h"
#include "OGSError.h"

// conversion includes
#include "ProjectData.h"
#include "GEOObjects.h"
#include "OGSIOVer4.h"
#include "XmlIO/XmlCndInterface.h"
#include "XmlIO/XmlGmlInterface.h"
#include "StringTools.h"

OGSFileConverter::OGSFileConverter(QWidget* parent)
	: QDialog(parent)
{
	setupUi(this);
}

OGSFileConverter::~OGSFileConverter()
{
}

void OGSFileConverter::convertGML2GLI(const QStringList input, const QString output)
{
	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects;
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder = createFileFinder();
	std::string schemaName(fileFinder.getPath("OpenGeoSysGLI.xsd"));
	FileIO::XmlGmlInterface xml(&project, schemaName);

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
		xml.readFile(*it);

	FileIO::writeAllDataToGLIFileV4(output.toStdString(), *geo_objects);
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertGLI2GML(const QStringList input, const QString output)
{
	ProjectData project;
	GEOLIB::GEOObjects* geo_objects = new GEOLIB::GEOObjects;
	project.setGEOObjects(geo_objects);

	std::vector<std::string> merge_list;
	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
	{
		std::string unique_name;
		std::vector<std::string> errors;
		if (! FileIO::readGLIFileV4(it->toStdString(), geo_objects, unique_name, errors)) 
		{
			for (size_t k(0); k<errors.size(); k++)
				OGSError::box(QString::fromStdString(errors[k]));
		}
		else
			merge_list.push_back(unique_name);
	}

	if (!merge_list.empty())
	{
		std::string merged_geo_name (merge_list[0]);
		if (merge_list.size()>1)
		{
			merged_geo_name = getFileNameFromPath(output.toStdString());
			geo_objects->mergeGeometries(merge_list, merged_geo_name);
		}
		FileFinder fileFinder = createFileFinder();
		std::string schemaName(fileFinder.getPath("OpenGeoSysGLI.xsd"));
		FileIO::XmlGmlInterface xml(&project, schemaName);
		xml.setNameForExport(merged_geo_name);
		xml.writeToFile(output.toStdString());
	}
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertCND2BC(const QStringList input, const QString output)
{
	ProjectData project;
	FileFinder fileFinder = createFileFinder();
	std::string schemaName(fileFinder.getPath("OpenGeoSysGLI.xsd"));
	FileIO::XmlCndInterface xml(&project, schemaName);

	std::vector<FEMCondition*> conditions;

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
		xml.readFile(conditions, *it);

	//now write file based on extension (bc, ic, st) and write only conditions matching that type
	OGSError::box("Not yet implemented");
}

void OGSFileConverter::convertBC2CND(const QStringList input, const QString output)
{
	ProjectData project;
	std::vector<FEMCondition*> conditions;
	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
		ConversionTools::getFEMConditionsFromASCIIFile(*it, conditions);

	if (!conditions.empty())
	{
		project.addConditions(conditions);
		FileFinder fileFinder = createFileFinder();
		std::string schemaName(fileFinder.getPath("OpenGeoSysCND.xsd"));
		FileIO::XmlCndInterface xml(&project, schemaName);
		xml.writeToFile(output.toStdString());
	}
	OGSError::box("File conversion finished");
}

FileFinder OGSFileConverter::createFileFinder()
{
	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));
	return fileFinder;
}

void OGSFileConverter::on_gml2gliButton_pressed()
{
	FileListDialog dlg(FileListDialog::GML, FileListDialog::GLI);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertGML2GLI(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_gli2gmlButton_pressed()
{
	FileListDialog dlg(FileListDialog::GLI, FileListDialog::GML);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertGLI2GML(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_bc2cndButton_pressed()
{
	FileListDialog dlg(FileListDialog::BC, FileListDialog::CND);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertBC2CND(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_cnd2bcButton_pressed()
{
	FileListDialog dlg(FileListDialog::CND, FileListDialog::BC);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertCND2BC(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_closeDialogButton_pressed()
{
	this->close();
}

