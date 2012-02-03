/**
 * \file XmlLutReader.h
 * 2011/01/30 KR Initial implementation
 */

#ifndef XMLLUTREADER_H
#define XMLLUTREADER_H

#include "VtkColorLookupTable.h"

#include <QFile>
#include <QtXml/QDomDocument>

#include <iostream>

/**
 * \brief Reader for vtk-Lookup-Tables (in XML / ParaView Format)
 */
class XmlLutReader
{
public:
	static VtkColorLookupTable* readFromFile(const QString &fileName)
	{
		VtkColorLookupTable* lut = VtkColorLookupTable::New();

		QFile* file = new QFile(fileName);
		if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
		{
			std::cout << "XmlLutReader::readFromFile() - Can't open xml-file." << std::endl;
			delete file;
			return NULL;
		}

		QDomDocument doc("ColorMap");
		doc.setContent(file);
		QDomElement docElement = doc.documentElement();
		if (docElement.nodeName().compare("ColorMap"))
		{
			std::cout << "XmlLutReader::readFromFile() - Unexpected XML root." << std::endl;
			delete file;
			return NULL;
		}

		if (docElement.hasAttribute("interpolation"))
		{
			if (docElement.attribute("interpolation").compare("Linear") == 0)
				lut->setInterpolationType(VtkColorLookupTable::LINEAR);
			else if (docElement.attribute("interpolation").compare("Exponential") == 0)
				lut->setInterpolationType(VtkColorLookupTable::EXPONENTIAL);
			else 
				lut->setInterpolationType(VtkColorLookupTable::NONE);
		}
		else // default
			lut->setInterpolationType(VtkColorLookupTable::NONE);

		QDomElement point = docElement.firstChildElement();
		QDomElement last_point = docElement.lastChildElement();

		double range_start(0), range(1);
		if (point.hasAttribute("x") && last_point.hasAttribute("x"))
		{
			range_start = strtod((point.attribute("x")).toStdString().c_str(),0);
			range = strtod((last_point.attribute("x")).toStdString().c_str(),0) - range_start;
		}
		else return NULL;
		
		lut->SetTableRange(range_start, range_start+range);

		while (!point.isNull())
		{
			if ((point.nodeName().compare("Point") == 0 ) 
				&& point.hasAttribute("x")
				&& point.hasAttribute("r") 
				&& point.hasAttribute("g") 
				&& point.hasAttribute("b"))
			{
				double value = strtod((point.attribute("x")).toStdString().c_str(),0);
				char r = static_cast<int>(255 * strtod((point.attribute("r")).toStdString().c_str(),0));
				char g = static_cast<int>(255 * strtod((point.attribute("g")).toStdString().c_str(),0));
				char b = static_cast<int>(255 * strtod((point.attribute("b")).toStdString().c_str(),0));
				char o = static_cast<int>(255 * (point.hasAttribute("o") ? strtod((point.attribute("o")).toStdString().c_str(),0) : 1));

				unsigned char a[4] = { r, g, b, o };
				lut->setColor(value, a);
			}
			point = point.nextSiblingElement();
		}

		delete file;

		return lut;
	};


};

#endif // XMLLUTREADER_H
