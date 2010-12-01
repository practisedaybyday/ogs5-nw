/*
 * ConvertSHPToGLI.cpp
 *
 *  Created on: May 3, 2010
 *      Author: fischeth
 */

//ShapeLib includes
#include "shapefil.h"

#include <iostream>
#include <fstream>

void convertPoints (DBFHandle dbf_handle, std::ofstream &out)
{
	if (!out) return;

	int n_records (DBFGetRecordCount (dbf_handle));
	std::cout << "writing " << n_records << " records" << std::endl;

	out << "#POINTS" << std::endl;

	for (int k(0); k<n_records; k++) {
		// *** points Kaenozoikum
		double x (DBFReadDoubleAttribute( dbf_handle, k, 0));
		double y (DBFReadDoubleAttribute( dbf_handle, k, 1));
//		double z0 (DBFReadDoubleAttribute( dbf_handle, k, 3));
		double z1 (DBFReadDoubleAttribute( dbf_handle, k, 4));
		std::string name (DBFReadStringAttribute( dbf_handle, k, 2));
		out << k << " " << x << " " << y << " " << z1 << " $NAME " << name.c_str() << std::endl;

		// *** points Keuper, Muschelkalk, Buntsandstein, Zechstein
//		double x (DBFReadDoubleAttribute( dbf_handle, k, 0));
//		double y (DBFReadDoubleAttribute( dbf_handle, k, 1));
//		double z (DBFReadDoubleAttribute( dbf_handle, k, 2));
//		std::string name (DBFReadStringAttribute( dbf_handle, k, 3));
//		out << k << " " << x << " " << y << " " << z << " $NAME " << name << std::endl;
	}
}

int main (int argc, char *argv[])
{
	if (argc == 1) return -1;

	std::string fname (argv[1]);
	std::cout << "open file " << fname.c_str() << " ... " << std::flush;

	DBFHandle dbf_handle = DBFOpen(fname.c_str(),"rb");
	if(dbf_handle) {
		std::cout << "opened" << std::endl;
		char *field_name (new char[12]);
		int width(0), n_decimals(0);
		size_t n_fields (DBFGetFieldCount(dbf_handle));
		std::cout << "file contains " << n_fields << " fields" << std::endl;

		for (size_t field_idx (0); field_idx < n_fields; field_idx++) {
			std::cout << "get field info " << field_idx << " ... " << std::flush;
			DBFGetFieldInfo( dbf_handle, field_idx, field_name, &width, &n_decimals);
			std::cout << "ok" << std::endl;
			std::cout << "DBFGetFieldInfo: " << field_idx << ", FieldName: " << field_name << ", Width: " << width << ", n_decimal: " << n_decimals << std::endl;
		}
		delete [] field_name;

		fname = fname.substr (0, fname.size()-4);
		fname += ".gli";

		std::ofstream out (fname.c_str());

		convertPoints (dbf_handle, out);

		DBFClose (dbf_handle);
	}

	return 0;
}
