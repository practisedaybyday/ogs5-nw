/*
 * OGSIOVer4.cpp
 *
 *  Created on: Jan 14, 2010
 *      Author: TF
 */

#include <sstream>
#include <iomanip>

// FileIO
#include "OGSIOVer4.h"
#include "GMSHInterface.h"

// Base
#include "StringTools.h"
#include "quicksort.h"

// GEO
#include "GEOObjects.h"
#include "Point.h"
#include "Polyline.h"
#include "Polygon.h"
#include "SimplePolygonTree.h"
#include "Triangle.h"
#include "Surface.h"

// for tests only
#include "PointVec.h"
#include "PolylineSmoother.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "EarClippingTriangulation.h"

using namespace GEOLIB;

namespace FileIO {
/**************************************************************************
 GeoLib- Funktion: readPoints
 Aufgabe: Lesen der GLI Points und schreiben in einen Vector
 08/2005 CC Implementation
 01/2010 TF big modifications
 **************************************************************************/
/** reads the points inclusive their names from input stream in
 * using the OGS-4 file format */
std::string readPoints(std::istream &in, std::vector<Point*>* pnt_vec,
		bool &zero_based_indexing, std::map<std::string,size_t>* pnt_id_name_map)
{
	std::string line;
	size_t cnt(0);

	getline(in, line);
	// geometric key words start with the hash #
	// while not found a new key word do ...
	while (line.find("#") == std::string::npos && !in.eof() && !in.fail()) {
		// read id and point coordinates
		std::stringstream inss(line);
		size_t id;
		double x, y, z;
		inss >> id >> x >> y >> z;
		if (!inss.fail ()) {
			if (cnt == 0) {
				if (id == 0)
					zero_based_indexing = true;
				else
					zero_based_indexing = false;
			}
			pnt_vec->push_back(new Point(x, y, z));

			// read mesh density
			if (line.find("$MD") != std::string::npos) {
				double mesh_density;
				size_t pos1(line.find_first_of("M"));
				inss.str(line.substr(pos1 + 2, std::string::npos));
				inss >> mesh_density;
			}

			// read name of point
			size_t pos (line.find("$NAME"));
			if (pos != std::string::npos) { //OK
				size_t end_pos ((line.substr (pos+6)).find(" "));
				if (end_pos != std::string::npos) {
					(*pnt_id_name_map)[line.substr (pos+6, end_pos)] = id;
//					std::cout << "* name: " << line.substr (pos+6, end_pos) << ", id: " << id << std::endl;
				} else {
					(*pnt_id_name_map)[line.substr (pos+6)] = id;
//					std::cout << "name: " << line.substr (pos+6) << ", id: " << id << std::endl;
				}

			}

			size_t id_pos (line.find("$ID"));
			if (id_pos != std::string::npos) {
				std::cout << "WARNING / ERROR: found tag $ID - please use tag $NAME for reading point names" << cnt << std::endl;
			}
			cnt++;
		}
		getline(in, line);
	}

	return line;
}

/** reads points from a vector */
void readPolylinePointVector(const std::string &fname,
		std::vector<Point*>& pnt_vec, Polyline *ply, const std::string &path)
{
	// open file
	std::ifstream in((path+fname).c_str());
	if (!in)
		std::cerr << "error opening stream from " << fname << std::endl;

	double x, y, z;
	while (in) {
		in >> x >> y >> z;
		size_t pnt_id(pnt_vec.size());
		pnt_vec.push_back(new Point(x, y, z));
		ply->addPoint(pnt_id);
	}
}

/**************************************************************************
 GeoLib-Method: Read
 Task: Read polyline data from file
 Programing:
 03/2004 CC Implementation
 09/2004 OK file path for PLY files
 07/2005 CC PLY id
 08/2005 CC parameter
 09/2005 CC itoa - convert integer to string
 01/2010 TF cleaned method from unused variables
 **************************************************************************/
/** read a single Polyline from stream in into the ply_vec-vector */
std::string readPolyline(std::istream &in,
		std::vector<Polyline*>* ply_vec, std::map<std::string,size_t>& ply_vec_names,
		std::vector<Point*>& pnt_vec, bool zero_based_indexing, const std::vector<size_t>& pnt_id_map,
		const std::string &path)
{
	std::string line, name_of_ply;
	Polyline *ply(new Polyline(pnt_vec));
	size_t type=2; // need an initial value

	// Schleife ueber alle Phasen bzw. Komponenten
	do {
		in >> line;
		if (line.find("$ID") != std::string::npos) { // subkeyword found CC
			in >> line; // read value
			//			id = strtol(line_string.data(), NULL, 0);
		}
		//....................................................................
		if (line.find("$NAME") != std::string::npos) { // subkeyword found
			in >> line;
			name_of_ply = line.substr(0); // read value
		}
		//....................................................................
		if (line.find("$TYPE") != std::string::npos) { // subkeyword found
			in >> line; // read value
			type = static_cast<size_t> (strtol(line.c_str(), NULL, 0));
		}
		//....................................................................
		if (line.find("$EPSILON") != std::string::npos) { // subkeyword found
			in >> line; // read value
		}
		//....................................................................
		if (line.find("$MAT_GROUP") != std::string::npos) { // subkeyword found
			in >> line; // read value
		}
		//....................................................................
		if (line.find("$POINTS") != std::string::npos) { // subkeyword found
			// read the point ids
			in >> line;
			if (type != 100) {
				while (!in.eof() && line.size() != 0 && (line.find("#")
						== std::string::npos) && (line.find("$")
						== std::string::npos)) {
					size_t pnt_id(str2number<size_t> (line));
					if (!zero_based_indexing)
						pnt_id--; // one based indexing
					size_t ply_size (ply->getNumberOfPoints());
					if (ply_size > 0) {
						if (ply->getPointID (ply_size-1) != pnt_id_map[pnt_id])
							ply->addPoint(pnt_id_map[pnt_id]);
					} else {
						ply->addPoint(pnt_id_map[pnt_id]);
					}
					in >> line;
				}
			} else {
				std::cerr << "*** polyline is an arc *** reading not implemented" << std::endl;
			}
			// empty line or the keyword or subkeyword or end of file
		}
		//....................................................................
		if (line.find("$POINT_VECTOR") != std::string::npos) { // subkeyword found
			in >> line; // read file name
			line = path + line;
			readPolylinePointVector(line, pnt_vec, ply, path);
		} // subkeyword found
	} while (line.find("#") == std::string::npos && line.size() != 0 && in);

	if (type != 100) {
		ply_vec_names.insert (std::pair<std::string,size_t>(name_of_ply, ply_vec->size()));
		ply_vec->push_back(ply);
	}

	return line;
}

/**************************************************************************
 GEOLib-Function:
 Task: polyline read function
 Programming:
 03/2004 CC Implementation
 05/2004 CC Modification
 04/2005 CC Modification calculate the minimal distance between points reference for
 mesh density of line element calculation
 07/2005 CC read ID of polyline
 08/2005 CC parameter
 01/2010 TF changed signature of function
 **************************************************************************/
/** reads polylines */
std::string readPolylines(std::istream &in, std::vector<Polyline*>* ply_vec,
		std::map<std::string,size_t>& ply_vec_names, std::vector<Point*>& pnt_vec,
		bool zero_based_indexing, const std::vector<size_t>& pnt_id_map,
		const std::string &path)
{
	if (!in) {
		std::cerr << "*** readPolylines input stream error " << std::endl;
		return std::string("");
	}
	std::string tag("#POLYLINE");

	while (!in.eof() && tag.find("#POLYLINE") != std::string::npos) {
		tag = readPolyline(in, ply_vec, ply_vec_names, pnt_vec,
				zero_based_indexing, pnt_id_map, path);
	}

	return tag;
}

void readTINFile(const std::string &fname, Surface* sfc,
		std::vector<Point*> &pnt_vec)
{
	// open file
	std::ifstream in(fname.c_str());
	if (!in)
		std::cerr << "readTINFile error opening stream from " << fname
				<< std::endl;

	size_t id;
	double x, y, z;
	while (in) {
		// read id
		in >> id;
		// determine size
		size_t pnt_pos(pnt_vec.size());
		// read first point
		in >> x >> y >> z;
		pnt_vec.push_back(new Point(x, y, z));
		// read second point
		in >> x >> y >> z;
		pnt_vec.push_back(new Point(x, y, z));
		// read third point
		in >> x >> y >> z;
		pnt_vec.push_back(new Point(x, y, z));
		// create new Triangle
		sfc->addTriangle(pnt_pos, pnt_pos+1, pnt_pos+2);
	}
}

/**************************************************************************
 GeoLib-Method: readSurface
 Task: Read surface data from input stream
 Programing:
 03/2004 OK Implementation
 05/2005 OK EPSILON
 09/2005 CC file_path_base
 01/2010 TF signatur modification, reimplementation
 **************************************************************************/
/** read a single Surface */
std::string readSurface(std::istream &in,
		std::vector<Polygon*> &polygon_vec,
		std::vector<Surface*> &sfc_vec, std::map<std::string,size_t>& sfc_names,
		const std::vector<Polyline*> &ply_vec, const std::map<std::string, size_t>& ply_vec_names,
		std::vector<Point*> &pnt_vec,
		const std::string &path)
{
	std::string line;
	Surface *sfc(NULL);

	int type (-1);
	std::string name;
	size_t ply_id (0); // std::numeric_limits<size_t>::max());

	do {
		in >> line;
		if (line.find("$ID") != std::string::npos) { // subkeyword found CC
			in >> line; // read value
			//			id = strtol(line_string.data(), NULL, 0);
		}
		//....................................................................
		if (line.find("$NAME") != std::string::npos) { // subkeyword found
			in >> line; // read value
			name = line.substr(0);
		}
		//....................................................................
		if (line.find("$TYPE") != std::string::npos) { // subkeyword found
			in >> line; // read value
			type = strtol(line.c_str(), NULL, 0);
		}
		//....................................................................
		if (line.find("$EPSILON") != std::string::npos) { // subkeyword found
			in >> line; // read value
		}
		//....................................................................
		if (line.find("$TIN") != std::string::npos) { // subkeyword found
			in >> line; // read value (file name)
			line = path + line;
//			if (type == 1) std::cerr << "reading tin file " << line << " ... " << std::flush;
			sfc = new Surface(pnt_vec);
			readTINFile(line, sfc, pnt_vec);
//			std::cout << "ok" << std::endl;
		}
		//....................................................................
		if (line.find("$MAT_GROUP") != std::string::npos) { // subkeyword found
			in >> line; // read value
		}
		//....................................................................
		if (line.find("$POLYLINES") != std::string::npos) { // subkeyword found
			// read the name of the polyline(s)
			in >> line;
			while (!in.eof() && line.size() != 0
					&& (line.find("#")== std::string::npos)
					&& (line.find("$") == std::string::npos)) {
				// we did read the name of a polyline -> search the id for polyline
				std::map<std::string,size_t>::const_iterator it (ply_vec_names.find (line));
				if (it != ply_vec_names.end())
					ply_id = it->second;
				else
					ply_id = ply_vec.size();

				if (ply_id == ply_vec.size())
					std::cerr << "polyline for surface not found!" << std::endl;
				else {
					if (type == 3) {
						std::cerr << "surface type 3: flat surface with any normal direction - - reading not implemented"
								<< std::endl;
					}

					if (type == 2) {
						std::cerr << "vertical surface - reading not implemented"
												<< std::endl;
					}
//					if (type == 0 || type == -1)
//						std::cerr << "reading Polygon " << ply_vec_names[ply_id] << std::endl;
					else std::cerr << "unknown polyline type " << type << "... " << std::endl;
				}
				in >> line;
			}
			// empty line or a keyword is found
		}
	} while (line.find("#") == std::string::npos && line.size() != 0 && in);

	if (!name.empty())
		sfc_names.insert(std::pair<std::string,size_t>(name,sfc_vec.size()));

	if (sfc) {
		// surface create by TIN
		sfc_vec.push_back (sfc);
	} else {
		// surface created by polygon
		if (ply_id != std::numeric_limits<size_t>::max()) {
			if (ply_vec[ply_id]->isClosed()) {
				polygon_vec.push_back (new Polygon (*(ply_vec[ply_id]), true));
			}
		}
	}

	return line;
}

/**************************************************************************
 GEOLib-Method:
 Task: Surface read function
 Programming:
 03/2004 OK Implementation
 05/2004 CC Modification
 01/2010 TF changed signature of function, big modifications
 **************************************************************************/
std::string readSurfaces(std::istream &in,
		std::vector<Surface*> &sfc_vec, std::map<std::string, size_t>& sfc_names,
		const std::vector<Polyline*> &ply_vec, const std::map<std::string,size_t>& ply_vec_names,
		std::vector<Point*> &pnt_vec,
		const std::string &path)
{
	if (!in.good()) {
		std::cerr << "*** readSurfaces input stream error " << std::endl;
		return std::string("");
	}
	std::string tag("#SURFACE");

	std::vector<Polygon*> polygon_vec;

	while (!in.eof() && tag.find("#SURFACE") != std::string::npos) {
		size_t n_polygons (polygon_vec.size());
		tag = readSurface(in, polygon_vec, sfc_vec, sfc_names, ply_vec, ply_vec_names, pnt_vec, path);
		if (n_polygons < polygon_vec.size()) {
			// subdivide polygon in simple polygons
			GEOLIB::Surface *sfc (GEOLIB::Surface::createSurface (*(dynamic_cast<GEOLIB::Polyline*>(polygon_vec[polygon_vec.size()-1]))));
			sfc_vec.push_back (sfc);
		}
	}
	for (size_t k(0); k<polygon_vec.size(); k++) {
		delete polygon_vec[k];
	}
	std::cout << "readSurfaces: number of read polygons  " << polygon_vec.size() << std::endl;

//	// subdivide all polygons in simple polygons
//	for (std::vector<GEOLIB::Polygon*>::iterator polygon_it (polygon_vec.begin());
//			polygon_it != polygon_vec.end(); polygon_it++) {
//		// compute list of simple polygons
//		(*polygon_it)->computeListOfSimplePolygons ();
//	}

//	// forest consist of (hierarchy) trees
//	std::list<SimplePolygonTree*> polygon_forest;
//	// create polygon forest
//	for (std::vector<GEOLIB::Polygon*>::iterator polygon_it (polygon_vec.begin());
//				polygon_it != polygon_vec.end(); polygon_it++) {
//		// get the list and insert the elements as SimplePolygonTree items into the forest
//		const std::list<Polygon*> simple_polygon_list ((*polygon_it)->getListOfSimplePolygons());
//		for (std::list<Polygon*>::const_iterator simple_polygon_it (simple_polygon_list.begin());
//			simple_polygon_it != simple_polygon_list.end(); simple_polygon_it++) {
//			SimplePolygonTree *spt (new SimplePolygonTree (*simple_polygon_it));
//			polygon_forest.push_back (spt);
//		}
//	}
//	std::cout << "readSurfaces: \"Polygon forest\" consists of " << polygon_forest.size() << " trees" << std::endl;
//
//	// create the hierarchy
//	createPolygonTree (polygon_forest);
//	std::cout << "readSurfaces: \"Polygon forest\" consists of " << polygon_forest.size() << " trees" << std::endl;
//
//	std::string out_fname ("GMSHTest.geo");
//	std::cout << "writing input file for GMSH " << out_fname << " ... " << std::flush;
//	GMSHInterface gmsh_io (out_fname);
//	// writing points
//	gmsh_io.writeGMSHPoints(pnt_vec);
//	// writing simple polygon tree
//	for (std::list<GEOLIB::SimplePolygonTree*>::const_iterator polygon_tree_it (polygon_forest.begin());
//	polygon_tree_it != polygon_forest.end(); polygon_tree_it++) {
//		(*polygon_tree_it)->visitAndProcessNodes (gmsh_io);
//	}
//	std::cout << "done" << std::endl;

	// create surfaces from simple polygons
//	for (std::vector<GEOLIB::Polygon*>::iterator polygon_it (polygon_vec.begin());
//		polygon_it != polygon_vec.end(); polygon_it++) {
//
//		const std::list<GEOLIB::Polygon*>& list_of_simple_polygons ((*polygon_it)->getListOfSimplePolygons());
//
//		for (std::list<GEOLIB::Polygon*>::const_iterator simple_polygon_it (list_of_simple_polygons.begin());
//			simple_polygon_it != list_of_simple_polygons.end(); simple_polygon_it++) {
//
//			std::list<GEOLIB::Triangle> triangles;
//			MathLib::earClippingTriangulationOfPolygon(*simple_polygon_it, triangles);
//			std::cout << "done - " << triangles.size () << " triangles " << std::endl;
//
//			Surface *sfc(new Surface(pnt_vec));
//			// add Triangles to Surface
//			std::list<GEOLIB::Triangle>::const_iterator it (triangles.begin());
//			while (it != triangles.end()) {
//				sfc->addTriangle ((*it)[0], (*it)[1], (*it)[2]);
//				it++;
//			}
//			sfc_vec.push_back (sfc);
//		}
//	}

	return tag;
}

void readGLIFileV4(const std::string& fname, GEOObjects* geo)
{
	std::cout << "GEOLIB::readGLIFile open stream from file " << fname
			<< " ... " << std::flush;
	std::ifstream in(fname.c_str());
	if (!in)
		std::cerr << "error opening stream from " << fname << std::endl;
	std::cout << "done" << std::endl;

	std::string tag;
	while (tag.find("#POINTS") == std::string::npos && !in.eof()) {
		getline (in, tag);
	}

	// read names of points into vector of strings
	std::map<std::string,size_t>* pnt_id_names_map (new std::map<std::string,size_t>);
	bool zero_based_idx(true);
	std::vector<Point*> *pnt_vec(new std::vector<Point*>);
	std::cout << "read points from stream ... " << std::flush;
	tag = readPoints(in, pnt_vec, zero_based_idx, pnt_id_names_map);
	std::cout << " ok, " << pnt_vec->size() << " points read" << std::endl;

	std::string unique_name(fname);
	if (!pnt_vec->empty())
		geo->addPointVec(pnt_vec, unique_name, pnt_id_names_map); // KR: insert into GEOObjects if not empty

	// extract path for reading external files
	size_t pos(fname.rfind("/"));
	std::string path(fname.substr(0, pos + 1));

	// read names of plys into temporary string-vec
	std::map<std::string,size_t> *ply_names (new std::map<std::string,size_t>);
	std::vector<Polyline*> *ply_vec(new std::vector<Polyline*>);
	if (tag.find("#POLYLINE") != std::string::npos && in) {
		std::cout << "read polylines from stream ... " << std::flush;
		tag = readPolylines(in, ply_vec, *ply_names, *pnt_vec,
				zero_based_idx, geo->getPointVecObj(unique_name)->getIDMap(), path);
		std::cout << " ok, " << ply_vec->size() << " polylines read"
				<< std::endl;
	} else
		std::cerr
				<< "tag #POLYLINE not found or input stream error in GEOObjects"
				<< std::endl;

	std::vector<Surface*> *sfc_vec(new std::vector<Surface*>);
	std::map<std::string,size_t> *sfc_names (new std::map<std::string,size_t>);
	if (tag.find("#SURFACE") != std::string::npos && in) {
		std::cout << "read surfaces from stream ... " << std::flush;
		tag = readSurfaces(in, *sfc_vec, *sfc_names, *ply_vec, *ply_names, *pnt_vec, path);
		std::cout << " ok, " << sfc_vec->size() << " surfaces read"
				<< std::endl;
	} else
		std::cerr
				<< "tag #SURFACE not found or input stream error in GEOObjects"
				<< std::endl;
	in.close();

	if (!ply_vec->empty())
		geo->addPolylineVec(ply_vec, unique_name, ply_names); // KR: insert into GEOObjects if not empty
	if (!sfc_vec->empty())
		geo->addSurfaceVec(sfc_vec, unique_name, sfc_names); // KR: insert into GEOObjects if not empty
}

void writeGLIFileV4 (const std::string& fname, const std::string& geo_name, const GEOLIB::GEOObjects& geo)
{
	const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec(geo_name));
	std::ofstream os (fname.c_str());
	if (pnts) {
		std::cout << "writing " << pnts->size () << " points to file " << fname << std::endl;
		os << "#POINTS" << std::endl;
		os.precision (20);
		for (size_t k(0); k<pnts->size(); k++) {
			os << k << " " << *((*pnts)[k]) << std::endl;
		}
	}

	std::cout << "writing " << std::flush;
	const GEOLIB::PolylineVec* plys_vec (geo.getPolylineVecObj (geo_name));
	if (plys_vec) {
		const std::vector<GEOLIB::Polyline*>* plys (plys_vec->getVector());
		std::cout << plys->size () << " polylines to file " << fname << std::endl;
//		const std::vector<size_t>& pnt_id_map (geo.getPointVecObj(geo_name)->getIDMap());
		for (size_t k(0); k<plys->size(); k++) {
			os << "#POLYLINE" << std::endl;
			os << " $NAME " << std::endl << "  " << k << std::endl; // plys_vec->getNameOfElement ((*plys)[k]) << std::endl;
			os << " $POINTS" << std::endl;
			for (size_t j(0); j<(*plys)[k]->getNumberOfPoints(); j++) {
//				os << "  " << pnt_id_map[((*plys)[k])->getPointID(j)] << std::endl;
				os << "  " << ((*plys)[k])->getPointID(j) << std::endl;
			}
		}
	}

	std::cout << "writing " << std::flush;
	if (plys_vec) {
		const std::vector<GEOLIB::Polyline*>* plys (plys_vec->getVector());
		std::cout << plys->size () << " closed polylines as surfaces to file " << fname << std::endl;
		for (size_t k(0); k<plys->size(); k++) {
			if ((*plys)[k]->isClosed()) {
				os << "#SURFACE" << std::endl;
				os << " $NAME " << std::endl << "  " << k << std::endl; //plys_vec->getNameOfElement ((*plys)[k]) << std::endl;
				os << " $TYPE " << std::endl << "  0" << std::endl;
				os << " $POLYLINES" << std::endl << "  " << k << std::endl; //plys_vec->getNameOfElement ((*plys)[k]) << std::endl;
			}
		}
	}

	os << "#STOP" << std::endl;
	os.close ();
}

void writeAllDataToGLIFileV4 (const std::string& fname, const GEOLIB::GEOObjects& geo)
{
	std::vector<std::string> geo_names;
	geo.getGeometryNames (geo_names);

	std::ofstream os (fname.c_str());

	size_t pnts_offset (0);
	std::vector<size_t> pnts_id_offset;
	pnts_id_offset.push_back (0);

	// writing all points
	os << "#POINTS" << std::endl;
	for (size_t j(0); j<geo_names.size(); j++) {
		os.precision (20);
		const std::vector<GEOLIB::Point*>* pnts (geo.getPointVec(geo_names[j]));
		if (pnts) {
			for (size_t k(0); k<pnts->size(); k++) {
				os << k+pnts_offset << " " << *((*pnts)[k]) << std::endl;
			}
			pnts_offset += pnts->size();
			pnts_id_offset.push_back (pnts_offset);
		}
	}

	std::cout << "wrote " << pnts_offset << " points" << std::endl;

	// writing all stations
	std::vector<std::string> stn_names;
	geo.getStationNames (stn_names);
	for (size_t j(0); j<stn_names.size(); j++) {
		os.precision (20);
		const std::vector<GEOLIB::Point*>* pnts (geo.getStationVec(stn_names[j]));
		if (pnts) {
			for (size_t k(0); k<pnts->size(); k++) {
				os << k+pnts_offset << " " << *((*pnts)[k]) << " $NAME " << static_cast<GEOLIB::Station*>((*pnts)[k])->getName() << std::endl;
			}
			pnts_offset += pnts->size();
			pnts_id_offset.push_back (pnts_offset);
		}
	}

	size_t plys_cnt (0);

	// writing all polylines
	for (size_t j(0); j<geo_names.size(); j++) {
		const GEOLIB::PolylineVec* plys_vec (geo.getPolylineVecObj (geo_names[j]));
		if (plys_vec) {
			const std::vector<GEOLIB::Polyline*>* plys (plys_vec->getVector());
			for (size_t k(0); k<plys->size(); k++) {
				os << "#POLYLINE" << std::endl;
				std::string ply_name;
				if (plys_vec->getNameOfElementByID (plys_cnt, ply_name)) {
					os << " $NAME " << std::endl << "  " << ply_name << std::endl;
				} else {
					os << " $NAME " << std::endl << "  " << plys_cnt << std::endl;
				}
				os << " $POINTS" << std::endl;
				for (size_t l(0); l<(*plys)[k]->getNumberOfPoints(); l++) {
					os << "  " << pnts_id_offset[j] + ((*plys)[k])->getPointID(l) << std::endl;
				}
				plys_cnt++;
			}
		}
	}

	os << "#STOP" << std::endl;
	os.close ();
}

} // end namespace
