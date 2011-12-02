//1. Eclipse starten (die Files liegen alle vorbereitet vor
//2. Daten auslesen -> Tool ecl2gs
//3. Datenumwandeln (�bertragung auf Knoten, Geschwindigkeit berechnen)
//4. Geosys daten schreiben
//5. Berechnung Geosys
//6. Eclipse restart

// declarations for ecl_interface ECLIPSE - GeoSys

//#include "stdafx.h"
#include "Eclipse.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>  // Datei streams
#include <iostream>  // Bildschirmausgabe
#include <sstream>
#include <sstream>   // string streams (in)
#include <string>
#include <sys/types.h>
#include <vector>
//#include "dirent.h"
#include "display.h"
#include <errno.h>
#include <math.h>
//#include "Windows.h"
#include "mathlib.h"
#include <cfloat>
#include <ctime>
#include <sys/stat.h> //for check if files exist

using namespace std;

CECLIPSEBlock::CECLIPSEBlock(long Nodelength, long Facelength)
{
	(void)Facelength; // unused
	this->row = 0;
	this->column = 0;
	this->index = 0;
	this->layer = 0;
	this->x_coordinates.resize(Nodelength);
	this->y_coordinates.resize(Nodelength);
	this->z_coordinates.resize(Nodelength);
	this->x_barycentre = 0;
	this->y_barycentre = 0;
	this->z_barycentre = 0;
	this->active = 0;
	this->NeighbourElement.resize(6);
	this->volume = 0;
}

CECLIPSEBlock::~CECLIPSEBlock(void)
{
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CReadTextfiles
   Task: Read textfiles in a vector of strings
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
CReadTextfiles_ECL::CReadTextfiles_ECL(void)
{
	this->NumberOfRows = 0;
	this->Data.clear();
}

CReadTextfiles_ECL::~CReadTextfiles_ECL(void)
{
	this->Data.clear();
}

bool CReadTextfiles_ECL::Read_Text(std::string Filename)
{
	char Line[MAX_ZEILEN];
	bool error = false;
	bool abort = false;
	std::string tempstring;

	// .data ... provides the filename as it is necessary for C, ios::in ... reads the file
	ifstream in_file (Filename.data(),ios::in);
	if (in_file == 0)
	{
		error = true;
		cout << "The file " << Filename.data() << " can not be opened!" << endl;
	}
	else
	{
		error = false;
		while ((abort == false) && (in_file.eof() == false))
		{
			tempstring.clear();
			in_file.getline(Line,MAX_ZEILEN);
			tempstring = Line;
			// You could basically use this later to avoid line lenghts of pre-defined length only
			// getline(in_file,tempstring);
			if(tempstring.length() == MAX_ZEILEN - 1)
			{
				cout <<
				" Error - increase MAX_ZEILEN in order to read ECLIPSE data file "
				     << endl;
				cout << " Or shorten the line in " << Filename.data() << ": " <<
				tempstring << " to " << MAX_ZEILEN << " characters " << endl;
				exit(0);
			}
			if (tempstring.compare("#STOP") != 0)
			{
				this->Data.push_back(Line);
				this->NumberOfRows = this->NumberOfRows + 1;
			}
			else
				abort = true;
		}
	}
	return error;
}

CWriteTextfiles_ECL::CWriteTextfiles_ECL(void)
{
}

CWriteTextfiles_ECL::~CWriteTextfiles_ECL(void)
{
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Write_Text
   Task: Writes a vector of strings to a text file
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CWriteTextfiles_ECL::Write_Text(std::string Filename, vector<std::string> Text)
{
	ofstream textfile;
	textfile.open(Filename.data(), ios::out);

	for (int i = 0; i < int(Text.size()); i++)
	{
		if (i == int(Text.size()) - 1)
			textfile << Text[i].data();
		else
			textfile << Text[i].data() << endl;
		//file.write(Text[i]);
	}
	textfile.close();
}

/*-------------------------------------------------------------------------
   Constructor and Destructor of the class CECLIPSEData
   -------------------------------------------------------------------------*/
CECLIPSEData::CECLIPSEData(void)
{
	this->rows = 0;
	this->columns = 0;
	this->elements = 0;
	this->layers = 0;
	this->times = 0;
	this->numberOutputParameters = 0;
	this->activeCells = 0;
	this->NodeData.clear(); //resize(1); // warum das?
	this->Data = NULL;
	this->ProcessIndex_CO2inLiquid = -1;
	this->Radial_I = false;
	this->Radial_J = false;
	this->RadialModellIpos = false;
	this->RadialModellJpos = false;
}

CECLIPSEData::~CECLIPSEData(void)
{
}
/*-------------------------------------------------------------------------
   GeoSys - Function: GetVariableIndex
   Task: Returns the index of the given variable
   Return: Variableindex
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
int CECLIPSEData::GetVariableIndex(std::string Variablename)
{
	int index_Variable = -1;

	for (unsigned int i = 0; i < this->Variables.size(); i++)
		if (this->Variables[i].compare(Variablename) == 0)
			index_Variable = i;
	return index_Variable;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: SplitStrings
   Task: Separate a string with a delimiter
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::SplitStrings(const std::string str, std::string delimiter)
{
	this->SplittedString.clear();

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiter, lastPos);

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		this->SplittedString.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter, lastPos);
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Round
   Task: Round numbers
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
double CECLIPSEData::Round(double Number, int Decimalplaces)
{
	Number *= pow( 10, double(Decimalplaces));
	Number = floor(Number + 0.5);
	Number *= pow(10, double(-Decimalplaces));
	return Number;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Round
   Task: Round numbers
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
typeExponentialNumber CECLIPSEData::RoundEXP(double Number, int Decimalplaces)
{
	typeExponentialNumber Result;
	int Exponent;
	double tempNumber;
	//WW int sign;

	//WW sign = 1;
	//WW if (Number < 0)
	//WW	sign = -1;

	Exponent = 0;
	tempNumber = fabs(Number);
	if (tempNumber != 0)
	{
		if (tempNumber >= 1)
		{
			do {
				tempNumber /= 10;
				Exponent -= 1;
			} while(tempNumber >= 1);
		}
		else if (tempNumber < 0.1)
		{
			do {
				tempNumber *= 10;
				Exponent += 1;
			} while(tempNumber < 0.1);
		}
		Number *= pow( 10, double(Exponent + Decimalplaces));
		Number = floor(Number + 0.5);
		Number *= pow(10, double(-Decimalplaces));
	}
	else
		Exponent = -0;
	Result.Exponent = -Exponent;
	Result.Number = Number;
	return Result;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: AddZero
   Task: Adds zero before or after the given number
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
std::string CECLIPSEData::AddZero(double Number, int Places, bool before)
{
	std::string tempstring;
	std::string Result;
	stringstream NumberString, tempNumberString;
	size_t position, position_e;
	double tempNumber;
	int precision;

	precision = Places + 5;
	NumberString.precision(precision);
	NumberString << Number;
	tempstring = NumberString.str();
	Result = tempstring;

	if (before == true)
		while (int(Result.length()) < Places)
		{
			if (Result.substr(0, 1) == "-")
				Result = "-0" + Result.substr(1, Result.length());
			else
				Result = "0" + Result;
		}
	else
	{
		position = Result.find(".");
		// if the string is longer than the string is cutted and rounded up
		if (Result.length() >= Places + position + 1)
		{
			tempNumber = this->Round(Number, Places);
			tempstring = "";
			precision = Places + position;
			tempNumberString.precision(precision);
			tempNumberString << tempNumber;
			tempstring = tempNumberString.str();
			Result = tempstring;
		}
		//else {
		while (size_t(Result.length()) < Places + position + 1)
		{
			position = Result.find(".");
			if (position < Result.length())
			{
				position_e = Result.find("e");
				if (position_e < Result.length())
					Result = Result.substr(0, position_e) + "0" + Result.substr(
					        position_e,
					        Result.length() - (position_e));
				else
					Result = Result + "0";
			}
			else
			{
				position_e = Result.find("e");
				if (position_e < Result.length())
					Result = Result.substr(0, position_e) + "." + Result.substr(
					        position_e,
					        Result.length() - (position_e));
				else
					Result = Result + ".";
			}
		}
		//}
	}
	return Result;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CheckIfFileExists
   Task: Check if the file exists
   Return: true if exist
   Programming: 11/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CheckIfFileExists(string strFilename)
{
// code source: http://www.techbytes.ca/techbyte103.html (6.11.2009) no restriction for use

	struct stat stFileInfo;
	bool blnReturn;
	int intStat;

	// Attempt to get the file attributes
	intStat = stat(strFilename.c_str(),&stFileInfo);
	if(intStat == 0)
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	else
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;

	return blnReturn;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReadEclipseGrid
   Task: Reads the Eclipse model grid
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadEclipseGrid(std::string Filename)
{
	std::string tempstring;
	CECLIPSEBlock* m_eclipseblock = NULL;
	CReadTextfiles_ECL* TextFile;
	int Corners_row = 0;
	std::stringstream in;
	bool error = false;
	vector<std::string> files = vector<std::string>();
	clock_t start,finish;
	double time;

	start = clock();

	cout << "        ReadEclipseGrid()";

	TextFile = new CReadTextfiles_ECL;
	error = TextFile->Read_Text(Filename);

	if (error == true)
	{
		cout << "The program is canceled" << endl;
		//system("Pause");
		exit(0);
	}

	//Read in the Grid
	for (long i = 0; i < TextFile->NumberOfRows - 1; i++)
	{
		//cout << TextFile->Data[i] << endl;
		tempstring = TextFile->Data[i].substr(2,8);

		//Check if the grid is radial or not
		if (tempstring.compare("RADIAL  ") == 0)
		{
			if (TextFile->Data[i + 1].substr(2,8) == "FALSE   ")
				Radialmodell = false;
			else if (TextFile->Data[i + 1].substr(2,8) == "TRUE    ")
				Radialmodell = true;
			else
			{
				cout <<
				"The phrase below the keyword RADIAL does not contain TRUE or FALSE! The program is aborted!"
				     << endl;
				system("Pause");
				exit(0);
			}
		}

		if (tempstring.compare("COORDS  ") == 0)
		{
			m_eclipseblock = new CECLIPSEBlock(8, 6); // new member of eclblock
			in.str((string) TextFile->Data[i + 1]);
			in >> m_eclipseblock->column >> m_eclipseblock->row >>
			m_eclipseblock->layer >> m_eclipseblock->index >> m_eclipseblock->active;
			in.clear();
			// SB und KB

			//cout << m_eclipseblock->column << endl;

			//this->SplitStrings(TextFile->Daten[i+1]," ");
			//if (this->SplittedString.size()<=4) {
			//	cout << endl << "The grid output format (GRIDFILE) for Eclipse should be set to 2 0 instead of 1 1! The run is terminated now!" << endl;
			//	system("Pause");
			//	exit(0);
			//}
			//m_eclipseblock->column = atoi(this->SplittedString[0].data());
			//m_eclipseblock->row = atoi(this->SplittedString[1].data());
			//m_eclipseblock->layer = atoi(this->SplittedString[2].data());
			//m_eclipseblock->index = atoi(this->SplittedString[3].data());
			//m_eclipseblock->active = atoi(this->SplittedString[4].data());
			//index of Eclipse start with 1 -> change to 0
			//m_eclipseblock->index = m_eclipseblock->index - 1;

			// inactive cells are not red anymore
			if (m_eclipseblock->active == 1)
			{
				if (this->columns < m_eclipseblock->column)
					this->columns = m_eclipseblock->column;
				if (this->rows < m_eclipseblock->row)
					this->rows = m_eclipseblock->row;
				// layer number increases with depth -> layer 1 = top layer
				if (this->layers < m_eclipseblock->layer)
					this->layers = m_eclipseblock->layer;
				this->elements += 1;
				if (m_eclipseblock->active == 1)
					this->activeCells += 1;

				m_eclipseblock->index = this->elements - 1;

				//read in alle corner coordinates
				//Debug.Print(Strings.Mid(Gridfile.Daten(0, i + 2), 3, 8))

				if (TextFile->Data[i + 2].substr(2,8).compare("CORNERS ") == 0)
					Corners_row = 0;
				else
				{
					if (TextFile->Data[i + 3].substr(2,
					                                 8).compare("CORNERS ") ==
					    0)
						Corners_row = 1;
					else
						cout <<
						"Error in the file structure of the grid file!" <<
						endl;
				}
				//Reads the corner points; order of the points (x1<x2, y1<y2, z1<z2):
				//	0..x1, y1, z2; 1..x2, y1, z2; 2..x1, y2, z2; 3..x2, y2, z2
				//	4..x1, y1, z1; 5..x2, y1, z1; 6..x1, y2, z1; 7..x2, y2, z1
				in.str((string) TextFile->Data[i + 3 + Corners_row]);
				in >> m_eclipseblock->x_coordinates[0] >>
				m_eclipseblock->y_coordinates[0] >>
				m_eclipseblock->z_coordinates[0] >>
				m_eclipseblock->x_coordinates[1];
				in.clear();

				in.str((string) TextFile->Data[i + 4 + Corners_row]);
				in >> m_eclipseblock->y_coordinates[1] >>
				m_eclipseblock->z_coordinates[1] >>
				m_eclipseblock->x_coordinates[2] >>
				m_eclipseblock->y_coordinates[2];
				in.clear();

				in.str((string) TextFile->Data[i + 5 + Corners_row]);
				in >> m_eclipseblock->z_coordinates[2] >>
				m_eclipseblock->x_coordinates[3] >>
				m_eclipseblock->y_coordinates[3] >>
				m_eclipseblock->z_coordinates[3];
				in.clear();

				in.str((string) TextFile->Data[i + 6 + Corners_row]);
				in >> m_eclipseblock->x_coordinates[4] >>
				m_eclipseblock->y_coordinates[4] >>
				m_eclipseblock->z_coordinates[4] >>
				m_eclipseblock->x_coordinates[5];
				in.clear();

				in.str((string) TextFile->Data[i + 7 + Corners_row]);
				in >> m_eclipseblock->y_coordinates[5] >>
				m_eclipseblock->z_coordinates[5] >>
				m_eclipseblock->x_coordinates[6] >>
				m_eclipseblock->y_coordinates[6];
				in.clear();

				in.str((string) TextFile->Data[i + 8 + Corners_row]);
				in >> m_eclipseblock->z_coordinates[6] >>
				m_eclipseblock->x_coordinates[7] >>
				m_eclipseblock->y_coordinates[7] >>
				m_eclipseblock->z_coordinates[7];
				in.clear();

				//m_eclipseblock->CalcBarycentre(); // For block to node interpolation
				//m_eclipseblock->CalculateFaceCentres();

				//cout << m_eclipseblock->x_coordinates[0] << " " << m_eclipseblock->x_coordinates[1] << " " << m_eclipseblock->x_coordinates[2] << endl;
				// adds 1 dataset of typ eclblocks to vec_eclblocks
				this->eclgrid.push_back(m_eclipseblock);
			}
		}
	}

	//recalculate coordinates if it is a radial grid based on angle and radius
	double alpha = 0;
	double radius = 0;
	if (Radialmodell == true)
	{
		for (long i = 0; i < this->elements; i++)
			for (int j = 0; j < int(this->eclgrid[i]->x_coordinates.size()); j++)
			{
				//coordinates are defined as radius, angle (clockwise in degree), heigth
				//recalculate alpha from degree to radian (Bogenma�)
				//cout << "Element: " << i << " Point: " << j << " alpha: " << this->eclgrid[i]->y_coordinates[j] << " radius: " << this->eclgrid[i]->x_coordinates[j];
				alpha = this->eclgrid[i]->y_coordinates[j] * PI / 180;
				radius = this->eclgrid[i]->x_coordinates[j];
				this->eclgrid[i]->x_coordinates[j] = sin(alpha) * radius;
				this->eclgrid[i]->y_coordinates[j] = cos(alpha) * radius;
				//cout << " radius: " << radius << " x: " << this->eclgrid[i]->x_coordinates[j] << " y: " << this->eclgrid[i]->y_coordinates[j] << endl;
			}
		//this->eclgrid[i]->CalcBarycentre();
		double sum = 0;
		// check if radial model is in positiv I (EAST) direction  (KB)
		for (long i = 0; i < this->elements; i++)
		{
			for (int k = 0; k < int(this->eclgrid[i]->y_coordinates.size()); k++)
				sum = sum + this->eclgrid[i]->y_coordinates[k];
			for (int j = 0; j < int(this->eclgrid[i]->x_coordinates.size()); j++)
			{
				if ((int(this->eclgrid[i]->x_coordinates[j]) >= 0) &&
				    ( sum > -0.000000001) && ( sum < 0.00000001))
					this->RadialModellIpos = true;
				else
				{
					cout <<
					"The Radialmodell is not perpendicular to the positive I axis (East direction)!!"
					     << endl;
					cout.flush();
					exit(0);
				}
			}
		}
	}

	//Release memory
	delete (TextFile);

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                     Time: " << time << " seconds." << endl;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: DetermineNeighbourElements
   Task: Checks the consitency of both grids
   Return: true or false
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::DetermineNeighbourElements(string Projectname)
{
	string Filename;
	string tempstring;
	CReadTextfiles_ECL* TextFile;
	std::stringstream in;
	//double Element; double X; double Y; double Z; double N1; double N2; double N3; double N4; double N5; double N6;
	double X, Y, Z;
	long Element,  N[6];
	bool error = false;
	clock_t start,finish;
	double time;
	double fractpart, param, intpart;
	start = clock();

	cout << "        DetermineNeighbourElements()";
	//check if file *.neighbours exists, which stores the neighbours of each element
	Filename = Projectname + ".neighbours";
	TextFile = new CReadTextfiles_ECL;
	error = TextFile->Read_Text(Filename);
	bool Error_NeighbourFile = false;
	double epsilon = 1e-3;
	if (error == false)
	{
		//Read the neighbours of each element
		if ((TextFile->NumberOfRows - 2) == this->elements)
			for (long i = 0; i < TextFile->NumberOfRows - 2; i++)
			{
				in.str ((string) TextFile->Data[i + 1]);
				in >> Element >> X >> Y >> Z;
				for(int j = 0; j < 6; j++)
					in >> N[j];
				in.clear();

				//string Neighbours[6];
				//Neighbours[0] = N1;
				//Neighbours[1] = N2;
				//Neighbours[2] = N3;
				//Neighbours[3] = N4;
				//Neighbours[4] = N5;
				//Neighbours[5] = N6;

				//this->SplittedString.clear();
				//this->SplitStrings(TextFile->Daten[i + 1], ";");

				//if (atoi(this->SplittedString[0].data()) == i) {#
				if (Element == i)
				{
					//cout << i << ": " << atof(this->SplittedString[1].data()) << " " << this->eclgrid[i]->x_barycentre << " " << atof(this->SplittedString[2].data()) << " " << this->eclgrid[i]->y_barycentre << " " << atof(this->SplittedString[3].data()) << " " << this->eclgrid[i]->z_barycentre << endl;
					//cout << i << ": " << X << " " << this->eclgrid[i]->x_barycentre << " " << Y << " " << this->eclgrid[i]->y_barycentre << " " << Z << " " << this->eclgrid[i]->z_barycentre << endl;
					//cout << epsilon << " " << abs(X - this->eclgrid[i]->x_barycentre) << " " << abs(Y - this->eclgrid[i]->y_barycentre) << " " << abs(Z - this->eclgrid[i]->z_barycentre) << endl;
					//if ((abs(atof(this->SplittedString[1].data()) - this->eclgrid[i]->x_barycentre) < epsilon) && (abs(atof(this->SplittedString[2].data()) - this->eclgrid[i]->y_barycentre) < epsilon) && (abs(atof(this->SplittedString[3].data()) - this->eclgrid[i]->z_barycentre) < epsilon)) {
					if ((abs(X - this->eclgrid[i]->x_barycentre) < epsilon) &&
					    (abs(Y - this->eclgrid[i]->y_barycentre) < epsilon) &&
					    (abs(Z - this->eclgrid[i]->z_barycentre) < epsilon))
						//if ((atof(this->SplittedString[1].data()) == this->eclgrid[i]->x_barycentre) && (atof(this->SplittedString[2].data()) == this->eclgrid[i]->y_barycentre) && (atof(this->SplittedString[3].data()) == this->eclgrid[i]->z_barycentre)) {
						for (int j = 0; j < 6; j++)
							//this->eclgrid[i]->NeighbourElement[j] = atoi(this->SplittedString[j + 4].data());
							//for (int k = 1; k < Neighbours[6].size(); k++){
							this->eclgrid[i]->NeighbourElement[j] =
							        N[j];
					//atoi(Neighbours[j].data());
					else
						Error_NeighbourFile = true;
				}
				else
					Error_NeighbourFile = true;
			}
		else
			Error_NeighbourFile = true;

		//if (Error_NeighbourFile == true) {
		//	cout << "The program is canceled because the *.neighbours file is not correct!" << endl;
		//	system("Pause");
		//	exit(0);
		//}
		//Release memory
		delete (TextFile);
	}

	if ((error == true) || (Error_NeighbourFile == true))
	{
		cout << "           Create new *.neighbours file." << endl;
		//Set neighbourElements to -1
		for (unsigned long i = 0; i < this->eclgrid.size(); i++)
			for (unsigned int j = 0; j < this->eclgrid[i]->NeighbourElement.size(); j++)
				this->eclgrid[i]->NeighbourElement[j] = -1;

		//Determine neighboured Elements
		//Order of the faces: 0..left(x); 1..right(x); 2..front(y); 3..back(y); 4..top(z); 5..bottom(z)
		for (size_t i = 0; i < this->eclgrid.size(); i++)
			for (size_t j = i + 1; j < this->eclgrid.size(); j++)
			{
				param = static_cast<double>(j + (i * this->eclgrid.size()));
				param /= 1000;
				fractpart = modf (param, &intpart);
				if (fractpart == 0)
				{
					param = static_cast<double>(100 * (j + i * this->eclgrid.size()) /
					        (this->eclgrid.size() * this->eclgrid.size()));
					cout << "              " << Round(param, 3) << " %" << endl;
				}
				if (this->eclgrid[i]->layer == this->eclgrid[j]->layer)
				{
					if (this->eclgrid[i]->row == this->eclgrid[j]->row)
					{
						//right neighbour of the current element i
						if (this->eclgrid[i]->column ==
						    this->eclgrid[j]->column - 1)
						{
							this->eclgrid[i]->NeighbourElement[1] =
							        this->eclgrid[j]->index;
							this->eclgrid[j]->NeighbourElement[0] =
							        this->eclgrid[i]->index;
						}
						//left neighbour of the current element i
						if (this->eclgrid[i]->column ==
						    this->eclgrid[j]->column + 1)
						{
							this->eclgrid[i]->NeighbourElement[0] =
							        this->eclgrid[j]->index;
							this->eclgrid[j]->NeighbourElement[1] =
							        this->eclgrid[i]->index;
						}
					}

					if (this->eclgrid[i]->column == this->eclgrid[j]->column)
					{
						//Back neighbour of the current element i
						if (this->eclgrid[i]->row ==
						    this->eclgrid[j]->row - 1)
						{
							this->eclgrid[i]->NeighbourElement[3] =
							        this->eclgrid[j]->index;
							this->eclgrid[j]->NeighbourElement[2] =
							        this->eclgrid[i]->index;
						}
						//Front neighbour of the current element i
						if (this->eclgrid[i]->row ==
						    this->eclgrid[j]->row + 1)
						{
							this->eclgrid[i]->NeighbourElement[2] =
							        this->eclgrid[j]->index;
							this->eclgrid[j]->NeighbourElement[3] =
							        this->eclgrid[i]->index;
						}
					}
				}

				if (this->eclgrid[i]->layer == this->eclgrid[j]->layer - 1)
					//bottom neighbour of the current element i
					if ((this->eclgrid[i]->row == this->eclgrid[j]->row) &&
					    (this->eclgrid[i]->column == this->eclgrid[j]->column))
					{
						this->eclgrid[i]->NeighbourElement[5] =
						        this->eclgrid[j]->index;
						this->eclgrid[j]->NeighbourElement[4] =
						        this->eclgrid[i]->index;
					}
				if (this->eclgrid[i]->layer == this->eclgrid[j]->layer + 1)
					//top neighbour of the current element i
					if ((this->eclgrid[i]->row == this->eclgrid[j]->row) &&
					    (this->eclgrid[i]->column == this->eclgrid[j]->column))
					{
						this->eclgrid[i]->NeighbourElement[4] =
						        this->eclgrid[j]->index;
						this->eclgrid[j]->NeighbourElement[5] =
						        this->eclgrid[i]->index;
					}
			}
		//check
		//cout << this->eclgrid[i]->index << ": " << this->eclgrid[i]->NeighbourElement[0] << " " << this->eclgrid[i]->NeighbourElement[1] << " " << this->eclgrid[i]->NeighbourElement[2] << " " << this->eclgrid[i]->NeighbourElement[3] << " " << this->eclgrid[i]->NeighbourElement[4] << " " << this->eclgrid[i]->NeighbourElement[5] << endl;

		//data output
		vector <string> vec_string;
		ostringstream temp;
		tempstring = "Element  X  Y  Z  N1  N2  N3  N4  N5  N6";
		vec_string.push_back(tempstring);
		// Loop over all elements
		for (long i = 0; i < this->elements; i++)
		{
			temp.precision(12);
			temp.str("");
			temp.clear();
			temp << this->eclgrid[i]->index;
			tempstring = temp.str();
			temp.str("");
			temp.clear();
			temp << this->eclgrid[i]->x_barycentre;
			tempstring += "  " + temp.str();
			temp.str("");
			temp.clear();
			temp << this->eclgrid[i]->y_barycentre;
			tempstring += "  " + temp.str();
			temp.str("");
			temp.clear();
			temp << this->eclgrid[i]->z_barycentre;
			tempstring += "  " + temp.str();
			for (int j = 0; j < 6; j++)
			{
				temp.str("");
				temp.clear();
				temp << this->eclgrid[i]->NeighbourElement[j];
				tempstring += "  " + temp.str();
			}
			vec_string.push_back(tempstring);
		} // end element loop

		// Test Output
		Filename = Projectname + ".neighbours";
		ofstream aus;
		aus.open(Filename.data(),ios::out);
		for (unsigned int i = 0; i < vec_string.size(); i++)
			aus << vec_string[i] << endl;
		aus.close();
	}

	// Does not work if inactive cells are neglected
	//for (unsigned long i = 0; i < this->eclgrid.size(); i++) {
	//	//Calculating neighboured element in x-direction
	//	if (this->eclgrid[i]->column > 1) {
	//		this->eclgrid[i]->NeighbourElement[0] = this->eclgrid[i]->index - 1;
	//	}
	//	if (this->eclgrid[i]->column < this->columns) {
	//		this->eclgrid[i]->NeighbourElement[1] = this->eclgrid[i]->index + 1;
	//	}
	//	//Calculating neighboured element in y-direction
	//	if (this->eclgrid[i]->row > 1) {
	//		this->eclgrid[i]->NeighbourElement[2] = this->eclgrid[i]->index - this->columns;
	//	}
	//	if (this->eclgrid[i]->row < this->rows) {
	//		this->eclgrid[i]->NeighbourElement[3] = this->eclgrid[i]->index + this->columns;
	//	}
	//	//Calculating neighboured element in z-direction
	//	if (this->eclgrid[i]->layer > 1) {
	//		this->eclgrid[i]->NeighbourElement[4] = this->eclgrid[i]->index - (this->columns*this->rows);
	//	}
	//	if (this->eclgrid[i]->layer < this->layers) {
	//		this->eclgrid[i]->NeighbourElement[5] = this->eclgrid[i]->index + (this->columns*this->rows);
	//	}
	//	//check
	//	//cout << this->eclgrid[i]->index << ": " << this->eclgrid[i]->NeighbourElement[0] << " " << this->eclgrid[i]->NeighbourElement[1] << " " << this->eclgrid[i]->NeighbourElement[2] << " " << this->eclgrid[i]->NeighbourElement[3] << " " << this->eclgrid[i]->NeighbourElement[4] << " " << this->eclgrid[i]->NeighbourElement[5] << endl;
	//}
	//cout << "i" << endl;

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "          Time: " << time << " seconds." << endl;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadDataFromInputFile
   Task: Read the necessary data from .data- File: density of CO2 at surface conditions,
   necessary for calculating co2 concentration from RS value, Well Data
   Return: true or false
   Programming: 01/2011 BG
   Modification: 02/2011 KB Change from ReadSurfaceDensity to ReadDataFromInputFile
   -------------------------------------------------------------------------*/
int CECLIPSEData::ReadDataFromInputFile(std::string Filename)
{
	std::string tempstring;
	//CBoundaryConditions* m_BoundaryConditions = NULL; // unused
	CReadTextfiles_ECL* TextFile;
	std::stringstream in;
	double density;
	bool error;
	//WW bool success = false;
	std::string dummy_rate,dummy_zeile,rest_zeile, name, phase, open_flag, control_mode;
	int jj = 0;

	density = -1;
	TextFile = new CReadTextfiles_ECL;
	error = TextFile->Read_Text(Filename);

	if (error == true)
	{
		cout << "The file: " << Filename <<
		"could not been read! The program is canceled" << endl;
		system("Pause");
		exit(0);
	}

	//Read in the input file and search for the keyword DENSITY
	for (long l = 0; l < TextFile->NumberOfRows; l++)
	{
		tempstring = TextFile->Data[l].substr(0,7);
		long zeilen = 0;
		if (tempstring.compare("DENSITY") == 0)
		{
			do {
				zeilen = zeilen + 1;
			}
			while (TextFile->Data[l + zeilen] == ""); // look for possible empty rows between the Keyword and the items

			//zeilen = zeilen + 1;//End of the block is characterised with "" or "/"

			while ((TextFile->Data[l + zeilen] != "") &&
			       (TextFile->Data[l + zeilen] != "/"))
			{
				string delimiter = " ";
				SplittedString.clear();
				SplitStrings(TextFile->Data[l + zeilen], delimiter);
				if (TextFile->Data[l + zeilen].substr(0, 2) != "--") // ignore comments
				{
					density = atof(SplittedString[2].data());
					this->SurfaceCO2Density = density;
					//check if the number of bc is plausibel
					if ((density > 800) || (density < 0))
					{
						cout <<
						"The surface density of CO2 seems not to be correct: "
						     << density << endl;
						return 0;
					}
					else
						break;
				}
				zeilen = zeilen + 1;
			}
		}
		// Read Well Data
		if (this->existWells == true)
			if (tempstring.compare("WCONINJE") == 0)
			{
				//WW			success = true;
				long zeilen = 0;

				do {
					zeilen = zeilen + 1;
				}
				while (TextFile->Data[l + zeilen] != ""); // look for possible empty rows between the Keyword and the items
				jj = 0;
				for(long unsigned j = l + zeilen + 1;
				    j < (l + this->ecl_well.size() + zeilen + 1); j++)
				{
					if (this->actual_time == 0)
					{
						in.str(TextFile->Data[j]);
						in >> name;
						in >> this->ecl_well[jj]->phase;
						in >> this->ecl_well[jj]->open_flag;
						in >> this->ecl_well[jj]->control_mode;
						in >> dummy_rate;
						in >> dummy_zeile;
						if (in != 0)
							rest_zeile += dummy_zeile;

						in.clear();
					}
					jj++;
				}
			}
	}
	return 1;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadWellData
   Task: Read the well information and the time schedule for injection
   Return: true or false
   Programming: 02/2011 KB
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadWellData(std::string Filename_Wells)
{
	structWell* ecl_single_well;
	char Line[MAX_ZEILEN];
	std::string tempstring;
	std::stringstream line;
	//WW bool error = false;
	streampos position;
	std::string dollar("$");
	std::string hash("#");
	//WW bool new_subkeyword = false;
	// bool new_keyword = false; // unused

	std::string tempstring_name;
	std::string tempvalue_rate;
	double tempvalue_time;

	clock_t start,finish;
	double time;

	start = clock();

	cout << "        ReadWellSection()";

	ifstream in(Filename_Wells.data(),ios::in);

	if (!in)
		std::cout << "File not found." << std::endl;
	//WW	error = true;
	else
		while (!in.eof())
		{
			in.getline(Line, MAX_ZEILEN);
			position = in.tellg();
			tempstring = Line;
			if (tempstring.find("#STOP") != std::string::npos)
				std::cout << "ok" << std::endl;

			//if(tempstring.size() < 1) break;
			//if(tempstring.find(hash)!=string::npos) {
			//	new_keyword = true;
			//	break;
			//}

			//NAME
			if(tempstring.find("$NAME") != string::npos) //subkeyword found
			{
				ecl_single_well = new structWell;
				in >> tempstring_name;
				ecl_single_well->name = tempstring_name;
				in.clear();

				//continue;
			}

			//TIMECURVE
			if(tempstring.find("$TIMECURVE") != string::npos) //subkeyword found
			{
				while (tempstring.find("#") != 0)
				{
					position = in.tellg();
					in.getline(Line, MAX_ZEILEN);
					tempstring = Line;
					if (tempstring.find("#") != string::npos)
						break;
					if(tempstring.find("$") != string::npos)
						//WW					new_subkeyword = true;
						break;
					in.seekg(position);
					in >> tempvalue_time;
					ecl_single_well->time.push_back(tempvalue_time);
					in >> tempvalue_rate;
					ecl_single_well->rate.push_back(tempvalue_rate);
					//if ((double(ecl_single_well->time.begin())) != "0"){
					//	cout << "Warning: first item of time in timecurve is not 0!" << endl;
					//	break;
					//}
					in.clear();
					in.getline(Line, MAX_ZEILEN);
					tempstring = Line;

					//position = in.tellg();
					//in.getline(Line, MAX_ZEILEN);
					//tempstring = Line;
					//string delimiter=" ";
					//SplittedString.clear();
					//
					//SplitStrings(tempstring, delimiter);

					//for (unsigned int m = 0; m < SplittedString.size(); m++) {
					//		if (m == 0)ecl_single_well->time = atof(SplittedString[m]);
					//		if (m == 1)ecl_single_well->time.push_back(SplittedString[m]);
					//}

					//position = in.tellg();
					//in.getline(Line, MAX_ZEILEN);
					//tempstring = Line;
					//in >> tempvalue_time;
					//ecl_single_well->time.push_back(tempvalue_time);
					//in >> tempvalue_rate;
					//ecl_single_well->rate.push_back(tempvalue_rate);
					//in.clear();
					//in.seekg(position);
					//in.getline(Line, MAX_ZEILEN);

					//continue;
				}
				if (!tempstring_name.empty())
					// Kontrolle, dass Name nicht leer!
					this->ecl_well.push_back(ecl_single_well);
			}
		}
	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "          Time: " << time << " seconds." << endl;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReadPositionBoundaryCondition
   Task: Read the position of boundary conditions
   Return: true or false
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReadPositionBoundaryCondition(std::string Filename)
{
	std::string tempstring;
	CBoundaryConditions* m_BoundaryConditions = NULL;
	CReadTextfiles_ECL* TextFile;
	std::stringstream in;
	bool error = false;
	int zeile;
	int number;
	int i1;
	int i2;
	int j1;
	int j2;
	int k1;
	int k2;
	int index;
	std::string boundary_position;

	TextFile = new CReadTextfiles_ECL;
	error = TextFile->Read_Text(Filename);

	if (error == true)
	{
		cout << "The file: " << Filename <<
		"could not been read! The program is canceled" << endl;
		system("Pause");
		exit(0);
	}

	//Read in the input file and search for the keyword AQUANCON
	for (long l = 0; l < TextFile->NumberOfRows; l++)
	{
		tempstring = TextFile->Data[l].substr(0,8);
		zeile = 1;
		if (tempstring.compare("AQUANCON") == 0)
		{
			//End of the block is characterised with "" or "/"

			while ((TextFile->Data[l + zeile] != "") &&
			       (TextFile->Data[l + zeile] != "/"))
			{
				std::string delimiter = " ";
				SplittedString.clear();
				SplitStrings(TextFile->Data[l + zeile], delimiter);

				if (TextFile->Data[l + zeile].substr(0, 2) != "--") // ignore comments
				{
					number = atoi(SplittedString[0].data());
					//check if the number of bc is plausibel
					if (number > 200)
					{
						cout <<
						"There were more than 200 bc found in the Eclipse model. Please make sure that the definitions at the keyword AQUANCON are correct!"
						     << endl;
						return false;
					}
					i1 = atoi(SplittedString[1].data());
					i2 = atoi(SplittedString[2].data());
					j1 = atoi(SplittedString[3].data());
					j2 = atoi(SplittedString[4].data());
					k1 = atoi(SplittedString[5].data());
					k2 = atoi(SplittedString[6].data());
					boundary_position = SplittedString[7].substr(1,2);
					index = -1;
					//order of cells in output files: z1, y1, x1..xn; z1, y2, x1..xn; z2, y1, x1..xn
					for (int k = k1; k <= k2; k++)
					{
						for (int j = j1; j <= j2; j++)
							for (int i = i1; i <= i2; i++)
							{
								index = index + 1;
								m_BoundaryConditions =
								        new CBoundaryConditions(); // new member of eclblock
								m_BoundaryConditions->index = index;
								m_BoundaryConditions->number =
								        number;
								m_BoundaryConditions->
								boundary_position =
								        boundary_position;
								for (long m = 0; m < this->elements;
								     m++)
									if ((this->eclgrid[m]->
									     column == i) &&
									    (this->eclgrid[m]->row
									     ==
									     j) &&
									    (this->eclgrid[m]->
									     layer ==
									     k))
										m_BoundaryConditions
										->connected_element
										        = this->
										          eclgrid[m
										          ]->
										          index;

								this->BC.push_back(
								        m_BoundaryConditions);
								long index_boundary =
								        long(this->BC.size() - 1);
								this->eclgrid[m_BoundaryConditions
								              ->connected_element]
								->
								ConnectedBoundaryCondition.
								push_back(
								        index_boundary);
							}
					}
					if (m_BoundaryConditions->connected_element == -1)
						return false;
				}
				zeile = zeile + 1;
			}
			l = TextFile->NumberOfRows;
		}
	}
	return true;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadBoundaryData
   Task: Read the flow data for boundary conditions
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReadBoundaryData(int index_boundary, vector <string> Data)
{
	int counter;
	int number_values = 4;
	int number_cell;
	int number_value;

	//4 numbers are given for each boundary cell (order of the cells at the moment unknown!!!)
	counter = -1;
	for (unsigned long i = 0; i < Data.size(); i++)
	{
		this->SplittedString.clear();
		if (i < (Data.size() - 1))
		{
			std::string SplitBoundaryData[3];

			//this->SplitBoundaryData.push_back(Data[i].substr(3,20));
			//this->SplitBoundaryData.push_back(Data[i].substr(26,20));
			//this->SplitBoundaryData.push_back(Data[i].substr(49,20));

			SplitBoundaryData[0] = Data[i].substr(3,20);
			SplitBoundaryData[1] = Data[i].substr(26,20);
			SplitBoundaryData[2] = Data[i].substr(49,20);

			//cout << SplitBoundaryData.size() << endl;

			//for (unsigned int j = 0; j < this->SplittedString.size(); j++){
			for (unsigned int j = 0; j < 3; j++)
			{
				counter = counter + 1;
				number_cell = int(counter / number_values);
				number_value = counter - (number_cell * number_values);
				//search for the corresponding BC (the same boundary number and the index indentical to the counter
				for (unsigned int k = 0; k < this->BC.size(); k++)
					if ((this->BC[k]->number == index_boundary) &&
					    (this->BC[k]->index == number_cell))
					{
						//this->BC[k]->value[number_value] = atof(this->SplittedString[j].data());
						this->BC[k]->value[number_value] = atof(
						        SplitBoundaryData[j].data());
						k = int(this->BC.size());
					}
			}
		}
		else
		{
			this->SplittedString.clear();
			this->SplitStrings(Data[i]," ");
			if (i == 0)
				if (SplittedString.size() != 3)
					return false;
			//string SplitBoundaryData[2];

			//SplitBoundaryData[0] = Data[i].substr(3,20);
			//SplitBoundaryData[1] = Data[i].substr(26,20);

			//for (unsigned int j = 0; j < 2; j++){
			for (unsigned int j = 0; j < this->SplittedString.size(); j++)
			{
				counter = counter + 1;
				number_cell = int(counter / number_values);
				number_value = counter - (number_cell * number_values);
				//search for the corresponding BC (the same boundary number and the index indentical to the counter
				for (unsigned int k = 0; k < this->BC.size(); k++)
					if ((this->BC[k]->number == index_boundary) &&
					    (this->BC[k]->index == number_cell))
					{
						this->BC[k]->value[number_value] = atof(
						        this->SplittedString[j].data());
						//this->BC[k]->value[number_value] = atof(SplitBoundaryData[j].data());
						k = int(this->BC.size());
					}
			}
		}
	}

	return true;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ReadEclipseData
   Task: Read the data of the current time step of the Eclipse model from output file
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::ReadEclipseData(std::string Filename, long Timestep)
{
	(void)Timestep; // unused
	vector<string> files = vector<string>();
	CReadTextfiles_ECL* TextFile;
	std::string tempstring;
	bool saturation_water, saturation_gas, saturation_oil;
	clock_t start,finish;
	double time;

	start = clock();

	saturation_water = saturation_gas = saturation_oil = false;

	cout << "        ReadEclipseData() ";
	//get memory for data structures
	if(Data == NULL) // allocate first time
	{
		Data = (double**) malloc((this->elements) * sizeof(double*));
		for(long i = 0; i < (this->elements); i++)
			Data[i] = (double*) malloc(this->numberOutputParameters * sizeof(double));
		//for(long i = 0; i < (this->elements); i++)
		//	for(long j = 0; j < this->times; j++)
		//		Data[i][j] = (double *) malloc(this->numberOutputParameters*sizeof(double));

		//Deallokieren
		/*Data2 = (double **) malloc(1000000 * sizeof(double*));
		   for(long i = 0; i < (1000000); i++)
		   Data2[i] = (double *) malloc(100 * sizeof(double));

		   for(long i=0;i<1000000;i++) free(Data2[i]);

		   free(Data2);*/
	}

	//initialise data structure
	for(long i = 0; i < (this->elements); i++)
		for(int k = 0; k < this->numberOutputParameters; k++)
			Data[i][k] = 0;

	//Reads the text file (usually 60 bytes in one row of the eclipse output
	bool Error;
	TextFile = new CReadTextfiles_ECL;
	Error = TextFile->Read_Text(Filename);
	if (Error == true)
	{
		cout << "The program is canceled" << endl;
		system("Pause");
		exit(0);
	}
	for (long j = 0; j < TextFile->NumberOfRows; j++)
	{
		//if (j==565 || j ==1384) {
		//	cout << TextFile->Daten[j] << endl;
		//}
		double Multiplier = 1;
		if (TextFile->Data[j].length() > 0 && TextFile->Data[j].substr(1,1) == "'")
		{
			std::string tempstring = TextFile->Data[j].substr(2,8);
			tempstring = tempstring.substr(0,tempstring.find_first_of(" "));
			for (unsigned int k = 0; k < this->Variables.size(); k++)
				if (tempstring.compare(this->Variables[k]) == 0)
				{
					if (this->Variables[k].compare("ACAQNUM") == 0)
					{
						int index_boundary = atoi(
						        TextFile->Data[j + 1].substr(0,12).data());
						this->SplittedString.clear();
						this->SplitStrings(TextFile->Data[j + 2], " ");
						int temprows =
						        int(atoi(this->SplittedString[2].data()) /
						            3) + 1;
						vector <string> temp_Daten;
						for (int l = 0; l < temprows; l++)
							temp_Daten.push_back(TextFile->Data[j + 3 +
							                                    l]);
						if (this->ReadBoundaryData(index_boundary,
						                           temp_Daten) == false)
						{
							cout <<
							"Error while reading boundary data." <<
							endl;
							system("Pause");
							exit(0);
						}
					}
					else
					{
						//convert pressure units from bar to Pa
						if ((this->Variables[k].compare("PRESSURE") ==
						     0) ||
						    (this->Variables[k].compare("PCOW") == 0) ||
						    (this->Variables[k].compare("PCOG") == 0) ||
						    (this->Variables[k].compare("PWAT") == 0) ||
						    (this->Variables[k].compare("PGAS") == 0) ||
						    (this->Variables[k].compare("POIL") == 0))
							Multiplier = 100000.0;
						// consider different orientation of the z-axis in Eclipse and GeoSys
						if ((this->Variables[k].compare("FLOWATK+") ==
						     0) ||
						    (this->Variables[k].compare("FLOOILK+") ==
						     0) ||
						    (this->Variables[k].compare("FLOGASK+") == 0))
							Multiplier = -1.0;

						if (this->Variables[k].compare("RS") == 0)
							Multiplier = 1.0;  // Unit: m� gas per m� oil

						if (this->Variables[k].compare("GAS_DEN") == 0)
							Multiplier = 1.0;  // Unit: kg/m�

						//set a flag of which phase saturation was red
						if (this->Variables[k].compare("SWAT") == 0)
							saturation_water = true;
						if (this->Variables[k].compare("SGAS") == 0)
							saturation_gas = true;
						if (this->Variables[k].compare("SOIL") == 0)
							saturation_oil = true;

						// read number of datapoints
						double tempNumber =
						        atoi(TextFile->Data[j].substr(12,13).data());
						long temprows = 0;
						if (tempNumber == this->elements)
						{
							temprows = int(ceil(tempNumber / 4.0));
							//Read data for identified variable
							long rowindex = 0;
							for (long l = j + 1; l < j + temprows + 1;
							     l++)
							{
								if (l < (j + temprows))
								{
									std::string Dataline[4]; // Definition einen strings mit 4 Arrays, Arrays immer in eckigen Klammern definieren
									Dataline[0] =
									        TextFile->Data[l].
									        substr(2,
									               15);
									Dataline[1] =
									        TextFile->Data[l].
									        substr(19,
									               15);
									Dataline[2] =
									        TextFile->Data[l].
									        substr(36,
									               15);
									Dataline[3] =
									        TextFile->Data[l].
									        substr(53,
									               15);
									for (unsigned int m = 0;
									     m < 4; m++)
									{
										rowindex =
										        rowindex +
										        1; //z�hlt die gesplitteten Zeilen in der Zeile, in Variante KB4 entspricht das 4
										this->Data[rowindex
										           - 1][k]
										        = atof(
										        Dataline[m]
										        .data())
										          *
										          Multiplier; //schreibt die Zeile aus Dataline mit allen m und dem Keyword k in das double Data rein
										//cout << "Element: " << rowindex << " Variable: " << this->Variables[k] << " Value: " << this->Data[rowindex-1][k] << endl;
									}
								}
								else
								{
									std::string delimiter = " ";
									SplittedString.clear();
									SplitStrings(
									        TextFile->Data[l],
									        delimiter);
									//cout << TextFile->Daten[l] << endl;
									for (unsigned int m = 0;
									     m <
									     SplittedString.size();
									     m++)
									{
										rowindex =
										        rowindex +
										        1;
										//cout << rowindex-1 << " " << k << " " << " " << this->eclgrid[rowindex-1]->x_barycentre <<  " " << this->eclgrid[rowindex-1]->y_barycentre  << " " << this->eclgrid[rowindex-1]->z_barycentre  << " " << atof(SplittedString[m].data()) << endl;
										this->Data[rowindex
										           - 1][k]
										        = atof(
										        SplittedString
										        [m].
										        data()) *
										          Multiplier;
										//cout << "Element: " << rowindex << " Variable: " << this->Variables[k] << " Value: " << this->Data[rowindex-1][k] << endl;
									}
								}
							}
							j = j + temprows;
							//a[k][1] = j + 1;
						}
						else
						{
							// there are inactive cells
							//cout << "Error while reading the data file. The number of data points doesn't fit to the grid." << endl;
							temprows = int(ceil(tempNumber / 4.0));
							//Read data for identified variable
							long rowindex = 0;
							for (long l = j + 1; l < j + temprows + 1;
							     l++)
							{
								if ( l < (j + temprows))
								{
									std::string Dataline[4]; // Definition einen strings mit 4 Arrays, Arrays immer in eckigen Klammern definieren
									Dataline[0] =
									        TextFile->Data[l].
									        substr(3,
									               15);
									Dataline[1] =
									        TextFile->Data[l].
									        substr(20,
									               15);
									Dataline[2] =
									        TextFile->Data[l].
									        substr(38,
									               15);
									Dataline[3] =
									        TextFile->Data[l].
									        substr(55,
									               15);
									//string delimiter=" ";
									//SplittedString.clear();
									//SplitStrings(TextFile->Daten[l], delimiter);
									//cout << TextFile->Daten[l] << endl;
									for (unsigned int m = 0;
									     m < 4; m++)
									{
										rowindex =
										        rowindex +
										        1;
										while (this->
										       eclgrid[
										               rowindex
										               -
										               1]->
										       active ==
										       0 &&
										       rowindex <
										       this->
										       elements)
											rowindex =
											        rowindex
											        + 1;
										if (rowindex <
										    this->elements)
											this->Data[
											        rowindex
											        -
											        1][
											        k]
											        =
											                atof(
											                        Dataline
											                        [
											                                m
											                        ]
											                        .
											                        data())
											                *
											                Multiplier;
										else
											cout <<
											"1 data point couldn't be allocated to a grid point"
											     <<
											endl;
									}
								}
								else
								{
									std::string delimiter = " ";
									SplittedString.clear();
									SplitStrings(
									        TextFile->Data[l],
									        delimiter);
									//cout << TextFile->Daten[l] << endl;
									for (unsigned int m = 0;
									     m <
									     SplittedString.size();
									     m++)
									{
										rowindex =
										        rowindex +
										        1;
										//cout << rowindex-1 << " " << timestep << " " << k << " " << " " << this->eclgrid[rowindex-1]->x_barycentre <<  " " << this->eclgrid[rowindex-1]->y_barycentre  << " " << this->eclgrid[rowindex-1]->z_barycentre  << " " << atof(SplittedString[m].data()) << endl;
										this->Data[rowindex
										           - 1][k]
										        = atof(
										        SplittedString
										        [m].
										        data()) *
										          Multiplier;
										//	cout << "Element: " << rowindex << " Variable: " << this->Variables[k] << " Value: " << this->Data[rowindex-1][timestep][k] << endl;
									}
								}
							}
							j = j + temprows;
						}
					}
				}
		}
	}

	// Release Textfile object - sb
	//delete(TextFile->Daten);
	delete (TextFile);

	if (this->E100 == true)
	{
		//--------------------------------------------------------------------------------------------
		// calculating phase pressure in the case of Eclipse E100 from pressure and capillary pressure
		if (int(this->Phases.size()) == 1)
		{
			int index_pwat, index_pgas, index_poil, index_pressure;
			int index_aim = -1;
			index_pwat = this->GetVariableIndex("PWAT");
			index_pgas = this->GetVariableIndex("PGAS");
			index_poil = this->GetVariableIndex("POIL");
			index_pressure = this->GetVariableIndex("PRESSURE");

			if (index_pwat >= 0)
				index_aim = index_pwat;
			if (index_pgas >= 0)
				index_aim = index_pgas;
			if (index_poil >= 0)
				index_aim = index_poil;

			// Transfer the pressure
			for (int i = 0; i < this->elements; i++)
				this->Data[i][index_aim] = this->Data[i][index_pressure];
		}
		//assumption: if water and oil -> pressure = oil pressure
		if (int(this->Phases.size()) == 2)
		{
			// water and oil
			if ((this->Phases[0] == "WATER") && (this->Phases[1] == "OIL"))
			{
				// Get the variable index
				int index_pwat, index_poil, index_pressure, index_pcap;
				index_pwat = this->GetVariableIndex("PWAT");
				index_poil = this->GetVariableIndex("POIL");
				index_pressure = this->GetVariableIndex("PRESSURE");
				index_pcap = this->GetVariableIndex("PCOW");

				// Transfer the pressure
				for (int i = 0; i < this->elements; i++)
				{
					this->Data[i][index_pwat] = this->Data[i][index_pressure];
					this->Data[i][index_poil] = this->Data[i][index_pwat] +
					                            this->Data[i][index_pcap];
				}
			}
			if ((this->Phases[0] == "WATER") && (this->Phases[1] == "GAS"))
			{
				cout << endl;
				cout <<
				"GAS-WATER System can not be considered with ECLIPSE E100 and GeoSys"
				     << endl;
				system("Pause");
				exit(0);
			}
			if ((this->Phases[0] == "OIL") && (this->Phases[1] == "GAS"))
			{
				// Get the variable index
				int index_pgas, index_poil, index_pressure, index_pcap;
				index_pgas = this->GetVariableIndex("PGAS");
				index_poil = this->GetVariableIndex("POIL");
				index_pressure = this->GetVariableIndex("PRESSURE");
				index_pcap = this->GetVariableIndex("PCOG");

				// Transfer the pressure
				for (int i = 0; i < this->elements; i++)
				{
					this->Data[i][index_poil] = this->Data[i][index_pressure];
					this->Data[i][index_pgas] = this->Data[i][index_poil] +
					                            this->Data[i][index_pcap];
				}
			}
		}
		if (int(this->Phases.size()) == 3)
		{
			cout << endl;
			cout <<
			"Currently not more than 2 phases are considered for reading eclipse pressure"
			     << endl;
			cout <<
			"Three phases are only valid if a water saturation exists only at the boundaries of the model domain!!";
			//system("Pause");
			//exit(0);

			int index_pwat, index_poil, index_pgas, index_pressure,
			    index_pcap_oil_water, index_pcap_oil_gas;
			index_pwat = this->GetVariableIndex("PWAT");
			index_poil = this->GetVariableIndex("POIL");
			index_pgas = this->GetVariableIndex("PGAS");
			index_pressure = this->GetVariableIndex("PRESSURE");
			index_pcap_oil_water = this->GetVariableIndex("PCOW");
			index_pcap_oil_gas = this->GetVariableIndex("PCOG");

			// Transfer the pressure
			for (int i = 0; i < this->elements; i++)
			{
				this->Data[i][index_pwat] = this->Data[i][index_pressure];
				this->Data[i][index_poil] = this->Data[i][index_pwat] +
				                            this->Data[i][index_pcap_oil_water];
				this->Data[i][index_pgas] = this->Data[i][index_poil] +
				                            this->Data[i][index_pcap_oil_gas];
			}
		}

		//--------------------------------------------------------------------------------------------
		// calculating phase saturation for all phases

		// 2 existing phases
		if (this->Phases.size() == 2)
		{
			int index_aim = -1, index_source = -1;
			int index_swat, index_sgas, index_soil;
			index_swat = this->GetVariableIndex("SWAT");
			index_sgas = this->GetVariableIndex("SGAS");
			index_soil = this->GetVariableIndex("SOIL");

			if ((index_swat >= 0) && (saturation_water == true))
			{
				index_source = index_swat;
				if (index_sgas >= 0)
					index_aim = index_sgas;
				if (index_soil >= 0)
					index_aim = index_soil;
			}
			if ((index_sgas >= 0) && (saturation_gas == true))
			{
				index_source = index_sgas;
				if (index_swat >= 0)
					index_aim = index_swat;
				if (index_soil >= 0)
					index_aim = index_soil;
			}
			if ((index_soil >= 0) && (saturation_oil == true))
			{
				index_source = index_soil;
				if (index_swat >= 0)
					index_aim = index_swat;
				if (index_sgas >= 0)
					index_aim = index_sgas;
			}

			// Calculating the phase saturation
			for (int i = 0; i < this->elements; i++)
				this->Data[i][index_aim] = 1 - this->Data[i][index_source];
		}
		// 3 existing phases
		if (this->Phases.size() == 3)
		{
			int index_aim = -1, index_source1 = -1, index_source2 = -1;
			int index_swat, index_sgas, index_soil;
			index_swat = this->GetVariableIndex("SWAT");
			index_sgas = this->GetVariableIndex("SGAS");
			index_soil = this->GetVariableIndex("SOIL");

			if (saturation_water == false)
			{
				index_aim = index_swat;
				index_source1 = index_sgas;
				index_source2 = index_soil;
			}
			if (saturation_gas == false)
			{
				index_aim = index_sgas;
				index_source1 = index_swat;
				index_source2 = index_soil;
			}
			if (saturation_oil == false)
			{
				index_aim = index_soil;
				index_source1 = index_sgas;
				index_source2 = index_swat;
			}
			// Calculating the phase saturation
			for (int i = 0; i < this->elements; i++)
				this->Data[i][index_aim] = 1 - this->Data[i][index_source1] -
				                           this->Data[i][index_source2];
			//cout << "Element: " << i << " RS: " << this->Data[i][this->GetVariableIndex("RS")] << " SWAT: " << this->Data[i][index_swat];
			//cout  << " SOIL: " << this->Data[i][index_soil] << " SGAS: " << this->Data[i][index_sgas] << endl;
		}

		if (this->Phases.size() > 3)
		{
			cout <<
			"Currently not more than 3 phases are considered for reading eclipse data"
			     << endl;
			system("Pause");
			exit(0);
		}
	}
	//--------------------------------------------------------------------------------------------
	// output of vertical flow data
	//for (unsigned long i = 0; i < this->elements; i++) {
	//	cout << "Q(k) " << this->Data[i][0][2] << endl;
	//}
	//cout << "Fertig";

	////Test output
	//vector <string> vec_string;
	//ostringstream temp;
	//tempstring = "Element; X; Y; Z";
	//for (unsigned int k = 0; k < this->Variables.size(); k++) {
	//	tempstring = tempstring + "; " + this->Variables[k];
	//}
	//vec_string.push_back(tempstring);
	//// Loop over all elements
	//for (long i = 0; i < this->elements; i++){
	//	CECLIPSEBlock *elem = this->eclgrid[i];
	//	temp.str(""); temp.clear(); temp << i; tempstring = temp.str();
	//	temp.str(""); temp.clear(); temp << elem->x_barycentre; tempstring += "; " + temp.str();
	//	temp.str(""); temp.clear(); temp << elem->y_barycentre; tempstring += "; " + temp.str();
	//	temp.str(""); temp.clear(); temp << elem->z_barycentre; tempstring += "; " + temp.str();

	//	for (unsigned int k = 0; k < this->Variables.size(); k++) {
	//		temp.str(""); temp.clear(); temp.precision(12); temp << this->Data[i][k]; tempstring += "; " + temp.str();
	//	}
	//	//if (i == 470) {
	//	//	cout << i;
	//	//}
	//	vec_string.push_back(tempstring);
	//}  // end element loop

	////Test Output
	//int position = Filename.find_last_of("\\");
	//string path = Filename.substr(0,position);
	//position = path.find_last_of("\\");
	//path = path.substr(0,position);
	////temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
	//string aus_file = path + "\\CheckDataRedIn_0.csv";
	//ofstream aus;
	//aus.open(aus_file.data(),ios::out);
	//for (unsigned int i = 0; i < vec_string.size(); i++) {
	//	aus << vec_string[i] << endl;
	//}
	//aus.close();

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                    Time: " << time << " seconds." << endl;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CalculateRSfromMassFraction_E300
   Task: Uses component flux (E300) to calculate phase flux and calculate the dissolved CO2 concentration in water
   Return: nothing
   Programming: 01/2011 BG
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::CalculateRSfromMassFraction_E300()
{
	//CElem* m_ele = NULL;
	//CFEMesh* m_msh = fem_msh_vector[0];
	clock_t start,finish;
	double time;
	double xw_CO2_liquid;
	double RS;
	//WW int index_xw_CO2_liquid;
	int index_density_liquid;
	double rho_liquid;

	start = clock();

	cout << "        CalculatePhaseFlux() ";

	//WW index_xw_CO2_liquid = this->GetVariableIndex("XMF2");
	index_density_liquid = this->GetVariableIndex("DENW");

	if ((index_density_liquid < 0))
	{
		cout <<
		"Not all variables are existing in the eclipse output files that are necessary to write data back to e300!"
		     << endl;
		system("Pause");
		exit(0);
	}

	// Calculate the phase flux
	for (long i = 0; i < this->elements; i++)
	{
		// mole fraction of component c at phase alpha
		xw_CO2_liquid = this->Data[i][this->GetVariableIndex("XFW2")];
		// porosity
		rho_liquid = this->Data[i][this->GetVariableIndex("DENW")];
		// calculate pseudo RS value from the mass fraction of each component on the liquid phase
		//RS [m�_CO2 / m�_liquid] = (V_liquid * rho_liquid * xw_CO2_liquid) / (rho_CO2_surface * V_liquid) = rho_liquid [kg_liq / m�_liq] * xw_CO2_liquid [kg_CO2 / kg_liq] / rho_CO2_surface [kg_CO2 / m�_CO2]
		RS = rho_liquid * xw_CO2_liquid / this->SurfaceCO2Density;

		// store RS in E100 variables
		this->Data[i][this->GetVariableIndex("RS")] = RS;
	}
	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                    Time: " << time << " seconds." << endl;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReplaceASectionInFile
   Task: Determine corresponding nodes and elements between Geosys and Eclipse
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReplaceASectionInFile(std::string Filename,
                                         string Keyword,
                                         vector <std::string> Data,
                                         bool CheckLengthOfSection)
{
	CReadTextfiles_ECL* TextFile;
	std::string tempstring;
	bool success = false;

	//Reads the text file (usually 60 bytes in one row of the eclipse output
	bool Error;
	TextFile = new CReadTextfiles_ECL;
	Error = TextFile->Read_Text(Filename);
	if (Error == true)
	{
		cout << "The program is canceled" << endl;
		system("Pause");
		exit(0);
	}

	//scan file for the keyword and replace the following data
	for (long i = 0; i < TextFile->NumberOfRows; i++)
		if (TextFile->Data[i].substr(0, Keyword.length()) == Keyword)
		{
			//count data rows in the file and compare it to the rows in the new Data array, if equal than replace data
			//count data rows in the file and compare it to the rows in the new Data array, if equal than replace data
			success = true;
			long zeilen = 0;
			if (CheckLengthOfSection == true)
			{
				do {
					zeilen = zeilen + 1;
				}
				while (TextFile->Data[i + zeilen] != "");
				zeilen = zeilen - 1;
			}
			else
				zeilen = long(Data.size());
			if (zeilen == long(Data.size()))
				//Replace Data in TextFile variable
				for (long j = 0; j < zeilen; j++)
					TextFile->Data[i + 1 + j] = Data[j];

			else
			{
				cout <<
				"Replacing a Section in the Eclipse input file is not possible because the section length doesn't fit!"
				     << endl;
				success = false;
				return success;
			}
		}
	//Rewrite the file
	CWriteTextfiles_ECL* OutputFile;
	OutputFile = new CWriteTextfiles_ECL;
	OutputFile->Write_Text(Filename, TextFile->Data);

	//Release memory
	delete (TextFile);

	return success;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ReplaceWellRate
   Task: Replace Injection Rate with Data from .well-File
   Return: nothing
   Programming: 02/2011 KB
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ReplaceWellRate(std::string Filename, std::string Keyword_well)
{
	CReadTextfiles_ECL* TextFile;
	std::string tempstring;
	bool success = false;
	stringstream in;
	std::string dummy_rate,dummy_zeile,rest_zeile, name, phase, open_flag, control_mode;
	int jj = 0;

	//Reads the text file (usually 60 bytes in one row of the eclipse output
	bool Error;
	TextFile = new CReadTextfiles_ECL;
	Error = TextFile->Read_Text(Filename);
	if (Error == true)
	{
		cout << "The program is canceled" << endl;
		system("Pause");
		exit(0);
	}

	//scan file for the keyword and replace the following data
	for (long i = 0; i < TextFile->NumberOfRows; i++)
		if (TextFile->Data[i].substr(0, Keyword_well.length()) == Keyword_well)
		{
			//count data rows in the file and compare it to the rows in the new Data array, if equal than replace data
			success = true;
			long zeilen = 0;

			do {
				zeilen = zeilen + 1;
			}
			while (TextFile->Data[i + zeilen] != ""); // look for possible empty rows between the Keyword and the items
			jj = 0;
			for(long unsigned j = i + zeilen + 1;
			    j < (i + this->ecl_well.size() + zeilen + 1); j++)
			{
				if (this->actual_time == 0)
				{
					in.str(TextFile->Data[j]);

					in >> name;
					in >> this->ecl_well[jj]->phase;
					in >> this->ecl_well[jj]->open_flag;
					in >> this->ecl_well[jj]->control_mode;
					in >> dummy_rate;
					in >> dummy_zeile;
					if (in != 0)
						rest_zeile += dummy_zeile;

					in.clear();

					std::string outline = this->ecl_well[jj]->name + "   " +
					                      this->ecl_well[jj]->phase + "   " +
					                      this->ecl_well[jj]->open_flag +
					                      "   " +
					                      this->ecl_well[jj]->control_mode +
					                      "   " +
					                      this->WellRates[jj] + "   " + " / ";
					TextFile->Data[j] = outline;
					jj++;
				}
				else
				{
					std::string outline = this->ecl_well[jj]->name + "   " +
					                      this->ecl_well[jj]->phase + "   " +
					                      this->ecl_well[jj]->open_flag +
					                      "   " +
					                      this->ecl_well[jj]->control_mode +
					                      "   " +
					                      this->WellRates[jj] + "   " + " / ";
					TextFile->Data[j] = outline;
					jj++;
				}
			}
		}

	//Rewrite the file
	CWriteTextfiles_ECL* OutputFile;
	OutputFile = new CWriteTextfiles_ECL;
	OutputFile->Write_Text(Filename, TextFile->Data);

	//Release memory
	delete (TextFile);

	return success;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CorrespondingElements
   Task: Determine corresponding nodes and elements between Geosys and Eclipse
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CorrespondingElements()
{
	MeshLib::CElem* m_ele = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];

	//check if number of elements is equal
	if (long(m_msh->ele_vector.size()) != this->elements)
	{
		cout << "Error: The number of elements between Geosys and Eclipse is not equal!" <<
		endl;
		return false;
	}

	//set element vectors to the start value -1
	CorrespondingEclipseElement.resize(m_msh->ele_vector.size());
	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
		CorrespondingEclipseElement[i] = -1;
	CorrespondingGeosysElement.resize(m_msh->ele_vector.size());
	for (long j = 0; j < this->elements; j++)
		CorrespondingGeosysElement[j] = -1;

	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
		for (long j = 0; j < this->elements; j++)
		{
			m_ele = m_msh->ele_vector[i];
			double const* grav_c(m_ele->GetGravityCenter());
			//check if coordinates of the gravity centre are equal
			if ((grav_c[0] == this->eclgrid[j]->x_barycentre) &&
			    (grav_c[1] == this->eclgrid[j]->y_barycentre) &&
			    (grav_c[2] == -this->eclgrid[j]->z_barycentre))
			{
				CorrespondingEclipseElement[i] = j;
				CorrespondingGeosysElement[j] = i;
			}
		}

	//check if all values in the correspondingElementVector are larger than -1
	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
	{
		cout << "   CorrespondingEclipseElement[" << i << "] " <<
		CorrespondingEclipseElement[i] << endl;
		if (CorrespondingEclipseElement[i] < 0)
		{
			cout << "Error: No Eclipse element linked to the Geosys element!" << endl;
			return 0;
		}
	}
	for (long j = 0; j < this->elements; j++)
		if (CorrespondingGeosysElement[j] < 0)
		{
			cout << "Error: No Geosys element linked to the Eclipse element!" << endl;
			return 0;
		}
	return 1;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: CompareElementsGeosysEclipse
   Task: Check whether the Eclipse and the Geosys grid is identical and has the same order
   Return: bool value if there occured an error
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CompareElementsGeosysEclipse()
{
	MeshLib::CElem* m_ele = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];
	Math_Group::vec<MeshLib::CNode*> ele_nodes(8);
	clock_t start,finish;
	double time;
	double epsilon = 1e-7;

	start = clock();

	cout << "        CompareElements() ";

	//check if number of elements is equal
	if (long(m_msh->ele_vector.size()) != this->elements)
	{
		cout << "Error: The number of elements between Geosys and Eclipse is not equal!" <<
		endl;
		return false;
	}

	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
	{
		m_ele = m_msh->ele_vector[i];
		m_ele->GetNodes(ele_nodes);
		bool ElementIsEqual = true;
		// Test output
		//for(int j=0;j<8;j++){
		//	a_node = ele_nodes[j];
		//	cout << endl;
		//	cout << a_node->X() << " " << a_node->Y() << " " << a_node->Z() << " ECL " << this->eclgrid[i]->x_coordinates[j] << " " << this->eclgrid[i]->y_coordinates[j]  << " " << this->eclgrid[i]->z_coordinates[j] << endl;
		//}
		//A loop can't be used because the order of nodes is not equal

		double const* pnt (ele_nodes[0]->getData());
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[0]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[0]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[0]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[1]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[1]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[1]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[1]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[2]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[3]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[3]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[3]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[3]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[2]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[2]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[2]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[4]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[4]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[4]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[4]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[5]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[5]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[5]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[5]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[6]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[7]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[7]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[7]) > epsilon)
			ElementIsEqual = false;

		pnt = ele_nodes[7]->getData();
		if (fabs(pnt[0] - this->eclgrid[i]->x_coordinates[6]) > epsilon ||
		    fabs(pnt[1] - this->eclgrid[i]->y_coordinates[6]) > epsilon  ||
		    fabs(pnt[2] - -this->eclgrid[i]->z_coordinates[6]) > epsilon)
			ElementIsEqual = false;

		if (ElementIsEqual == false)
			return 0;

		//Get the gravity centre of the element from geosys
		double const* gc(m_ele->GetGravityCenter());
		this->eclgrid[i]->x_barycentre = gc[0];
		this->eclgrid[i]->y_barycentre = gc[1];
		this->eclgrid[i]->z_barycentre = gc[2];
		this->eclgrid[i]->volume = m_ele->GetVolume();
	}
	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                    Time: " << time << " seconds." << endl;

	return 1;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: CreateFaces
   Task: Create all Faces and connect them to nodes and elements
   Return: bool value if there occured an error
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
// TF commented since we can use sqrt (MathLib::sqrDist(pnt1, pnt2))
// this has three advantages:
// (1) we use existing code and do not invent the distance computation again
// (2) since the proposed function takes const arguments you can use
//     it in const methods, which is not the case with this method
// (3) imho the distance computation has nothing to do with the eclipse interface
//double CECLIPSEData::CalculateDistanceBetween2Points(double Point1[3], double Point2[3])
//{
//	double distance;
//
//	distance = sqrt(pow(Point1[0] - Point2[0],2) + pow(Point1[1] - Point2[1],2) + pow(Point1[2] - Point2[2],2));
//
//	if (distance >= 0)
//		return distance;
//	else
//		return -1;
//};

/*-------------------------------------------------------------------------
   GeoSys - Function: CreateFaces
   Task: Create all Faces and connect them to nodes and elements
   Return: bool value if there occured an error
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CreateFaces(void)
{
	CFaces* m_face = NULL;
	MeshLib::CElem* m_element = NULL;
	vector <long> element_indices;
	CFEMesh* m_msh = fem_msh_vector[0];
	Math_Group::vec <MeshLib::CNode*> element_nodes(8);
	//CNode* Node1=NULL;
	//CNode* Node2=NULL;
	//CNode* Node3=NULL;
	//CNode* Node4=NULL;
	vector <MeshLib::CNode*> vec_face_nodes;
	clock_t start,finish;
	double time;

	start = clock();

	//order of nodes for each element in geosys: 0..x1,y1,zoben; 1..x2,y1,zoben; 2..x2,y2,zoben; 3..x1,y2,zoben;
	//4..x1,y1,zu; 5..x2,y1,zu; 6..x2,y2,zu; 7..x1,y2,zu

	// in case of radial models i+...direction of radius, j+...direction of theta -> the faces are build with the same points as for rectangular models, but its not necessarily the same direction
	cout << "        CreateFaces() ";

	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
	{
		m_element = m_msh->ele_vector[i];
		m_element->GetNodes(element_nodes);

		//Test output
		//double *tmp_coord;
		//for (int j = 0; j < 8; j++) {
		//	tmp_coord = element_nodes[j]->GetCoordinates();
		//	cout << "Punkt " << j << " " << tmp_coord[0] << " " << tmp_coord[1] << " " << tmp_coord[2] << endl;
		//}
		//if ((i == 608) || (i == 645))
		//	cout << i << endl;
		//non-radial elements: Right face of the Element (i-direction) (Nodes: 1,5,2,6) -> with this order of the nodes the normal vector shows automatically to the positve i-direction
		//radial elements: Right face of the Element (i-direction) (Nodes: 6,5,2,1) -> with this order of the nodes the normal vector shows automatically to the positve i-direction
		m_face = new CFaces(int(this->Phases.size()));
		m_face->index = long(faces.size());
		m_face->model_axis = "I+";
		//store all nodes of the face in a vector
		vec_face_nodes.clear();
		if (Radialmodell == false)
		{
			vec_face_nodes.push_back(element_nodes[1]);
			vec_face_nodes.push_back(element_nodes[5]);
			vec_face_nodes.push_back(element_nodes[2]);
			vec_face_nodes.push_back(element_nodes[6]);
		}
		if (this->RadialModellIpos == true )
		{
			vec_face_nodes.push_back(element_nodes[6]);
			vec_face_nodes.push_back(element_nodes[5]);
			vec_face_nodes.push_back(element_nodes[2]);
			vec_face_nodes.push_back(element_nodes[1]);
		}

		//create the face
		if (m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
		                       vec_face_nodes[3]) == false)
		{
			cout << "Error at creating face: The given points form not a plane!" <<
			endl;
			system("Pause");
			return false;
		}
		//add connected elements to the face (first the left (x) element, second the right element)
		element_indices.clear();
		element_indices.push_back(i);
		if (this->eclgrid[i]->NeighbourElement[1] > -1)
			element_indices.push_back(this->eclgrid[i]->NeighbourElement[1]);
		m_face->SetElements(element_indices);
		//add the face to the vector
		faces.push_back(m_face);
		for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
		{
			//Connect nodes with faces
			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(
			        m_face->index);
			//Calculate distance between node and gravity centre of the face and store it in a vector
//			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->distance_to_connected_faces.push_back(this->CalculateDistanceBetween2Points(vec_face_nodes[j]->getData(), m_face->GetFaceGravityCentre()));
			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
			distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[
			                                                                    j
			                                                            ]->getData(),
			                                                            m_face->
			                                                            GetFaceGravityCentre())));
		}

		//non-radial elements: Upper face of the Element (j-direction) (Nodes: 3,2,7,6) -> with this order of the nodes the normal vector shows automatically to the positve j-direction
		//radial elements: Upper face of the Element (j-direction) (Nodes: 6,2,7,3) -> with this order of the nodes the normal vector shows automatically to the positve j-direction
		m_face = new CFaces(int(this->Phases.size()));
		m_face->index = long(faces.size());
		m_face->model_axis = "J+";
		//store all nodes of the face in a vector
		vec_face_nodes.clear();
		if (Radialmodell == false)
		{
			vec_face_nodes.push_back(element_nodes[3]);
			vec_face_nodes.push_back(element_nodes[2]);
			vec_face_nodes.push_back(element_nodes[7]);
			vec_face_nodes.push_back(element_nodes[6]);
		}
		if (this->RadialModellIpos == true )
		{
			vec_face_nodes.push_back(element_nodes[6]);
			vec_face_nodes.push_back(element_nodes[2]);
			vec_face_nodes.push_back(element_nodes[7]);
			vec_face_nodes.push_back(element_nodes[3]);
		}
		//create the face
		if (m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
		                       vec_face_nodes[3]) == false)
		{
			cout << "Error at creating face: The given points form not a plane!" <<
			endl;
			system("Pause");
			return false;
		}
		//add connected elements to the face (first the lower (y) element, second the upper element)
		element_indices.clear();
		element_indices.push_back(i);
		if (this->eclgrid[i]->NeighbourElement[3] > -1)
			element_indices.push_back(this->eclgrid[i]->NeighbourElement[3]);
		m_face->SetElements(element_indices);
		//add the face to the vector
		faces.push_back(m_face);
		for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
		{
			//Connect nodes with faces
			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(
			        m_face->index);
			//Calculate distance between node and gravity centre of the face and store it in a vector
			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
			distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[
			                                                                    j
			                                                            ]->getData(),
			                                                            m_face->
			                                                            GetFaceGravityCentre())));
		}

		//non-radial elements: Bottom face of the Element (k-direction) (Nodes: 4,5,7,6) -> with this order of the nodes the normal vector shows automatically to the positve k-direction
		//radial elements: Bottom face of the Element (k-direction) (Nodes: 6,5,7,4) -> with this order of the nodes the normal vector shows automatically to the positve k-direction
		m_face = new CFaces(int(this->Phases.size()));
		m_face->index = long(faces.size());
		m_face->model_axis = "K+";
		//store all nodes of the face in a vector
		vec_face_nodes.clear();
		if (Radialmodell == false)
		{
			vec_face_nodes.push_back(element_nodes[4]);
			vec_face_nodes.push_back(element_nodes[5]);
			vec_face_nodes.push_back(element_nodes[7]);
			vec_face_nodes.push_back(element_nodes[6]);
		}
		else
		{
			vec_face_nodes.push_back(element_nodes[6]);
			vec_face_nodes.push_back(element_nodes[5]);
			vec_face_nodes.push_back(element_nodes[7]);
			vec_face_nodes.push_back(element_nodes[4]);
		}
		//create the face
		if (m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1], vec_face_nodes[2],
		                       vec_face_nodes[3]) == false)
		{
			cout << "Error at creating face: The given points form not a plane!" <<
			endl;
			system("Pause");
			return false;
		}
		//add connected elements to the face (first the top (z) element, second the bottom element, because the order of z is in Eclipse opposit)
		element_indices.clear();
		element_indices.push_back(i);
		if (this->eclgrid[i]->NeighbourElement[5] > -1)
			element_indices.push_back(this->eclgrid[i]->NeighbourElement[5]);
		m_face->SetElements(element_indices);
		//add the face to the vector
		faces.push_back(m_face);
		for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
		{
			//Connect nodes with faces
			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.push_back(
			        m_face->index);
			//Calculate distance between node and gravity centre of the face and store it in a vector
			m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
			distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(vec_face_nodes[
			                                                                    j
			                                                            ]->getData(),
			                                                            m_face->
			                                                            GetFaceGravityCentre())));
		}

		//non-radial elements: Left face of the Element at the left (x) site of the grid (i-direction) (Nodes: 0,4,3,7) -> with this order of the nodes the normal vector shows automatically to the positve i-direction
		//radial elements: Left face of the Element at the left (x) site of the grid (i-direction) (Nodes: 7,4,3,0) -> with this order of the nodes the normal vector shows automatically to the positve i-direction
		//if (this->eclgrid[i]->column == 1){
		// if there is no element at the left site of the current element
		if (this->eclgrid[i]->NeighbourElement[0] == -1)
		{
			m_face = new CFaces(int(this->Phases.size()));
			m_face->index = long(faces.size());
			m_face->model_axis = "I-";
			//store all nodes of the face in a vector
			vec_face_nodes.clear();
			if (Radialmodell == false)
			{
				vec_face_nodes.push_back(element_nodes[0]);
				vec_face_nodes.push_back(element_nodes[4]);
				vec_face_nodes.push_back(element_nodes[3]);
				vec_face_nodes.push_back(element_nodes[7]);
			}
			if (this->RadialModellIpos == true )
			{
				vec_face_nodes.push_back(element_nodes[7]);
				vec_face_nodes.push_back(element_nodes[4]);
				vec_face_nodes.push_back(element_nodes[3]);
				vec_face_nodes.push_back(element_nodes[0]);
			}
			//create the face
			if (m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1],
			                       vec_face_nodes[2], vec_face_nodes[3]) == false)
			{
				cout <<
				"Error at creating face: The given points form not a plane!" <<
				endl;
				system("Pause");
				return false;
			}
			//add connected elements to the face
			element_indices.clear();
			element_indices.push_back(i);
			if (this->eclgrid[i]->NeighbourElement[0] > -1)
				return false;
			m_face->SetElements(element_indices);
			//add the face to the vector
			faces.push_back(m_face);
			for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
			{
				//Connect nodes with faces
				m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.
				push_back(m_face->index);
				//Calculate distance between node and gravity centre of the face and store it in a vector
				m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
				distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(
				                                                   vec_face_nodes[j
				                                                   ]->getData(),
				                                                   m_face->
				                                                   GetFaceGravityCentre())));
			}
		}

		//non-radial elements: Lower (y) face of the Element at the lower (y) site of the grid (j-direction) (Nodes: 0,1,4,5) -> with this order of the nodes the normal vector shows automatically to the positve j-direction
		//radial elements: Lower (y) face of the Element at the lower (y) site of the grid (j-direction) (Nodes: 5,1,4,0) -> with this order of the nodes the normal vector shows automatically to the positve j-direction
		//if (this->eclgrid[i]->row == 1){
		if (this->eclgrid[i]->NeighbourElement[2] == -1)
		{
			m_face = new CFaces(int(this->Phases.size()));
			m_face->index = long(faces.size());
			m_face->model_axis = "J-";
			//store all nodes of the face in a vector
			vec_face_nodes.clear();
			if (Radialmodell == false)
			{
				vec_face_nodes.push_back(element_nodes[0]);
				vec_face_nodes.push_back(element_nodes[1]);
				vec_face_nodes.push_back(element_nodes[4]);
				vec_face_nodes.push_back(element_nodes[5]);
			}
			if (this->RadialModellIpos == true )
			{
				vec_face_nodes.push_back(element_nodes[5]);
				vec_face_nodes.push_back(element_nodes[1]);
				vec_face_nodes.push_back(element_nodes[4]);
				vec_face_nodes.push_back(element_nodes[0]);
			}
			//create the face
			if (m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1],
			                       vec_face_nodes[2], vec_face_nodes[3]) == false)
			{
				cout <<
				"Error at creating face: The given points form not a plane!" <<
				endl;
				system("Pause");
				return false;
			}
			//add connected elements to the face
			element_indices.clear();
			element_indices.push_back(i);
			if (this->eclgrid[i]->NeighbourElement[2] > -1)
				return false;
			m_face->SetElements(element_indices);
			//add the face to the vector
			faces.push_back(m_face);
			for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
			{
				//Connect nodes with faces
				m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.
				push_back(m_face->index);
				//Calculate distance between node and gravity centre of the face and store it in a vector
				m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
				distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(
				                                                   vec_face_nodes[j
				                                                   ]->getData(),
				                                                   m_face->
				                                                   GetFaceGravityCentre())));
			}
		}

		//non-radial elements: Top (z) face of the Element at the top (z) site of the grid (k-direction) (Nodes: 0,1,3,2) -> with this order of the nodes the normal vector shows automatically to the positve k-direction
		//radial elements: Top (z) face of the Element at the top (z) site of the grid (k-direction) (Nodes: 2,1,3,0) -> with this order of the nodes the normal vector shows automatically to the positve k-direction
		//if (this->eclgrid[i]->layer == 1 /* SB->BG :changed - please check "if (this->eclgrid[i]->layer == this->layers"*/){
		if (this->eclgrid[i]->NeighbourElement[4] == -1)
		{
			m_face = new CFaces(int(this->Phases.size()));
			m_face->index = long(faces.size());
			m_face->model_axis = "K-";
			//store all nodes of the face in a vector
			vec_face_nodes.clear();
			if (Radialmodell == false)
			{
				vec_face_nodes.push_back(element_nodes[0]);
				vec_face_nodes.push_back(element_nodes[1]);
				vec_face_nodes.push_back(element_nodes[3]);
				vec_face_nodes.push_back(element_nodes[2]);
			}
			else
			{
				vec_face_nodes.push_back(element_nodes[2]);
				vec_face_nodes.push_back(element_nodes[1]);
				vec_face_nodes.push_back(element_nodes[3]);
				vec_face_nodes.push_back(element_nodes[0]);
			}
			//create the face
			if (m_face->CreateFace(vec_face_nodes[0], vec_face_nodes[1],
			                       vec_face_nodes[2], vec_face_nodes[3]) == false)
			{
				cout <<
				"Error at creating face: The given points form not a plane!" <<
				endl;
				system("Pause");
				return false;
			}
			//add connected elements to the face
			element_indices.clear();
			element_indices.push_back(i);
			if (this->eclgrid[i]->NeighbourElement[4] > -1)
				return false;
			m_face->SetElements(element_indices);
			//add the face to the vector
			faces.push_back(m_face);
			for (unsigned long j = 0; j < vec_face_nodes.size(); j++)
			{
				//Connect nodes with faces
				m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->connected_faces.
				push_back(m_face->index);
				//Calculate distance between node and gravity centre of the face and store it in a vector
				m_msh->nod_vector[vec_face_nodes[j]->GetIndex()]->
				distance_to_connected_faces.push_back(sqrt(MathLib::sqrDist(
				                                                   vec_face_nodes[j
				                                                   ]->getData(),
				                                                   m_face->
				                                                   GetFaceGravityCentre())));
			}
		}
	}

	// Test output faces
	//double *n_vec;
	//for(long i=0;i<faces.size();i++){
	//	cout << "  Face " << i << "   " ;
	//	for(int j=0; j< faces[i]->connected_elements.size();j++)
	//		cout << faces[i]->connected_elements[j] << " - " ;
	//	cout << "           " << faces[i]->model_axis ;
	//	cout << "            n_vec: " ;
	//	n_vec = faces[i]->PlaneEquation->GetNormalVector();
	//	for(int j=0;j<3;j++) cout << n_vec[j] << ", ";
	//	cout << "     Area: " << faces[i]->face_area;
	//	cout << endl;
	//}

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                        Time: " << time << " seconds." << endl;

	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CreateFacesAtElements
   Task: For each element, find all faces and store them in connected_faces
   Return: bool value if there occured an error
   Programming: 09/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::ConnectFacesToElements(void)
{
	//CECLIPSEBlock*m_block=NULL;
	CFaces* m_face = NULL;
	long ind_block, ind_face;
	vector <long> faces_at_block;
	clock_t start,finish;
	double time;

	start = clock();

	cout << "        ConnectFacesToElements() ";

	// go through all faces and add them to the eclipse blocks
	for(unsigned long i = 0; i < this->faces.size(); i++)
	{
		m_face = this->faces[i];
		ind_face = m_face->index;
		for(unsigned long j = 0; j < m_face->connected_elements.size(); j++)
		{
			ind_block = m_face->connected_elements[j];
			this->eclgrid[ind_block]->connected_faces.push_back(ind_face);
		}
	}
	// Test output
	//for(unsigned long i=0;i<eclgrid.size();i++){
	//	cout << "    Faces at block " << eclgrid[i]->index << " :  ";
	//	for(int j=0;j<eclgrid[i]->connected_faces.size();j++) cout << eclgrid[i]->connected_faces[j] << ", ";
	//	cout << endl;
	//}

	for(long i = 0; i < long(eclgrid.size()); i++)
	{
		// Test output: all faces at one block
		if(eclgrid[i]->connected_faces.size() != 6)
			cout << " Error - not 6 faces at one block " << endl;
		// initialize vector faces_at_block
		faces_at_block.clear();
		for(unsigned long j = 0; j < 6; j++)
			faces_at_block.push_back(-1);
		for(unsigned long j = 0; j < eclgrid[i]->connected_faces.size(); j++)
		{
			//cout << eclgrid[i]->connected_faces[j] << ", ";
			// Order faces +i, -i, +j. -j, +k , -k
			m_face = faces[eclgrid[i]->connected_faces[j]]; // get face
			if(m_face->model_axis == "I+")
			{
				if(m_face->connected_elements[0] == i)
					faces_at_block[0] = m_face->index;  // this is +i for this element
				else
					faces_at_block[1] = m_face->index;  // this is -i for this element
			}
			if(m_face->model_axis == "I-")
			{
				if(m_face->connected_elements[0] == i)
					faces_at_block[1] = m_face->index;  // for this element it is i-
				else
					faces_at_block[0] = m_face->index;  // for this element it is i+
			}
			if(m_face->model_axis == "J+")
			{
				if(m_face->connected_elements[0] == i)
					faces_at_block[2] = m_face->index;
				else
					faces_at_block[3] = m_face->index;
			}
			if(m_face->model_axis == "J-")
			{
				if(m_face->connected_elements[0] == i)
					faces_at_block[3] = m_face->index;
				else
					faces_at_block[2] = m_face->index;
			}
			if(m_face->model_axis == "K+")
			{
				if(m_face->connected_elements[0] == i)
					faces_at_block[4] = m_face->index;
				else
					faces_at_block[5] = m_face->index;
			}
			if(m_face->model_axis == "K-")
			{
				if(m_face->connected_elements[0] == i)
					faces_at_block[5] = m_face->index;
				else
					faces_at_block[4] = m_face->index;
			}
		}
		// copy sorted faces to ecl block
		for(unsigned long j = 0; j < eclgrid[i]->connected_faces.size(); j++)
			eclgrid[i]->connected_faces[j] = faces_at_block[j];
		//// Test output
		//cout << "    Faces at block " << eclgrid[i]->index << " :  ";
		//for(j=0;j<eclgrid[i]->connected_faces.size();j++) cout << eclgrid[i]->connected_faces[j] << ", ";
		//cout << endl;
	}
	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "             Time: " << time << " seconds." << endl;

	return true;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: GetFlowForFaces
   Task: Get's the flow from the corresponding elements
   Return: true or false
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::GetFlowForFaces(int phase_index)
{
	//CFEMesh* m_msh = fem_msh_vector[0];
	CFaces* m_face = NULL;
	long number_elements;
	long element_index;
	long variable_index;
	//WW long time_index;
	double vz = 1.0; // Sign
	clock_t start,finish;
	double time;

	start = clock();

	cout << "        GetFlowForFaces()";
	for (long i = 0; i < long(this->faces.size()); i++)
	{
		m_face = this->faces[i];
		//set initial values
		//if (i == 4225) {
		//	cout << i;
		//}
		m_face->phases[phase_index]->q_norm = 0.0;

		number_elements = long(m_face->connected_elements.size());

		if ((number_elements < 1) || (number_elements > 2))
		{
			cout << "There is an error in connecting faces and elements" << endl;
			system("Pause");
			exit(0);
		}

		// right face (x)
		if (m_face->model_axis == "I+")
		{
			element_index = m_face->connected_elements[0];
			//WW		time_index = 0;
			if (this->RadialModellJpos == true)
			{
				if (this->Phases[phase_index] == "WATER")
					variable_index = this->GetVariableIndex("FLOWATJ+");
				else
				{
					if (this->Phases[phase_index] == "GAS")
						variable_index = this->GetVariableIndex("FLOGASJ+");
					else
					{
						if (this->Phases[phase_index] == "OIL")
							variable_index = this->GetVariableIndex(
							        "FLOOILJ+");
						else
						{
							cout <<
							"This phase is not considered yet in GetFlowForFaces"
							     << endl;
							system("Pause");
							exit(0);
						}
					}
				}
			}
			else
			{
				if (this->Phases[phase_index] == "WATER")
					variable_index = this->GetVariableIndex("FLOWATI+");
				else
				{
					if (this->Phases[phase_index] == "GAS")
						variable_index = this->GetVariableIndex("FLOGASI+");
					else
					{
						if (this->Phases[phase_index] == "OIL")
							variable_index = this->GetVariableIndex(
							        "FLOOILI+");
						else
						{
							cout <<
							"This phase is not considered yet in GetFlowForFaces"
							     << endl;
							system("Pause");
							exit(0);
						}
					}
				}
			}
			if (variable_index < 0)
			{
				cout << "There are no variables!" << endl;
				system("Pause");
				exit(0);
			}
			//calculate q [m�/m�.s] for Geosys from Q [m�/d] from Eclipse
			m_face->phases[phase_index]->q_norm =
			        (this->Data[element_index][variable_index] /
			         m_face->face_area) / 86400.0;
			m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
			//cout << " Flow through face " << m_face->index << " element " << element_index << " :  " << this->Data[element_index][time_index][variable_index] << " ,  " << m_face->q_norm << endl;
		}
		// upper face (y)
		if (m_face->model_axis == "J+")
		{
			element_index = m_face->connected_elements[0];
			if (this->RadialModellJpos == true)
			{
				if (this->Phases[phase_index] == "WATER")
					variable_index = this->GetVariableIndex("FLOWATI+");
				else
				{
					if (this->Phases[phase_index] == "GAS")
						variable_index = this->GetVariableIndex("FLOGASI+");
					else
					{
						if (this->Phases[phase_index] == "OIL")
							variable_index = this->GetVariableIndex(
							        "FLOOILI+");
						else
						{
							cout <<
							"This phase is not considered yet in GetFlowForFaces"
							     << endl;
							system("Pause");
							exit(0);
						}
					}
				}
			}
			else
			{
				if (this->Phases[phase_index] == "WATER")
					variable_index = this->GetVariableIndex("FLOWATJ+");
				else
				{
					if (this->Phases[phase_index] == "GAS")
						variable_index = this->GetVariableIndex("FLOGASJ+");
					else
					{
						if (this->Phases[phase_index] == "OIL")
							variable_index = this->GetVariableIndex(
							        "FLOOILJ+");
						else
						{
							cout <<
							"This phase is not considered yet in GetFlowForFaces"
							     << endl;
							system("Pause");
							exit(0);
						}
					}
				}
			}
			//WW		time_index = 0;
			if (variable_index < 0)
			{
				cout << "There are no variables!" << endl;
				system("Pause");
				exit(0);
			}
			//calculate q [m�/m�.s] for Geosys from Q [m�/d] from Eclipse
			m_face->phases[phase_index]->q_norm =
			        (this->Data[element_index][variable_index] /
			         m_face->face_area) / 86400.0;
			m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
			//cout << " Flow through face " << m_face->index <<  " element " << element_index << " :  " << this->Data[element_index][time_index][variable_index] << " ,  " << m_face->q_norm << endl;
		}
		// bottom face (z) (Eclipse has the opposite z direction -> k = +
		if (m_face->model_axis == "K+")
		{
			element_index = m_face->connected_elements[0];
			if (this->Phases[phase_index] == "WATER")
				variable_index = this->GetVariableIndex("FLOWATK+");
			else
			{
				if (this->Phases[phase_index] == "GAS")
					variable_index = this->GetVariableIndex("FLOGASK+");
				else
				{
					if (this->Phases[phase_index] == "OIL")
						variable_index = this->GetVariableIndex("FLOOILK+");
					else
					{
						cout <<
						"This phase is not considered yet in GetFlowForFaces"
						     << endl;
						system("Pause");
						exit(0);
					}
				}
			}
			//WW		time_index = 0;
			if (variable_index < 0)
			{
				cout << "There are no variables!" << endl;
				system("Pause");
				exit(0);
			}
			//calculate q [m�/m�.s] for Geosys from Q [m�/d] from Eclipse
			m_face->phases[phase_index]->q_norm =
			        (this->Data[element_index][variable_index] /
			         m_face->face_area) / 86400.0;
			m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
			//cout << " Flow through face " << m_face->index <<  " element " << element_index << " :  " << this->Data[element_index][time_index][variable_index] << " ,  " << m_face->q_norm << endl;
		}
		// left, lower and top faces -> now flow over these faces
		if ((m_face->model_axis == "I-") || (m_face->model_axis == "J-") ||
		    (m_face->model_axis == "K-"))
		{
			m_face->phases[phase_index]->q_norm = 0;
			m_face->Calculate_components_of_a_vector(0, phase_index, Radialmodell);
		}
		//if ((element_index == 609) || (element_index == 610) || (element_index == 611) || (element_index == 646) || (element_index == 647) || (element_index == 648) || (element_index == 683) || (element_index == 684) || (element_index == 685)  || (element_index == 720)  ||  (element_index == 721)  ||  (element_index == 722))
		//	cout << "Zentrum" << endl;

		//
		//Get additional flow for border cells that contain boundary conditions
		if (m_face->connected_elements.size() == 1)
		{
			long ele_index = m_face->connected_elements[0];
			//check if there is one BC in the element
			if (this->eclgrid[ele_index]->ConnectedBoundaryCondition.size() > 1)
			{
				cout <<
				"There is more than 1 boundary condition assigned to the cell " <<
				ele_index << endl;
				system("Pause");
				exit(0);
			}
			if (this->eclgrid[ele_index]->ConnectedBoundaryCondition.size() == 1)
			{
				long bc_index =
				        this->eclgrid[ele_index]->ConnectedBoundaryCondition[0];
				//check if the boundary inflow fits to the current face
				if (m_face->model_axis == this->BC[bc_index]->boundary_position)
				{
					//check if q is 0
					if (m_face->phases[phase_index]->q_norm == 0)
					{
						//calculate q [m�/m�.s] for Geosys from Q [m�/d] from Eclipse
						//cout << " BC Flow through face " << m_face->index << " :  " << this->BC[bc_index]->value[0] << endl;
						vz = 1.0;
						//if(this->BC[bc_index]->value[0] < 0){ // outflow out of model area
						//The values in I+, J+ and K+ direction have to be counted with the opposit direction
						// because the flow from BC is given positive if fluid flows from the BC into the model
						if(m_face->model_axis == "I+")
							vz = -1.0;
						if(m_face->model_axis == "J+")
							vz = -1.0;
						if(m_face->model_axis == "K+")
							vz = -1.0;
						//}
						m_face->phases[phase_index]->q_norm =
						        (this->BC[bc_index]->value[0] * vz /
						         m_face->face_area) / 86400.0;
						m_face->Calculate_components_of_a_vector(
						        0,
						        phase_index,
						        Radialmodell);
						//cout << " Flow through face " << m_face->index <<  " element " << element_index << " :  " << this->BC[bc_index]->value[0] << " ,  " << m_face->q_norm << endl;
					}
					else
					{
						cout <<
						"There is already a flow assigned to the boundary face, which shouldn't happen "
						     << m_face->index << endl;
						system("Pause");
						exit(0);
					}
				}
			}
		}

		//cout << " Face " << m_face->index << ", " << m_face->connected_elements[0] ;
		//if(m_face->connected_elements.size() > 1) cout << ", " << m_face->connected_elements[1];
		//cout << ":           " << m_face->q_norm << ":   " << m_face->q[0] << ", " << m_face->q[1] << ", " << m_face->q[2] << endl;
	}

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                     Time: " << time << " seconds." << endl;

	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GetVelForFaces
   Task: Get's the velocities from the corresponding elements
   derived from GetFlowForFaces()
   Return: true or false
   Programming: 09/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::GetVelForFaces(void)
{
	//CFEMesh* m_msh = fem_msh_vector[0];
	CFaces* m_face = NULL;
	long number_elements;
	long element_index;
	long variable_index;
	//WW long time_index;

	for (unsigned long i = 0; i < this->faces.size(); i++)
	{
		m_face = this->faces[i];
		number_elements = long(m_face->connected_elements.size());

		if ((number_elements < 1) || (number_elements > 2))
			return false;

		// right face (x)
		if (m_face->model_axis == "I+")
		{
			element_index = m_face->connected_elements[0];
			//WW		time_index = 0;
			variable_index = this->GetVariableIndex("VELWATI+");
			if (variable_index < 0)
				return false;
			m_face->v_norm = this->Data[element_index][variable_index];
			//m_face->Calculate_components_of_a_vector(1);
		}
		// upper face (y)
		if (m_face->model_axis == "J+")
		{
			element_index = m_face->connected_elements[0];
			variable_index = this->GetVariableIndex("VELWATJ+");
			//WW		time_index = 0;
			if (variable_index < 0)
				return false;
			m_face->v_norm = this->Data[element_index][variable_index];
			//m_face->Calculate_components_of_a_vector(1);
		}
		// bottom face (z) (Eclipse has the opposite z direction -> k = +
		if (m_face->model_axis == "K+")
		{
			element_index = m_face->connected_elements[0];
			variable_index = this->GetVariableIndex("VELWATK+");
			//WW		time_index = 0;
			if (variable_index < 0)
				return false;
			m_face->v_norm = this->Data[element_index][variable_index];
			//m_face->Calculate_components_of_a_vector(1);
		}
		// left, lower and top faces -> now flow over these faces
		if ((m_face->model_axis == "I-") || (m_face->model_axis == "J-") ||
		    (m_face->model_axis == "K-"))
			m_face->v_norm = 0;
		//m_face->Calculate_components_of_a_vector(1);
		//

		cout << " Face " << m_face->index << ", " << m_face->connected_elements[0];
		if(m_face->connected_elements.size() > 1)
			cout << ", " << m_face->connected_elements[1];
		//cout << ": " << m_face->q_norm << endl;
		//if(fabs(m_face->v_norm) > MKleinsteZahl) cout << " q/v: " << m_face->q_norm/m_face->face_area/0.5 /m_face->v_norm ;
		//cout << endl;
	}
	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CalcBlockBudegt
   Task: Summ all budget terms from the faces for each block - corresponds to a mass balance in steady state
   This function is for checking purposes only
   Return: true or false
   Programming: 09/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CalcBlockBudget(int phase_index)
{
	CECLIPSEBlock* m_block = NULL;
	CFaces* m_face = NULL;
	double flow, flow_face, flow_max;
	double max_error = 0;

	// summ budget flow terms at feces for all blocks
	for(unsigned i = 0; i < this->eclgrid.size(); i++)
	{
		m_block = this->eclgrid[i];
		flow = flow_max = 0.0;
		//cout << "  Budget for block " << i << ": " << endl ;
		for(unsigned int j = 0; j < m_block->connected_faces.size(); j++)
		{
			m_face = faces[m_block->connected_faces[j]];
			flow_face = m_face->phases[phase_index]->q_norm * m_face->face_area;
			//cout << "flow_face " << flow_face << endl;
			if(fabs(flow_face) > flow_max)
				flow_max = fabs(flow_face);
			if((j == 1) || (j == 3) || (j == 5)) // minus faces: -i, -j or -k
				flow_face *= -1.0;
			//cout << " face " << m_face->index << ": " << flow_face << " ;   " << endl;
			flow += flow_face;
		}
		//cout << " total flow: " << flow << endl;
		if((fabs(flow / flow_max) > 1.0e-3) && (fabs(flow) > 1.0e-10))
		{
			if (max_error < abs(flow / flow_max))
				max_error = abs(flow / flow_max);
			cout << " Error in budget for block " << i << " :  Sum_flow: " << flow <<
			", max_flow : " << flow_max << ", rel_error: " << flow / flow_max <<
			", max_error: " << max_error << endl;
		}
		//cout << " max_error: " << max_error << endl;
		//if (max_error > 0.0043)
		//	cout << endl;
	}

	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateDataFromFacesToNodes
   Task: Interpolates data like phase flow velocity from faces (Eclipse) to nodes
   Return: true or false
   Programming: 09/2009 BG / SB
   Modification:
   -------------------------------------------------------------------------*/
//bool CECLIPSEData::MakeNodeVector(CRFProcess *m_pcs, std::string path, int timestep, int phase_index)
bool CECLIPSEData::MakeNodeVector(void)
{
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	//CFaces *m_face=NULL;
	//WW double weights_xyz[3];
	CPointData_ECL* m_NodeData = NULL;
	m_NodeData = new CPointData_ECL;
	vector <double> temp_q;

	//double coord_v_x[3] = {1.0,0.0,0.0};
	//double coord_v_y[3] = {0.0,1.0,0.0};
	//double coord_v_z[3] = {0.0,0.0,1.0};
	//double * normal_vec_face;
	//double distance;
	//double weight;
	//double sum_weights;
	//double val;

	for (int i = 0; i < 3; i++)
		temp_q.push_back(-1.0E+99);

	/* Go through GeoSys mesh, and creates CPointData - then Nodes in
	   CECLIPSEData::NodeData have the same order as in the Geosys mesh,
	   and can be used directly for interpolation to gauss points*/
	//this->NodeData.clear();
	if(this->NodeData.size() < 1)
		for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
		{
			// create new instance of CPointData
			m_NodeData = new CPointData_ECL (m_msh->nod_vector[i]->getData());
			//Get the node
// TF		m_node = m_msh->nod_vector[i];
// TF		m_NodeData->x = m_node->X();
// TF		m_NodeData->y = m_node->Y();
// TF		m_NodeData->z = m_node->Z();
			m_NodeData->phase_pressure.resize(3);
			m_NodeData->phase_saturation.resize(3);
			m_NodeData->phase_density.resize(3);
			//for (int j = 0; j < int(this->Phases.size()); j++)
			//	m_NodeData->q.push_back(temp_q);

			//Set variable to zero
			m_NodeData->pressure = -1.0E+99;
			m_NodeData->CO2inLiquid = -1.0E+99;
			m_NodeData->deltaDIC = -1.0E+99;
			for (long k = 0; k < 3; k++)
			{
				//m_NodeData->Flow[k] = 0.0;
				m_NodeData->phase_pressure[k] = -1.0E+99;
				m_NodeData->phase_saturation[k] = -1.0E+99;
				m_NodeData->phase_density[k] = -1.0E+99;
				//WW		weights_xyz[k] = 0.0;
			}

			// transfer Data to node
			this->NodeData.push_back(m_NodeData);

			// Test output
//		cout << " Node " << i << " (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", "<< m_NodeData->z << ") ";
//		for (long k = 0; k < 3; k++){
//			cout << m_NodeData->phase_density[k] << " " << m_NodeData->phase_pressure[k] << " " << m_NodeData->phase_saturation[k] << " ";
//		}
//		cout << endl;
//		// Test output
//		cout << " Node " << i << " (X,Y,Z): (" << this->NodeData[this->NodeData.size()-1]->x << ", " << this->NodeData[this->NodeData.size()-1]->y << ", "<< this->NodeData[this->NodeData.size()-1]->z << ") ";
//		for (long k = 0; k < 3; k++){
//			cout << this->NodeData[this->NodeData.size()-1]->phase_density[k] << " " << this->NodeData[this->NodeData.size()-1]->phase_pressure[k] << " " << this->NodeData[this->NodeData.size()-1]->phase_saturation[k] << " ";
//		}
//		cout << endl;
		}
	else
		for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
		{
			// Set variable to zero
			this->NodeData[i]->pressure = -1.0E+99;
			this->NodeData[i]->CO2inLiquid = -1.0E+99;
			for (long k = 0; k < 3; k++)
			{
				//this->NodeData[i]->Flow[k] = 0.0;
				this->NodeData[i]->phase_pressure[k] = -1.0E+99;
				this->NodeData[i]->phase_saturation[k] = -1.0E+99;
				this->NodeData[i]->phase_density[k] = -1.0E+99;
			}
		}

//				for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++) {
//	                                        //	//Get the node
//					m_NodeData = this->NodeData[i];
//					m_node = m_msh->nod_vector[i];
//					m_NodeData->x = m_node->X();
//					m_NodeData->y = m_node->Y();
//					m_NodeData->z = m_node->Z();
//					//sum_weights = 0;
//					// sum the distance weighted data from each connected face
//					for (unsigned int j = 0; j < m_node->connected_faces.size(); j++) {
//					  distance = weight = 0;
//					  m_face = this->faces[m_node->connected_faces[j]];
//					  distance = m_node->distance_to_connected_faces[j];
//					  //Weight of each face depending on distance
//					  weight = (1.0 / distance);
//					  // Sum of weights
//					  //sum_weights += weight;
//					  normal_vec_face = m_face->PlaneEquation->GetNormalVector();
//					  // Go through all three coordinates x, y, z and check, if face is perpendicular to axis
//					  for (int k = 0; k < 3; k++){
//						  if(k == 0) val = fabs(PointProduction(normal_vec_face,coord_v_x));
//						  if(k == 1) val = fabs(PointProduction(normal_vec_face,coord_v_y));
//						  if(k == 2) val = fabs(PointProduction(normal_vec_face,coord_v_z));
//						  if(val > MKleinsteZahl){ // face not perpendicular to ccordinate axis k
//							//m_NodeData->Flow[k] += m_face->phases[phase_index]->q[k] * weight;
//							weights_xyz[k] += weight;
//							//cout << " Node " << i << " contributed by face " << m_face->index << ", " << m_face->model_axis << "  , " << m_face->q_norm << ": ";
//							//for(int mm=0;mm<3;mm++) cout << m_face->q[mm] << ", ";
//							//cout << endl;
//						  }
//					  }
//					}
//
//					// normalize weighted sum by sum_of_weights sum_w
//					for (int k = 0; k < 3; k++) {
//					 // if(weights_xyz[k] > 0.0)
//						////m_NodeData->Flow[k] = m_NodeData->Flow[k] / weights_xyz[k];
//					 // else{
//						////m_NodeData->Flow[k] = m_NodeData->Flow[k];
//						//cout << " Warning - no faces for direction (I=0,J=1,K=2): " << k << " at node " << i << " with non-zero contributions! " << endl;
//					 // }
//					}
//				}
//				// TEst Output
//				std::string tempstring;
//				ostringstream temp;
//				vector <string> vec_string;
//				double v_geosys[3];
//				vec_string.push_back("Node; X; Y; Z; v_x_ECL; v_y_ECL; v_z_ECL; v_ECL; v_x_Geos; v_y_Geos; v_z_Geos; v_Geos");
//
//				for(unsigned long i=0; i< this->NodeData.size();i++){
//					m_NodeData = this->NodeData[i];
//					//cout << " NodeData["<< i << "]:  (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", " << m_NodeData->z << "):        " ;
//					tempstring="";
//					temp.str(""); temp.clear(); temp << i; tempstring = temp.str();
//					temp.str(""); temp.clear(); temp << m_NodeData->x; tempstring += "; " +temp.str();
//					temp.str(""); temp.clear(); temp << m_NodeData->y; tempstring += "; " + temp.str();
//					temp.str(""); temp.clear(); temp << m_NodeData->z; tempstring += "; " + temp.str();
//					for(int k=0;k< 3;k++)  {
//						//temp.str(""); temp.clear(); temp << m_NodeData->Flow[k]; tempstring += "; " + temp.str();
//						//cout << m_NodeData->Data_separated[k] << " ";
//					}
//					//cout << endl;
//					//calculate v
////					temp.str(""); temp.clear(); temp << sqrt(pow(m_NodeData->Flow[0],2) + pow(m_NodeData->Flow[1],2) + pow(m_NodeData->Flow[2], 2)); tempstring += "; " + temp.str();
//
//					//get node velocity of geosys
//					InterpolateGeosysVelocitiesToNodes(m_pcs, v_geosys, i);
//					for(int k=0;k< 3;k++)  {
//						temp.str(""); temp.clear(); temp << v_geosys[k]; tempstring += "; " + temp.str();
//					}
//					//calculate v
//					temp.str(""); temp.clear(); temp << sqrt(pow(v_geosys[0], 2) + pow(v_geosys[1], 2) + pow(v_geosys[2], 2)); tempstring += "; " + temp.str();
//					//write string
//					vec_string.push_back(tempstring);
//				}
//				int position = int(path.find_last_of("\\"));
//				path = path.substr(0,position);
//				position = int(path.find_last_of("\\"));
//				path = path.substr(0,position);
//				temp.str(""); temp.clear(); temp << timestep; tempstring = temp.str();
//				std::string aus_file = path + "\\CheckVelocityAtNodes_" + tempstring + ".csv";
//				ofstream aus;
//				aus.open(aus_file.data(),ios::out);
//				for (unsigned long i = 0; i < vec_string.size(); i++) {
//					aus << vec_string[i] << endl;
//				}
//				aus.close();

	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateDataFromFacesToNodes2
   Task: Interpolates phase velocity from faces (Eclipse) to nodes
   Return: true or false
   Programming: 10/2009 BG / SB
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::InterpolateDataFromFacesToNodes(long ele_nr,
                                                   double* n_vel_x,
                                                   double* n_vel_y,
                                                   double* n_vel_z,
                                                   int phase_index)
{
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	CFaces* m_face = NULL;
	MeshLib::CNode* m_node = NULL;
	double distance;
	double weight, weights_xyz[3];
	//double sum_weights;
	double val = 0;
	double coord_v_x[3] = {1.0,0.0,0.0};
	double coord_v_y[3] = {0.0,1.0,0.0};
	double coord_v_z[3] = {0.0,0.0,1.0};
	double data_separated[3];
	double* normal_vec_face;
	Math_Group::vec <long>nod_index(8);
	long nod_ind;
	bool choose;
	int face_flag = 1; // option: 0...take all faces at one node, 1...take 2 faces: only the faces which belong to the element
	//todo eingabeoption

	//cout << " InterpolateDataFromFacesToNodes2 " << endl << endl;

	/* Go through GeoSys mesh, and get each element (=Eclipse block). For each block, look at all corner points
	   (=Nodes) of the element, as phase velocities are needed there. For each node, look at all faces connected to the node.
	   From each face, get the flow accross the face, and make a weighted sum in the nodes.
	   Faces are not considered, if the are perpendicular to the flow direction, as they then would always contribute zero
	   Faces can be choosen with two options:
	   - All faces at one node, then also information of phase velocities from outside the element considered is
	   taken into account.
	   - Or only faces which are directly connected to the element considered.
	   The resulting phase velocities are stored component wise in the vectors n_vel_xyz and passed back.
	 */

	// Get element for which to calculate velocities
	m_ele = m_msh->ele_vector[ele_nr]; // get element
	m_ele->GetNodeIndeces(nod_index); // get node indexes of nodes on element

	//if ((ele_nr == 646) || (ele_nr == 647) || (ele_nr == 648) || (ele_nr == 683) || (ele_nr == 684) || (ele_nr == 685)  || (ele_nr == 720)  ||  (ele_nr == 721)  ||  (ele_nr == 722))
	//	cout << "Zentrum" << endl;

	for (long i = 0; i < long(nod_index.Size()); i++) // go through list of connected nodes of the element
	{ //Get the connected node
		nod_ind = nod_index[i];
		m_node = m_msh->nod_vector[nod_ind]; // m_node = Knoten, der gerade behandelt wird--> siehe Koordinaten
		//cout << " Get Node " << nod_ind << " on element " << ele_nr << endl;
		// Set variable to zero, initialization
		for (long k = 0; k < 3; k++)
		{
			data_separated[k] = 0.0;
			weights_xyz[k] = 0.0;
		}

		// sum the distance weighted data from each connected face
		for (unsigned int j = 0; j < m_node->connected_faces.size(); j++)
		{
			distance = weight = 0.0;
			m_face = this->faces[m_node->connected_faces[j]]; // Get face on node
			// Consider only if face is on element
			choose = false;
			if(m_face->connected_elements[0] == ele_nr)
				choose = true;         //check if face is part of the element
			if(m_face->connected_elements.size() > 1)
				if(m_face->connected_elements[1] == ele_nr)
					choose = true;
			if (face_flag == 0)
				choose = true;

			// for radial flow model
			if(this->Radial_I == true)
				// skip all J faces
				if(m_face->model_axis.find("J") == 0)
					choose = false;
			if(this->Radial_J == true)
				// skip all I faces
				if(m_face->model_axis.find("I") == 0)
					choose = false;

			if(choose)
			{
				//Testoutput
				//cout << endl;
				//cout << " Node " << nod_ind << " contributed by face " << m_face->index << ", " << m_face->model_axis << "  , " << m_face->phases[phase_index]->q_norm << ": ";
				//for(int mm=0;mm<3;mm++) cout << m_face->phases[phase_index]->q[mm] << ", ";
				//cout << endl;

				distance = m_node->distance_to_connected_faces[j];
				//Weight of each face depending on distance
				weight = (1.0 / distance);
				// Sum of weights
				//sum_weights += weight;
				normal_vec_face = m_face->PlaneEquation->GetNormalVector();
				// Go through all three coordinates x, y, z and check, if face is perpendicular to axis
				for (int k = 0; k < 3; k++)
				{
					if(k == 0)
						val =
						        fabs(PointProduction(normal_vec_face,
						                             coord_v_x));
					if(k == 1)
						val =
						        fabs(PointProduction(normal_vec_face,
						                             coord_v_y));
					if(k == 2)
						val =
						        fabs(PointProduction(normal_vec_face,
						                             coord_v_z));
					if(val > MKleinsteZahl) // face not perpendicular to ccordinate axis k
					{
						data_separated[k] +=
						        m_face->phases[phase_index]->q[k] * weight;
						weights_xyz[k] += weight;
					}
				}
			} // end if face is on element
		} // end loop faces at one node of element considered

		// normalize weighted sum by sum_of_weights sum_w
		for (int k = 0; k < 3; k++)
		{
			if(weights_xyz[k] > 0.0)
				data_separated[k] = data_separated[k] / weights_xyz[k];
			else
			{
				data_separated[k] = data_separated[k];
				if((k == 1) && (this->Radial_I == true))
				{
				}          // do nothing, Radial model perpendicular x axis
				else if ((k == 0) && (this->Radial_J == true))
				{
				}                // do nothing, Radial model perpendicular y axis
				else
					cout <<
					" Warning - no faces for direction (I=0,J=1,K=2): " << k <<
					" at node " << i << " with non-zero contributions! " <<
					endl;
			}
		}
		// transfer Data to node
		n_vel_x[i] = data_separated[0];
		n_vel_y[i] = data_separated[1];
		n_vel_z[i] = data_separated[2];

		//Test output
		//cout << " Node " << nod_ind << " (X,Y,Z): (" << m_node->X() << ", " << m_node->Y() << ", "<< m_node->Z() << ") " ;
		//cout << data_separated[0] << " " << data_separated[1] << " " << data_separated[2] << " " << endl;
	} // end loop connected nodes

	// Test Output
	//for(long i=0; i< nod_index.Size();i++){
	//	cout << " Node "<< nod_index[i] << ":  (X,Y,Z): (" << m_msh->nod_vector[nod_index[i]]->X() << ", " << m_msh->nod_vector[nod_index[i]]->Y() << ", "<< m_msh->nod_vector[nod_index[i]]->Z() << "):        " ;
	//	cout << n_vel_x[i] << " " << n_vel_y[i] << " " << n_vel_z[i] << " " << endl;
	//}/* */
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateDataFromBlocksToNodes
   Task: Interpolates data like phase pressure or saturation from blocks (Eclipse) to nodes (GeoSys)
   Return: true or false
   Programming: 10/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::InterpolateDataFromBlocksToNodes(CRFProcess* m_pcs,
                                                    std::string path,
                                                    int phase_index)
{
	(void)path; // unused
	clock_t start,finish;
	double time;
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	CECLIPSEBlock* m_block = NULL;
	MeshLib::CNode* m_node = NULL;
	//CFaces *m_face=NULL;
	//WW double distance;
	double volume;
	double weight;
	//double weights_xyz[3];
	double sum_weights, sum_weights_density;
	//WW double pressure,
	double saturation, sat, phase_pressure, phase_press, gas_dissolved, gas_dis;
	double den, density;
	//double press;
	//double velocity_x1, velo_x1, velocity_y1, velo_y1, velocity_z1, velo_z1;
	CPointData_ECL* m_NodeData = NULL;
	//m_NodeData = new CPointData;
	//WW long variable_index_pressure = -1, variable_index_phase_pressure = -1, variable_index_saturation = -1, time_index = 0, variable_index_Gas_dissolved = -1, variable_index_Gas_density = -1, variable_index_Water_density = -1, variable_index_Oil_density = -1;
	long variable_index_phase_pressure = -1, variable_index_saturation = -1,
	     variable_index_Gas_dissolved = -1, variable_index_Gas_density = -1,
	     variable_index_Water_density = -1, variable_index_Oil_density = -1;
	//long variable_index_velocity_x1, variable_index_velocity_y1, variable_index_velocity_z1;

	start = clock();

	cout << "        InterpolateDataFromBlocksToNodes()";

	//WW variable_index_pressure = this->GetVariableIndex("PRESSURE");

	//get saturation index if there are more than 1 phases
	if (this->Phases[phase_index] == "WATER")
	{
		if (int(this->Phases.size()) > 1)
			variable_index_saturation = this->GetVariableIndex("SWAT");
		variable_index_phase_pressure = this->GetVariableIndex("PWAT");
		if (this->E100 == true)
			variable_index_Water_density = this->GetVariableIndex("WAT_DEN");
		else
			variable_index_Water_density = this->GetVariableIndex("DENW");
	}
	else
	{
		if (this->Phases[phase_index] == "OIL")
		{
			if (int(this->Phases.size()) > 1)
				variable_index_saturation = this->GetVariableIndex("SOIL");
			variable_index_phase_pressure = this->GetVariableIndex("POIL");
			if (this->E100 == true)
				variable_index_Oil_density = this->GetVariableIndex("OIL_DEN");
			else
				variable_index_Oil_density = this->GetVariableIndex("DENO");
		}
		else
		{
			if (this->Phases[phase_index] == "GAS")
			{
				if (int(this->Phases.size()) > 1)
					variable_index_saturation = this->GetVariableIndex("SGAS");
				variable_index_phase_pressure = this->GetVariableIndex("PGAS");
				if (this->GetVariableIndex("RS") >= 0)
				{
					variable_index_Gas_dissolved = this->GetVariableIndex("RS");
					if (this->E100 == true)
						variable_index_Gas_density = this->GetVariableIndex(
						        "GAS_DEN");
					else
						variable_index_Gas_density = this->GetVariableIndex(
						        "DENG");
				}
			}
			else
			{
				cout << "This phase is not considered yet!" << endl;
				system("Pause");
				exit(0);
			}
		}
	}

	//WW time_index = 0;

	for (unsigned long i = 0; i < m_msh->nod_vector.size(); i++)
	{
		//Get the node
		m_node = m_msh->nod_vector[i];
		//WW	pressure = 0.0;
		phase_pressure = 0.0;
		saturation = 0.0;
		sum_weights = 0.0;
		gas_dissolved = 0.0;
		density = 0.0;
		sum_weights_density = 0.0;

		// sum the distance weighted data from each connected block
		for (unsigned int j = 0; j < m_node->getConnectedElementIDs().size(); j++)
		{
			//WW distance = 0.0;
			weight = 0.0;
			m_block  =  this->eclgrid[m_node->getConnectedElementIDs()[j]];

			//calculate representive volume of the considered node in each connected element for weighting
			volume = m_block->volume / 8.0; // ToDo volume weighting doesn't work if element is not simple!!!!!!!!! BG
			weight = (1.0 / volume);

			// Sum of weights
			sum_weights += weight;
			//Pressure
			//press = this->Data[m_block->index][variable_index_pressure];
			//pressure += press * weight;
			//Phase pressure
			phase_press = this->Data[m_block->index][variable_index_phase_pressure];
			phase_pressure += phase_press * weight;

			// saturation
			if (int(this->Phases.size()) > 1)
			{
				sat = this->Data[m_block->index][variable_index_saturation];
				saturation += sat * weight;
			}

			//dissolved gas
			if (variable_index_Gas_dissolved > -1)
			{
				gas_dis = this->Data[m_block->index][variable_index_Gas_dissolved];
				gas_dissolved += gas_dis * weight;
				den = this->Data[m_block->index][variable_index_Gas_density];
				if (den > 0)
				{
					density += den * weight;
					sum_weights_density += weight;
				}
			}

			if (variable_index_Water_density > -1)
			{
				den = this->Data[m_block->index][variable_index_Water_density];
				if (den > 0)
				{
					density += den * weight;
					sum_weights_density += weight;
				}
			}
			if (variable_index_Oil_density > -1)
			{
				den = this->Data[m_block->index][variable_index_Oil_density];
				if (den > 0)
				{
					density += den * weight;
					sum_weights_density += weight;
				}
			}

			// Test output
			//if ((m_block->z_barycentre == -2001.8)) {
			//	//cout << " Node " << i << endl;
			//	//cout << " :  connected block: " << m_node->connected_elements[j] << " (X,Y,Z): " << m_block->x_barycentre << ", " << m_block->y_barycentre << ", " << m_block->z_barycentre << " with saturation " << sat;
			//	//cout << "  distance: " << distance <<  " weight: " << weight << " sum of weights: " << sum_weights << endl;
			//	cout <<  "      Node " << i << " Element: " << m_node->connected_elements[j] << " SGAS: " << sat << endl;
			//}
		} // for j=connected_elements

		//pressure = pressure / sum_weights;
		//this->NodeData[i]->pressure = pressure;

		phase_pressure = phase_pressure / sum_weights;
		this->NodeData[i]->phase_pressure[phase_index] = phase_pressure;
		this->NodeData[i]->pressure = phase_pressure; // the last phase pressure is the highest pressure

		if (int(this->Phases.size()) > 1)
		{
			saturation = saturation / sum_weights;
			if ((saturation >= 0.0) && (saturation <= 1.0))
				this->NodeData[i]->phase_saturation[phase_index] = saturation;
			else
			{
				cout << "The saturation is not between 0 and 1!" << endl;
				if (saturation < 0.0)
					this->NodeData[i]->phase_saturation[phase_index] = 0.0;
				if (saturation > 1.0)
					this->NodeData[i]->phase_saturation[phase_index] = 1.0;
			}
		}
		//density
		if (sum_weights_density > 0)
			density = density / sum_weights_density;
		this->NodeData[i]->phase_density[phase_index] = density;

		//dissolved gas in oil (black oil mode)
		if (variable_index_Gas_dissolved > -1)
		{
			gas_dissolved = gas_dissolved / sum_weights;
			this->NodeData[i]->CO2inLiquid = gas_dissolved;
		}

		//Test output
		//if ((m_block->z_barycentre == -2001.8)) {
		//	//cout << " Node " << i << " (X,Y,Z): (" << m_node->X() << ", " << m_node->Y() << ", "<< m_node->X() << ") ";
		//	//cout << " Pressure at Node " << this->NodeData[i]->pressure << endl;
		//	//cout << " Node " << i << " (X,Y,Z): (" << m_node->X() << ", " << m_node->Y() << ", "<< m_node->X() << ") ";
		//	//cout << " Saturation at Node " << this->NodeData[i]->phase_saturation[phase_index] << endl;
		//	cout << " Node " << i << " SGAS: " <<  this->NodeData[i]->phase_saturation[phase_index] << endl;
		//}
	}

	// Test Output Nodes
	//if ((m_pcs->Tim->step_current==0) || (m_pcs->Tim->step_current==1) || (m_pcs->Tim->step_current==2) || (m_pcs->Tim->step_current==5) || (m_pcs->Tim->step_current==10) || (m_pcs->Tim->step_current==20) || (m_pcs->Tim->step_current==30) || (m_pcs->Tim->step_current==40) || (m_pcs->Tim->step_current==50) || (m_pcs->Tim->step_current==100) || (m_pcs->Tim->step_current==999) || (m_pcs->Tim->step_current==1000)) {
	vector <string> vec_string;
	std::string tempstring;
	ostringstream temp;
	double Val;
	if (this->Phases.size() == 1)
	{
		int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
		vec_string.push_back("Node; X; Y; Z; P-Geosys; P-Eclipse");
		for(unsigned long i = 0; i < this->NodeData.size(); i++)
		{
			m_NodeData = this->NodeData[i];
			//Get the Pressure from the Geosys run
			Val = m_pcs->GetNodeValue(i, nidx1);
			//cout << " NodeData[" << i << "]:  (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", " << m_NodeData->z << "):        " ;
			//cout << ", Geosys-P: " << Val << ", Eclipse: " << this->NodeData[i]->pressure << endl;;
			temp.str("");
			temp.clear();
			temp << i;
			tempstring = temp.str();
			temp.str("");
			temp.clear();
			temp << m_NodeData->x;
			tempstring += "; " + temp.str();
			temp.str("");
			temp.clear();
			temp << m_NodeData->y;
			tempstring += "; " + temp.str();
			temp.str("");
			temp.clear();
			temp << m_NodeData->z;
			tempstring += "; " + temp.str();
			temp.str("");
			temp.clear();
			temp.precision(12);
			temp << Val;
			tempstring += "; " + temp.str();
			temp.str("");
			temp.clear();
			temp.precision(12);
			temp << this->NodeData[i]->pressure;
			tempstring += "; " + temp.str();
			vec_string.push_back(tempstring);
		}
	}
	if (this->Phases.size() > 1)
		if (phase_index > 0)
		{
			int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
			int nidx2 = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;
			int nidx3 = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
			//vec_string.push_back("Node; X; Y; Z; P-Geosys_phase1; P-Geosys_phase2; Saturation_Geosys_phase1; P-Eclipse_phase1; P-Eclipse_phase2; Saturation_Eclipse_phase1");
			vec_string.push_back(
			        "Node; X; Y; Z; P-Geosys_phase1; P-Geosys_phase2; Saturation_Geosys_phase1; P-Eclipse_phase1; P-Eclipse_phase2; Saturation_Eclipse_phase1; q_x1; q_y1; q_z1; q_x2; q_y2; q_z2");
			for(unsigned long i = 0; i < this->NodeData.size(); i++)
			{
				m_NodeData = this->NodeData[i];
				//cout << " NodeData[" << i << "]:  (X,Y,Z): (" << m_NodeData->x << ", " << m_NodeData->y << ", " << m_NodeData->z << "):        " ;
				//cout << ", Geosys-P: " << Val << ", Eclipse: " << this->NodeData[i]->pressure << endl;;
				temp.str("");
				temp.clear();
				temp << i;
				tempstring = temp.str();
				temp.str("");
				temp.clear();
				temp << m_NodeData->x;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp << m_NodeData->y;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp << m_NodeData->z;
				tempstring += "; " + temp.str();

				//Get the Pressure1 from the Geosys run
				Val = m_pcs->GetNodeValue(i, nidx1);
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();
				//Get the Pressure2 from the Geosys run
				Val = m_pcs->GetNodeValue(i, nidx2);
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();
				//Get the Saturation1 from the Geosys run
				Val = m_pcs->GetNodeValue(i, nidx3);
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();

				//Eclipse data
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << this->NodeData[i]->phase_pressure[0];
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << this->NodeData[i]->phase_pressure[1];
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << this->NodeData[i]->phase_saturation[0];
				tempstring += "; " + temp.str();

				////Velocities
				//temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[0][0]; tempstring += "; " + temp.str();
				//temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[0][1]; tempstring += "; " + temp.str();
				//temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[0][2]; tempstring += "; " + temp.str();
				//temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[1][0]; tempstring += "; " + temp.str();
				//temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[1][1]; tempstring += "; " + temp.str();
				//temp.str(""); temp.clear(); temp.precision(12); temp << this->NodeData[i]->q[1][2]; tempstring += "; " + temp.str();

				vec_string.push_back(tempstring);
			}
		}
	//CWriteTextfiles *TextFile;
	//TextFile = new CWriteTextfiles;
	//int position = int(path.find_last_of("\\"));
	//std::string path_new;
	//path_new = path.substr(0,position);
	//position = int(path_new.find_last_of("\\"));
	//path_new = path_new.substr(0,position);
	//temp.str(""); temp.clear(); temp << m_pcs->Tim->step_current; tempstring = temp.str();
	//TextFile->Write_Text(path_new + "\\CheckNodeValues_" + tempstring + ".csv", vec_string);
	//} // End Test Output

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "    Time: " << time << " seconds." << endl;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: InterpolateGeosysVelocitiesToNodes
   Task: Returns the Nodevelocity of Geosys; v at the node as a inverse distance
   weighted mean of the connecting elements velocities
   Return: Velocity
   Programming: 11/2009 BG based on CB
   Modification:
   09/2011	TF changed access to coordinates of mesh node,
    - substituted access to mesh_element from pointer to direct access into the vector
    - made the mesh node a const pointer
    - made the pointer to the mesh const, made the mesh itself const
    - substituted pow(x,2) by x*x
   -------------------------------------------------------------------------*/
void CECLIPSEData::InterpolateGeosysVelocitiesToNodes(CRFProcess* m_pcs,
                                                      double* vel_nod, long node_number)
{
	CFEMesh const* const mesh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	size_t elem;
	double vel_ele[3];
	double distance, weight, sum_w (0.0);
	double PoreVel(0.0); //WW , theta;

	//double *vel_nod;
	//WW theta = m_pcs->m_num->ls_theta;

	// initialize data structures
	for (size_t i = 0; i < 3; i++)
		vel_nod[i] = vel_ele[i] = 0;

	// Get node coordinates
	MeshLib::CNode const* node(mesh->nod_vector[node_number]);
	double const* const coord(node->getData()); //Coordinates(coord);
	// get the indices of velocity of flow process
	long idxVx = m_pcs->GetElementValueIndex("VELOCITY1_X") + 1;
	long idxVy = m_pcs->GetElementValueIndex("VELOCITY1_Y") + 1;
	long idxVz = m_pcs->GetElementValueIndex("VELOCITY1_Z") + 1;

	for (size_t el = 0; el < node->getConnectedElementIDs().size(); el++)
	{
		distance = weight = 0; // initialize for each connected element
		elem = node->getConnectedElementIDs()[el];

		// get the velocity components of element elem
		// if idxVx = 0 it implies that the index for this parameter doesn't exist
		if (idxVx != 0)
			vel_ele[0] = m_pcs->GetElementValue(elem, idxVx);
		if (idxVy != 0)
			vel_ele[1] = m_pcs->GetElementValue(elem, idxVy);
		if (idxVz != 0)
			vel_ele[2] = m_pcs->GetElementValue(elem, idxVz);
		// calculate distance node <-> element center of gravity
		double const* grav_c(mesh->ele_vector[elem]->GetGravityCenter());
		for (size_t i = 0; i < 3; i++)
			distance += (coord[i] - grav_c[i]) * (coord[i] - grav_c[i]);  // TF pow((coord[i]-grav_c[i]),2);
		// linear inverse distance weight = 1/(distance)
		distance = sqrt(distance); // for quadratic interpolation uncomment this line
		weight = (1 / distance);
		sum_w += weight;
		for (size_t i = 0; i < 3; i++)
			// the 3 velocity components
			vel_nod[i] += vel_ele[i] * weight;
	}
	// normalize weighted sum by sum_of_weights sum_w
	const double sum_w_inv(1.0 / sum_w);
	for (size_t i = 0; i < 3; i++)
		vel_nod[i] *= sum_w_inv;
	// absolute value of velocity vector
	for (size_t i = 0; i < 3; i++)
		PoreVel += vel_nod[i] * vel_nod[i];  // TF pow(vel_nod[i],2);
	PoreVel = sqrt(PoreVel);
}

/*-------------------------------------------------------------------------
   GeoSys - Function: WriteDataToGeoSys
   Task: Uses the node velocities to calculate the velocity at gauss points
   Return: true or false
   Programming: 09/2009 SB
   Modification:
   -------------------------------------------------------------------------*/
void CECLIPSEData::WriteDataToGeoSys(CRFProcess* m_pcs, std::string path)
{
	//Daten�bergabe: pressure, saturation (ev. kf, n)
	long index_pressure1, index_pressure2, index_saturation1, index_water_density,
	     index_gas_density;
	//long index_saturation2, index_gas_density;
	int phase1, phase2;
	double value = 0.0;
	//double Molweight_CO2;

	// Pressure 1 is the wetting phase, Pressure 2 the non wetting phase
	//Pressure
	if (int(this->Phases.size()) == 1)
	{
		index_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
		index_water_density = m_pcs->GetNodeValueIndex("DENSITY1");
		for(unsigned long i = 0; i < m_pcs->nod_val_vector.size(); i++)
		{
			value = this->NodeData[i]->pressure;
			m_pcs->SetNodeValue(i, index_pressure1, value);
			value = this->NodeData[i]->phase_density[0];
			m_pcs->SetNodeValue(i, index_water_density, value);
		}
	}
	else if (int(this->Phases.size()) > 1)
	{
		if ((this->Phases[0] == "WATER") && (this->Phases[1] == "GAS") &&
		    (this->E100 == true))
		{
			cout <<
			"GAS-WATER System can not be considered with ECLIPSE E100 and GeoSys" <<
			endl;
			system("Pause");
			exit(0);
		}
		switch(int(this->Phases.size()))
		{
		case 2:
			phase1 = 0;
			phase2 = 1;
			break;
		case 3:
			// Assumption that there are 3 phases but water is only used for the boundaries -> the oil and gas phase are the relevant one for the exchange with OGS
			phase1 = 1;
			phase2 = 2;
			break;
		default:
			cout << "There are not more than 3 phases possible!" << endl;
			system("Pause");
			exit(0);
			break;
		}
		index_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
		index_pressure2 = m_pcs->GetNodeValueIndex("PRESSURE2") + 1; //+1... new time level
		index_saturation1 = m_pcs->GetNodeValueIndex("SATURATION1") + 1; //+1... new time level
		index_water_density = m_pcs->GetNodeValueIndex("DENSITY1");
		index_gas_density = m_pcs->GetNodeValueIndex("DENSITY2");

		for(unsigned long i = 0; i < m_pcs->nod_val_vector.size(); i++)
		{
			//if ((this->NodeData[i]->x == 499.85001) && (this->NodeData[i]->y == 599.84998) && (this->NodeData[i]->z == -2875))
			//	cout << "Punkt 44" << endl;
			//calculating capillary pressure (difference of the two phase pressures)
			//cout << " Node " << i << " (X,Y,Z): (" << this->NodeData[i]->x << ", " << this->NodeData[i]->y << ", "<< this->NodeData[i]->z << ") ";
			//for (long k = 0; k < 3; k++){
			//	cout << this->NodeData[i]->phase_density[k] << " " << this->NodeData[i]->phase_pressure[k] << " " << this->NodeData[i]->phase_saturation[k] << " ";
			//}
			//cout << endl;

			value = this->NodeData[i]->phase_pressure[phase2] -
			        this->NodeData[i]->phase_pressure[phase1];
			m_pcs->SetNodeValue(i, index_pressure1, value);
			value = this->NodeData[i]->phase_pressure[phase2];
			m_pcs->SetNodeValue(i, index_pressure2, value);
			value = this->NodeData[i]->phase_density[phase1];
			m_pcs->SetNodeValue(i, index_water_density, value);
			value = this->NodeData[i]->phase_density[phase2];
			m_pcs->SetNodeValue(i, index_gas_density, value);

			if (this->Phases.size() < 3)
			{
				//Saturation 1 is the wetting phase
				value = this->NodeData[i]->phase_saturation[phase1];
				//cout << " Node: " << i << " saturation1: " << value << endl;
				m_pcs->SetNodeValue(i, index_saturation1, value);
			}
			else
			{
				// the water phase is only used for boundary conditions, the wetting phase is phase 2 and the non-wetting phase is phase 3
				value = 1 - this->NodeData[i]->phase_saturation[phase2];
				//cout << " Node: " << i << " saturation1: " << value << endl;
				m_pcs->SetNodeValue(i, index_saturation1, value);
			}
			//transfer the gas density
			//if (this->GetVariableIndex("RS") >= 0) {
			//	value = this->NodeData[i]->phase_density[phase2];
			//	m_pcs->SetNodeValue(i, index_gas_density, value);
			//}
		}
		//assign dissolved gas to GeoSys nodes
		if (this->GetVariableIndex("RS") >= 0)
		{
			CRFProcess* n_pcs = NULL;
			int indexConcentration;
			if (this->ProcessIndex_CO2inLiquid == -1)
				for(int i = 0; i < int(pcs_vector.size()); i++)
				{
					n_pcs = pcs_vector[i];
					// identify your process and store idx of pcs-vector
					if ((n_pcs->nod_val_name_vector[0] == "C(4)") ||
					    (n_pcs->nod_val_name_vector[0] == "CO2_w"))          // "C(4)"...Phreeqc, "CO2_w"...Chemapp
						this->ProcessIndex_CO2inLiquid = i;
				}
			if (this->ProcessIndex_CO2inLiquid == -1)
			{
				cout <<
				"In the model exists dissolved gas but there is no dissolved C in water defined. Please ad the mass transport for immobile C(4) or CO2_w!"
				     << endl;
				system("Pause");
				exit(0);
			}

			// get index of species concentration in nodevaluevector of this process
			indexConcentration =
			        pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
			                pcs_vector[
			                        this
			                        ->ProcessIndex_CO2inLiquid]->
			                pcs_primary_function_name[0]) + 1;                                                                                        // +1: new timelevel
			for(unsigned long i = 0;
			    i < pcs_vector[this->ProcessIndex_CO2inLiquid]->nod_val_vector.size();
			    i++)
			{
				//recalculate dissolve gas: c_CO2 [mol/m�] = RS [m�gas/m�water] * gas_density[kg/m�] / (Molweight_CO2 [g/mol] * 1e-3 [kg/g])
				//value = this->NodeData[i]->Gas_dissolved * this->NodeData[i]->phase_density[phase2] / Molweight_CO2;
				value = this->NodeData[i]->CO2inLiquid * this->SurfaceCO2Density /
				        (Molweight_CO2 * 1e-3);
				pcs_vector[this->ProcessIndex_CO2inLiquid]->SetNodeValue(
				        i,
				        indexConcentration,
				        value);
				//cout << "Node: " << i << " Druck: " << m_pcs->GetNodeValue(i, index_pressure2) << " RS: " << this->NodeData[i]->Gas_dissolved << " Dichte: " << this->SurfaceCO2Density << " C(4): " << value << endl;
				if (value < 0)
				{
					cout << "Node: " << i << " Druck: " << m_pcs->GetNodeValue(
					        i,
					        index_pressure2) << " RS: " <<
					this->NodeData[i]->CO2inLiquid << " Dichte: " <<
					this->SurfaceCO2Density << " C(4): " << value << endl;
					cout << "  Fehler in Berechnung von DIC: " << value << endl;
				}
			}
		}
	}

	//Geschwindigkeiten an den Gau�punkten interpolieren
	for (int i = 0; i < int(this->Phases.size()); i++)
		m_pcs->CalGPVelocitiesfromECLIPSE(path,
		                                  m_pcs->Tim->step_current,
		                                  i,
		                                  this->Phases[i]);

	// Test Output Elements
	CFEMesh* m_msh = fem_msh_vector[0]; //SB: ToDo hart gesetzt
	CECLIPSEBlock* m_block = NULL;

	if ((m_pcs->Tim->step_current == 0) || (m_pcs->Tim->step_current == 1) ||
	    (m_pcs->Tim->step_current == 2) || (m_pcs->Tim->step_current == 5) ||
	    (m_pcs->Tim->step_current == 10) || (m_pcs->Tim->step_current == 20) ||
	    (m_pcs->Tim->step_current == 30) || (m_pcs->Tim->step_current == 40) ||
	    (m_pcs->Tim->step_current == 50) || (m_pcs->Tim->step_current == 100) ||
	    (m_pcs->Tim->step_current == 999) || (m_pcs->Tim->step_current == 1000))
	{
		vector <string> vec_string;
		std::string tempstring;
		ostringstream temp;
		MeshLib::CElem* elem = NULL;
		Math_Group::vec<MeshLib::CNode*> Nodes(8);
		double Val;
		if (this->Phases.size() == 1)
		{
			int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
			vec_string.push_back("Element; X; Y; Z; P-Geosys");
			for(unsigned long i = 0; i <  m_msh->ele_vector.size(); i++)
			{
				elem = m_msh->ele_vector[i];
				m_block = this->eclgrid[i];
				//Calculate the average of the pressure for all nodes of 1 element
				Val = 0;
				elem->GetNodes(Nodes);
				for (long j = 0; j < long(Nodes.Size()); j++)
					Val +=
					        m_pcs->GetNodeValue(Nodes[j]->GetIndex(),
					                            nidx1) / Nodes.Size();

				temp.str("");
				temp.clear();
				temp << i;
				tempstring = temp.str();
				temp.str("");
				temp.clear();
				temp << m_block->x_barycentre;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp << m_block->y_barycentre;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp << m_block->z_barycentre;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();
				vec_string.push_back(tempstring);
			}
		}
		if (this->Phases.size() > 1)
		{
			int nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1; //+1... new time level
			int nidx2 = m_pcs->GetNodeValueIndex("PRESSURE2") + 1;
			int nidx3 = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
			vec_string.push_back(
			        "Element; X; Y; Z; P-Geosys_phase1; P-Geosys_phase2; Saturation_Geosys_phase1");
			for(unsigned long i = 0; i <  m_msh->ele_vector.size(); i++)
			{
				elem = m_msh->ele_vector[i];
				m_block = this->eclgrid[i];
				temp.str("");
				temp.clear();
				temp << i;
				tempstring = temp.str();
				temp.str("");
				temp.clear();
				temp << m_block->x_barycentre;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp << m_block->y_barycentre;
				tempstring += "; " + temp.str();
				temp.str("");
				temp.clear();
				temp << m_block->z_barycentre;
				tempstring += "; " + temp.str();
				Val = 0;
				elem->GetNodes(Nodes);
				for (long j = 0; j < long(Nodes.Size()); j++)
					Val +=
					        m_pcs->GetNodeValue(Nodes[j]->GetIndex(),
					                            nidx1) / Nodes.Size();
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();

				Val = 0;
				for (long j = 0; j < long(Nodes.Size()); j++)
					Val +=
					        m_pcs->GetNodeValue(Nodes[j]->GetIndex(),
					                            nidx2) / Nodes.Size();
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();

				Val = 0;
				for (long j = 0; j < long(Nodes.Size()); j++)
					Val +=
					        m_pcs->GetNodeValue(Nodes[j]->GetIndex(),
					                            nidx3) / Nodes.Size();
				temp.str("");
				temp.clear();
				temp.precision(12);
				temp << Val;
				tempstring += "; " + temp.str();

				vec_string.push_back(tempstring);
			}
		}
		//CWriteTextfiles *TextFile;
		//TextFile = new CWriteTextfiles;
		//int position = int(path.find_last_of("\\"));
		//std::string path_new;
		//path_new = path.substr(0,position);
		//position = int(path_new.find_last_of("\\"));
		//path_new = path_new.substr(0,position);
		//temp.str(""); temp.clear(); temp << m_pcs->Tim->step_current; tempstring = temp.str();
		//TextFile->Write_Text(path_new + "\\ElementValues_Geosys_" + tempstring + ".csv", vec_string);
	} // End Test Output
}
/*-------------------------------------------------------------------------
   GeoSys - Function: ExecuteEclipse
   Task: Starts Eclipse
   Return: Projectname (string)
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
std::string CECLIPSEData::ExecuteEclipse(long Timestep, CRFProcess* m_pcs, std::string folder)
{
	std::string tempstring;
	std::string EclipseExe; // = "c:\\programme\\ecl\\2008.1\\bin\\pc\\eclipse.exe";
	std::string projectname;
	std::string Filename;
	int position;
	std::string Keyword;
	std::string Keyword_well;
	ostringstream temp;
	vector <string> vecData;
	CReadTextfiles_ECL* TextFile;

	bool Error;
	bool Water_phase_exists = false, Oil_phase_exists = false, Gas_phase_exists = false;
	double timestep_length;
	size_t found_H2O, found_CO2, found_NaCl;
	std::string root_folder, geosys_folder;

	EclipseExe = m_pcs->simulator_path; // path to eclipse executable
	//check if file exists
	if (this->UsePrecalculatedFiles == false)
		if (CheckIfFileExists(EclipseExe) == false)
		{
			cout << "The ECLIPSE executable could not be found! (" << EclipseExe <<
			")" << endl;
			system("Pause");
			exit(0);
		}

	Filename = m_pcs->simulator_model_path; // path to eclipse input data
	//check if file exists
	if (CheckIfFileExists(Filename) == false)
	{
		cout << "The ECLIPSE project file could not be found! (" << Filename << ")" << endl;
		system("Pause");
		exit(0);
	}

	//check if eclipse folder is subfolder of the geosys folder
	if (Timestep == 1)
	{
		if (this->Windows_System == true)
			position = int(Filename.find_last_of("\\"));
		else
			position = int(Filename.find_last_of("/"));
		root_folder = Filename.substr(0, position);
		if (this->Windows_System == true)
			position = int(root_folder.find_last_of("\\"));
		else
			position = int(root_folder.find_last_of("/"));
		root_folder = Filename.substr(0, position);

		if (this->Windows_System == true)
			position = int(m_pcs->file_name_base.find_last_of("\\"));
		else
			position = int(m_pcs->file_name_base.find_last_of("/"));
		geosys_folder = m_pcs->file_name_base.substr(0, position);

		if (root_folder != geosys_folder)
		{
			cout <<
			"Warning: The Eclipse simulator model path is not part of the GeoSys model path!!!"
			     << endl;
			cout << root_folder << endl;
			cout << geosys_folder << endl;
			system("Pause");
		}
	}

	if (this->Windows_System == true)
		position = int(Filename.find_last_of("\\"));
	else
		position = int(Filename.find_last_of("/"));
	// check if filename is given with or without extension
	tempstring = Filename.substr(Filename.length() - 5,Filename.length());
	if((tempstring.compare(".data") == 0) || (tempstring.compare(".DATA") == 0))
		projectname = Filename.substr(0, Filename.length() - 5);
	else
		projectname = Filename;
	if (Timestep > 1)
		projectname = projectname + "_RESTART";

	//Read length of current timestep and recalculate it to days
	if (m_pcs->Tim->time_unit == "DAY")
		timestep_length = m_pcs->Tim->time_step_length;
	else
	{
		if (m_pcs->Tim->time_unit == "MINUTE")
			timestep_length = m_pcs->Tim->time_step_length / 60 / 24;
		else
		{
			if (m_pcs->Tim->time_unit == "SECOND")
				timestep_length = m_pcs->Tim->time_step_length / 60 / 60 / 24;
			else
			{
				cout << "This time unit was not considered yet" << endl;
				system("Pause");
				exit(0);
			}
		}
	}
	//Write the timestep in the eclipse file
	if (this->UsePrecalculatedFiles == false)
	{
		Keyword = "TSTEP";
		tempstring = "1*";
		temp.str("");
		temp.clear();
		temp << timestep_length;
		tempstring = tempstring + temp.str();
		tempstring = tempstring + " /";
		vecData.clear();
		vecData.push_back(tempstring);
		std::string File = projectname + ".data";
		if (ReplaceASectionInFile(File, Keyword, vecData, true) == false)
		{
			cout << "Replacing a section in the file: " << File <<
			" didn't work for Keyword TSTEP!" << endl;
			system("Pause");
			exit(0);
		}
		vecData.clear();

		// KB: replace injection/production rate
		if (this->ecl_well.size() != 0)
		{
			for (int ii = 0; ii < int(this->ecl_well.size()); ii++)
			{
				for (int kk = 0; kk < int(this->ecl_well[ii]->time.size()); kk++)

					if (this->ecl_well[ii]->time[kk] == this->actual_time)
					{
						this->WellRates.push_back(ecl_well[ii]->rate[kk]);

						break;
					}
				std::string outline = this->ecl_well[ii]->name + "   " +
				                      this->ecl_well[ii]->phase + "   " +
				                      this->ecl_well[ii]->open_flag +
				                      "   " + this->ecl_well[ii]->control_mode +
				                      "   " +
				                      this->WellRates[ii] + "   " + " / ";
				vecData.push_back(outline);
			}
			Keyword = "WCONINJE";
			std::string File = projectname + ".data";
			if (ReplaceASectionInFile(File, Keyword, vecData, true) == false)
			{
				cout << "Replacing a section in the file: " << File <<
				" didn't work!" << endl;
				system("Pause");
				exit(0);
			}
		}
		vecData.clear();
	}

	if (Timestep == 1)
	{
		//only if no phases exists yet
		this->E100 = false;
		if (int(this->Phases.size()) == 0)
		{
			// Read used phases from the Eclipse input file
			TextFile = new CReadTextfiles_ECL;
			Error = TextFile->Read_Text(Filename);
			if (Error == true)
			{
				cout << "The program is canceled" << endl;
				system("Pause");
				exit(0);
			}
			for (int i = 0; i < TextFile->NumberOfRows; i++)
			{
				this->SplittedString.clear();
				this->SplitStrings(TextFile->Data[i]," ");
				if (this->SplittedString.size() > 0)
				{
					// check if one of the phases is water
					if (this->SplittedString[0] == "WATER")
					{
						Water_phase_exists = true;
						this->E100 = true;
					}
					if (this->SplittedString[0] == "GAS")
					{
						Gas_phase_exists = true;
						this->E100 = true;
					}
					if (this->SplittedString[0] == "OIL")
					{
						Oil_phase_exists = true;
						this->E100 = true;
					}
					if (this->SplittedString[0] == "CNAMES")
					{
						this->E100 = false;
						found_H2O = TextFile->Data[i + 1].find("H2O");
						found_CO2 = TextFile->Data[i + 2].find("CO2");
						found_NaCl = TextFile->Data[i + 3].find("NACL");
						if ((found_H2O == string::npos) ||
						    (found_CO2 == string::npos) ||
						    (found_NaCl == string::npos))
						{
							cout <<
							"The OGS-ECL code can currently only consider the components H2O, CO2 and NaCl in this order in E300!"
							     << endl;
							system("Pause");
							exit(0);
						}
						i = TextFile->NumberOfRows;
					}
				}
			}
			//fill the phase vector with the wetting phase first
			if (this->E100 == true)
			{
				if (Water_phase_exists == true)
					this->Phases.push_back("WATER");
				if (Oil_phase_exists == true)
					this->Phases.push_back("OIL");
				if (Gas_phase_exists == true)
					this->Phases.push_back("GAS");
				if (Water_phase_exists == false)
				{
					// The software is not checked if water is not the wetting phase
					cout <<
					"The program is canceled because no water phase is defined in ECLIPSE!"
					     << endl;
					system("Pause");
					exit(0);
				}

				//Define necessary Variables
				//The total pressure is only defined once
				this->Variables.push_back("PRESSURE");
				this->Variables.push_back("ACAQNUM");
				for (unsigned int i = 0; i < this->Phases.size(); i++)
				{
					if (this->Phases[i].compare("WATER") == 0)
					{
						this->Variables.push_back("FLOWATI+");
						this->Variables.push_back("FLOWATJ+");
						this->Variables.push_back("FLOWATK+");
						this->Variables.push_back("PWAT");
						this->Variables.push_back("WAT_DEN");
						if (this->Phases.size() > 1)
							this->Variables.push_back("SWAT");
					}
					if (this->Phases[i].compare("OIL") == 0)
					{
						this->Variables.push_back("FLOOILI+");
						this->Variables.push_back("FLOOILJ+");
						this->Variables.push_back("FLOOILK+");
						this->Variables.push_back("POIL");
						this->Variables.push_back("OIL_DEN");
						if (this->Phases.size() > 1)
							this->Variables.push_back("SOIL");
					}
					if (this->Phases[i].compare("GAS") == 0)
					{
						this->Variables.push_back("FLOGASI+");
						this->Variables.push_back("FLOGASJ+");
						this->Variables.push_back("FLOGASK+");
						this->Variables.push_back("PGAS");
						this->Variables.push_back("GAS_DEN");
						if (this->Phases.size() > 1)
							this->Variables.push_back("SGAS");
					}
				}
				// variable for reading capillary pressures in the case of Eclipse E100
				if (int(this->Phases.size()) == 2)
				{
					if ((this->Phases[0] == "WATER") &&
					    (this->Phases[1] == "OIL"))
						this->Variables.push_back("PCOW");
					else
					{
						if ((this->Phases[0] == "OIL") &&
						    (this->Phases[1] == "GAS"))
						{
							this->Variables.push_back("PCOG");
							this->Variables.push_back("RS"); //gas saturation in oil in the black oil mode
							this->Variables.push_back("GAS_DEN"); //gas density in the black oil mode
						}
						else
						{
							cout <<
							"Currently only Oil-Water or Oil-Gas is available because there is no output option in Eclipse E100 for the capillary pressure!"
							     << endl;
							system("Pause");
							exit(0);
						}
					}
				}
				if (int(this->Phases.size()) == 3)
				{
					this->Variables.push_back("PCOW");
					this->Variables.push_back("PCOG");
					this->Variables.push_back("RS"); //gas saturation in oil in the black oil mode
					//this->Variables.push_back("GAS_DEN"); //gas density in the black oil mode
				}
			}
			else
			{
				//E300 mode -> it is assumed that water and CO2 phase are considered
				//components
				this->Components.push_back("H2O");
				this->Components.push_back("CO2");
				this->Components.push_back("NaCl");
				//phases
				this->Phases.push_back("WATER");
				this->Phases.push_back("GAS");
				//variables
				this->Variables.push_back("PRESSURE"); // Total pressure
				this->Variables.push_back("PWAT"); // water phase pressure
				this->Variables.push_back("PGAS"); // gas phase pressure

				this->Variables.push_back("SWAT"); // water phase saturation
				this->Variables.push_back("SGAS"); // gas phase saturation

				this->Variables.push_back("RS"); // aqueous CO2 concentration

				this->Variables.push_back("RPORV"); // pore volume of the cell

				this->Variables.push_back("DENW"); // water phase density
				this->Variables.push_back("DENG"); // gas phase density

				this->Variables.push_back("XFW1"); // mass fraction of component 1 in liquid phase
				this->Variables.push_back("XFW2"); // mass fraction of component 2 in liquid phase
				this->Variables.push_back("XFW3"); // mass fraction of component 3 in liquid phase

				//this->Variables.push_back("XMF1");				// mole fraction of component 1 in liquid phase
				//this->Variables.push_back("XMF2");				// mole fraction of component 2 in liquid phase
				//this->Variables.push_back("XMF3");				// mole fraction of component 3 in liquid phase

				//this->Variables.push_back("YMF1");				// mole fraction of component 1 in gas phase
				//this->Variables.push_back("YMF2");				// mole fraction of component 2 in gas phase
				//this->Variables.push_back("YMF3");				// mole fraction of component 3 in gas phase

				this->Variables.push_back("MLSC1"); // total molar density for components(moles per reservoir volume)
				this->Variables.push_back("MLSC2"); // total molar density for components(moles per reservoir volume)
				this->Variables.push_back("MLSC3"); // total molar density for components(moles per reservoir volume)

				//this->Variables.push_back("FLOC1I+");			// inter-block component flow for component 1
				//this->Variables.push_back("FLOC1J+");			// inter-block component flow for component 1
				//this->Variables.push_back("FLOC1K+");			// inter-block component flow for component 1

				//this->Variables.push_back("FLOC2I+");			// inter-block component flow for component 2
				//this->Variables.push_back("FLOC2J+");			// inter-block component flow for component 2
				//this->Variables.push_back("FLOC2K+");			// inter-block component flow for component 2

				//this->Variables.push_back("FLOC3I+");			// inter-block component flow for component 3
				//this->Variables.push_back("FLOC3J+");			// inter-block component flow for component 3
				//this->Variables.push_back("FLOC3K+");			// inter-block component flow for component 3

				this->Variables.push_back("FLOWATI+");
				this->Variables.push_back("FLOWATJ+");
				this->Variables.push_back("FLOWATK+");

				this->Variables.push_back("FLOGASI+");
				this->Variables.push_back("FLOGASJ+");
				this->Variables.push_back("FLOGASK+");
			}

			this->numberOutputParameters = int(Variables.size());
			this->times = 1;
		}

		//Run Eclipse the first time from the original *.data file

		// make string for external program call to ECLIPSE
		if (this->UsePrecalculatedFiles == false)
		{
			if (this->Windows_System == true)
			{
				tempstring = EclipseExe + " " + projectname;
				if (system(tempstring.c_str()))
				{
					DisplayMsgLn("Warnung: Eclipse doesn't run properly!!! ");
					return 0;
				}
			}
			else
			{
				tempstring = EclipseExe + " " + projectname;
				if (system(tempstring.c_str()))
				{
					DisplayMsgLn("Warnung: Eclipse doesn't run properly!!! ");
					return 0;
				}
			}
		}
	}
	else
	{
		// change Restart file
		// Changes are included in the restart File, because there it is possible to define new porosities ore permeabilties
		// every time at least the number of the time step which is used for restart has to be changed

		//Run Restart with results from the previous step
		tempstring = EclipseExe + " " + projectname;

		// increase the number of time steps used for restart
		Keyword = "RESTART";
		//tempstring = " d:\\Projekte\\MoPa\Programme\\Eclipse\\BGfirstmodel\\";
		tempstring = folder + "TemporaryResults ";
		//tempstring +=  "TemporaryResults ";
		temp.str("");
		temp.clear();
		temp << Timestep - 1;
		tempstring = tempstring + temp.str();
		tempstring = tempstring + " /";
		//tempstring = tempstring + " SAVE FORMATTED /";
		vecData.clear();
		vecData.push_back(tempstring);
		std::string File = projectname + ".data";
		if (this->UsePrecalculatedFiles == false)
			if (ReplaceASectionInFile(File, Keyword, vecData, true) == false)
			{
				cout << "Replacing a section in the file: " << File <<
				" didn't work for Keyword RESTART!" << endl;
				system("Pause");
				exit(0);
			}

		// Call ECLIPSE
		if (this->UsePrecalculatedFiles == false)
		{
			if (this->Windows_System == true)
			{
				tempstring = EclipseExe + " " + projectname;
				if (system(tempstring.c_str()))
				{
					DisplayMsgLn("Warnung: Eclipse doesn't run properly!!! ");
					return 0;
				}
			}
			else
			{
				tempstring = EclipseExe + " " + projectname;
				if (system(tempstring.c_str()))
				{
					DisplayMsgLn("Warnung: Eclipse doesn't run properly!!! ");
					return 0;
				}
			}
		}
	}
	return projectname;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: CleanUpEclipseFiles
   Task: Rename, Copy and Delete Files from the Eclipse Run
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
bool CECLIPSEData::CleanUpEclipseFiles(std::string folder, std::string projectname)
{
	std::string DOScommand;

	// Rename files for ECLIPSE restart
	//delete temporary result files

	DOScommand = "del " + folder + "TemporaryResult*.*";
	cout << DOScommand << endl;
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete temporary result files! ");
		return 0;
	}

	DOScommand = "del " + projectname + ".FSAVE /Q";
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete result files! ");
		return 0;
	}

	//copy original result files
	//DOScommand = "copy " + projectname + ".F* " + folder + "TemporaryResults.F*";
	DOScommand = "ren " + projectname + ".F* " + "TemporaryResults.F*";
	cout << DOScommand << endl;
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not copy the result files! ");
		return 0;
	}

	//DOScommand = "del " + projectname + ".F* /Q";
	//if (system(DOScommand.c_str())){
	//	DisplayMsgLn("Could not delete result files! ");
	//	return 0;
	//}
	DOScommand = "del " + projectname + ".A* /Q";
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete result files! ");
		return 0;
	}
	DOScommand = "del " + projectname + ".P* /Q";
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete result files! ");
		return 0;
	}
	DOScommand = "del " + projectname + ".R* /Q";
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete result files! ");
		return 0;
	}
	DOScommand = "del " + projectname + ".M* /Q";
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete result files! ");
		return 0;
	}
	DOScommand = "del " + projectname + ".DB* /Q";
	if (system(DOScommand.c_str()))
	{
		DisplayMsgLn("Could not delete result files! ");
		return 0;
	}

	return true;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: CleanUpEclipseFiles
   Task: Rename, Copy and Delete Files from the Eclipse Run
   Return: nothing
   Programming: 09/2009 BG
   Modification:
   -------------------------------------------------------------------------*/
int CECLIPSEData::WriteDataBackToEclipse(CRFProcess* m_pcs, std::string folder)
{
	std::string Filename;
	//WW int phase1, phase2;
	CRFProcess* n_pcs = NULL;
	//int indexProcess;
	int indexConcentration;
	MeshLib::CElem* m_element = NULL;
	CFEMesh* m_msh = fem_msh_vector[0];
	Math_Group::vec <MeshLib::CNode*> vec_element_nodes(8);
	clock_t start, finish;
	double time;
	double delta_gas_dis, delta_gas_dissolved;
	double weight;
	vector <string> vecString;
	vector <double> vecRS;
	vector <double> vec_T_CO2;
	double deltaT_CO2;
	std::string Keyword;
	ostringstream temp;
	std::string tempstring;
	int variable_index_Gas_dissolved, variable_index_porevolume,
	    variable_index_water_saturation, variable_index_Moles_CO2;
	typeExponentialNumber tempNumber;
	const double epsilon = 1e-7;
	int j_max;

	cout << "      Write DIC to Eclipse File " << endl;

	start = clock();

	switch(int(this->Phases.size()))
	{
	case 2:
		//WW	phase1 = 0;
		//WW	phase2 = 1;
		break;
	case 3:
		// Assumption that there are 3 phases but water is only used for the boundaries -> the oil and gas phase are the relevant one for the exchange with OGS
		//WW	phase1 = 1;
		//WW	phase2 = 2;
		break;
	default:
		cout << "There are not more than 3 phases possible!" << endl;
		system("Pause");
		return 0;
		break;
	}

	//clean node vector
	//this->MakeNodeVector(m_pcs);

	// get DIC and CO2 density at nodes and recalculate RS
	//index_gas_density = m_pcs->GetNodeValueIndex("DENSITY2");
	variable_index_Gas_dissolved = this->GetVariableIndex("RS");
	variable_index_porevolume = this->GetVariableIndex("RPORV");
	variable_index_water_saturation = this->GetVariableIndex("SWAT");
	variable_index_Moles_CO2 = this->GetVariableIndex("MLSC2");
	if (variable_index_Gas_dissolved < 0)
	{
		cout <<
		"Not all variables are existing in the eclipse output files that are necessary to write data back to e100!"
		     << endl;
		return 0;
	}
	if (((variable_index_porevolume < 0) || (variable_index_water_saturation < 0) ||
	     (variable_index_Moles_CO2 < 0)) && (this->E100 == false))
	{
		cout <<
		"Not all variables are existing in the eclipse output files that are necessary to write data back to e300!"
		     << endl;
		return 0;
	}

	// get index of species concentration in nodevaluevector of this process
	indexConcentration = pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
	        pcs_vector[this->ProcessIndex_CO2inLiquid]->pcs_primary_function_name[0]) + 1;                                                                    // +1: new timelevel
	for(unsigned long i = 0;
	    i < pcs_vector[this->ProcessIndex_CO2inLiquid]->nod_val_vector.size(); i++)
	{
		//Read CO2 density
		// check if concentration are positive
		if (pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(i,
		                                                             indexConcentration) <
		    0)
		{
			cout << "The CO2 concentration after the reactions is negativ!" << endl;
			cout << "Node: " << i << " RS_alt: " << this->NodeData[i]->CO2inLiquid <<
			" DIC: " << n_pcs->GetNodeValue(i, indexConcentration);
			system("Pause");
			return 0;
		}
		if (this->NodeData[i]->CO2inLiquid < 0)
		{
			cout << "The CO2 concentration from Eclipse is negativ!" << endl;
			cout << "Node: " << i << " RS_alt: " << this->NodeData[i]->CO2inLiquid <<
			" DIC: " << n_pcs->GetNodeValue(i, indexConcentration);
			system("Pause");
			return 0;
		}

		//this->NodeData[i]->phase_density[phase2] = m_pcs->GetNodeValue(i, index_gas_density);
		//recalculate dissolve gas: RS [m�gas/m�water] = c_CO2 [mol/m�water] * (Molweight_CO2 [g/mol] * 1e-3 [kg/g]) / gas_density[kg/m�gas]
		this->NodeData[i]->deltaDIC =
		        (pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(i,
		                                                                  indexConcentration)
		         *
		         (Molweight_CO2 *
		          1e-3) / this->SurfaceCO2Density) - this->NodeData[i]->CO2inLiquid;
		//cout << "Node: " << i << " RS_alt: " << this->NodeData[i]->Gas_dissolved << " DIC: " << n_pcs->GetNodeValue(i, indexConcentration);
		if (this->NodeData[i]->deltaDIC != 0)
		{
			//neglect difference if it is too small
			if (((fabs(this->NodeData[i]->deltaDIC) /
			      this->NodeData[i]->CO2inLiquid)) < epsilon)
				//cout << " deltaDIC_alt: " << this->NodeData[i]->deltaDIC << " Fehler! " << endl;;
				this->NodeData[i]->deltaDIC = 0;
			else
				this->NodeData[i]->CO2inLiquid =
				        (pcs_vector[this->ProcessIndex_CO2inLiquid]->GetNodeValue(i,
				                                                                  indexConcentration)
				         *
				         (Molweight_CO2 * 1e-3) / this->SurfaceCO2Density);
		}

		//cout << " RS_neu: " << this->NodeData[i]->Gas_dissolved << " deltaDIC: " << this->NodeData[i]->deltaDIC << endl;
		if (this->NodeData[i]->CO2inLiquid < 0)
		{
			cout << " RS_neu: " << this->NodeData[i]->CO2inLiquid << " deltaDIC: " <<
			this->NodeData[i]->deltaDIC << endl;
			cout << "  Fehler in Berechnung von DIC: " <<
			this->NodeData[i]->CO2inLiquid << endl;
		}
	}

	//interpolate change of RS from nodes to elements
	for (unsigned long i = 0; i < m_msh->ele_vector.size(); i++)
	{
		m_element = m_msh->ele_vector[i];
		m_element->GetNodes(vec_element_nodes);
		delta_gas_dissolved = 0.0;
		//cout << endl;
		for (long j = 0; j < long(vec_element_nodes.Size()); j++)
		{
			//cout << "     Node: " << vec_element_nodes[j]->GetIndex() << " RS: " << this->NodeData[vec_element_nodes[j]->GetIndex()]->deltaDIC << endl;
			weight = float(1.0 / vec_element_nodes.Size()); //arithmetic average
			delta_gas_dis = this->NodeData[vec_element_nodes[j]->GetIndex()]->deltaDIC;
			delta_gas_dissolved += delta_gas_dis * weight;
		}
		if (delta_gas_dissolved < 0)
			//check if delta_gas_dissolved is not larger than the available amount of CO2
			if (fabs(delta_gas_dissolved) >
			    (this->Data[i][variable_index_Gas_dissolved] - epsilon))
			{
				cout <<
				"The amount of CO2 removed from water is larger than the available amount of CO2! The value is corrected!"
				     << endl;
				for (long j = 0; j < long(vec_element_nodes.Size()); j++)
					cout << "     Node: " <<
					vec_element_nodes[j]->GetIndex() << " RS: " <<
					this->NodeData[vec_element_nodes[j]->GetIndex()]->deltaDIC
					<< endl;
				cout << "  Element: " << i << " deltaDIC: " <<
				delta_gas_dissolved << " RS_alt: " <<
				this->Data[i][variable_index_Gas_dissolved];
				delta_gas_dissolved = -this->Data[i][variable_index_Gas_dissolved];
			}
		if (delta_gas_dissolved != 0)
		{
			//cout << "  Element: " << i << " deltaDIC: " << delta_gas_dissolved << " RS_alt: " << this->Data[i][variable_index_Gas_dissolved];
			if ((fabs(delta_gas_dissolved) /
			     this->Data[i][variable_index_Gas_dissolved]) < epsilon)
			{
				vecRS.push_back(this->Data[i][variable_index_Gas_dissolved]);
				if (this->E100 == false)
					//stores the change of the total moles of CO2 if e300 is used
					vec_T_CO2.push_back(this->Data[i][variable_index_Moles_CO2]);
			}
			else
			{
				vecRS.push_back(
				        this->Data[i][variable_index_Gas_dissolved] +
				        delta_gas_dissolved);
				if (this->E100 == false)
				{
					//stores the change of the total moles of CO2 if e300 is used
					//1. step: calculate mole of dissolved CO2 in m� liquid
					//c_CO2inLiquid [mol/m�liquid] = RS [m�CO2 / m�liquid] * rho_CO2 [kgCO2 / m�CO2] * 1000 [gCO2 / kgCO2] / M_CO2 [gCO2 / molCO2]
					deltaT_CO2 = delta_gas_dissolved *
					             this->SurfaceCO2Density * 1000 /
					             this->Molweight_CO2;
					//2. step: calculate total moles of CO2
					//T_CO2 [mol] = c_CO2inLiquid [molCO2 / m�liquid] * s_CO2 [m�liquid / m�PV] * V_PV [m�PV]
					//deltaT_CO2 = deltaT_CO2 * this->Data[i][variable_index_water_saturation] * this->Data[i][variable_index_porevolume];
					//cout << " water saturation: " << this->Data[i][variable_index_water_saturation] << " pore volume: " << this->Data[i][variable_index_porevolume] << endl;
					//T_CO2 [mol / L PV] = c_CO2inLiquid [molCO2 / m�liquid] * s_CO2 [m�liquid / m�PV] * 0.001 [m�PV / L PV]
					deltaT_CO2 = deltaT_CO2 *
					             this->Data[i][variable_index_water_saturation]
					             * 0.001;
					if (deltaT_CO2 < 0)
						if (fabs(deltaT_CO2) >
						    (this->Data[i][variable_index_Moles_CO2] -
						     epsilon))
							deltaT_CO2 =
							        -this->Data[i][
							                variable_index_Moles_CO2];
					deltaT_CO2 = this->Data[i][variable_index_Moles_CO2] +
					             deltaT_CO2;
					vec_T_CO2.push_back(deltaT_CO2);
				}
			}
		}
		else
		{
			vecRS.push_back(this->Data[i][variable_index_Gas_dissolved]);
			if (this->E100 == false)
				//stores the change of the total moles of CO2 if e300 is used
				vec_T_CO2.push_back(this->Data[i][variable_index_Moles_CO2]);
		}

		if (vecRS[vecRS.size() - 1] < 0)
		{
			cout << "The new calculated dissolved amount of CO2 (RS) is negative!" <<
			endl;
			system("Pause");
			return 0;
		}
		//cout  << "  Element: " << i << " RS_neu: " << vecRS[vecRS.size()-1] << endl;
	}

	if (this->E100 == true)
	{
		//write RS into the Eclipse restart file
		Filename = folder + "TemporaryResults.F" + AddZero(m_pcs->Tim->step_current - 1,
		                                                   4,
		                                                   true);
		Keyword = " 'RS";

		//check consitency of data
		if (int(vecRS.size()) != this->elements)
		{
			cout <<
			"There was not the right number of RS data found during writing to Eclipse F-File!"
			     << endl;
			system("Pause");
			return 0;
		}
		for (long i = 0; i < this->elements; i = i + 4)
		{
			tempstring = "";
			j_max = 4;
			if ((this->elements - i) < 4)
				j_max = (this->elements - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += "   ";
				tempNumber = this->RoundEXP(vecRS[i + j], 8);
				tempstring += this->AddZero(tempNumber.Number, 8, false);
				tempstring += "E";
				if (tempNumber.Exponent >= 0)
				{
					tempstring += "+";
					tempstring += this->AddZero(tempNumber.Exponent, 2, true);
				}
				else
					tempstring += this->AddZero(tempNumber.Exponent, 3, true);
			}
			vecString.push_back(tempstring);
		}
		if (ReplaceASectionInFile(Filename, Keyword, vecString, false) == false)
		{
			cout << "Replacing a section in the file: " << Filename <<
			" didn't work!" << endl;
			system("Pause");
			return 0;
		}
	}

	if (this->E100 == false)
	{
		Filename = folder + "TemporaryResults.F" + AddZero(m_pcs->Tim->step_current - 1,
		                                                   4,
		                                                   true);
		Keyword = " 'MLSC2";
		vecString.clear();
		//check consitency of data
		if (int(vecRS.size()) != this->elements)
		{
			cout <<
			"There was not the right number of RS data found during writing to Eclipse F-File!"
			     << endl;
			system("Pause");
			return 0;
		}
		for (long i = 0; i < this->elements; i = i + 4)
		{
			tempstring = "";
			j_max = 4;
			if ((this->elements - i) < 4)
				j_max = (this->elements - i);
			for (int j = 0; j < j_max; j++)
			{
				tempstring += "   ";
				tempNumber = this->RoundEXP(vec_T_CO2[i + j], 8);
				tempstring += this->AddZero(tempNumber.Number, 8, false);
				tempstring += "E";
				if (tempNumber.Exponent >= 0)
				{
					tempstring += "+";
					tempstring += this->AddZero(tempNumber.Exponent, 2, true);
				}
				else
					tempstring += this->AddZero(tempNumber.Exponent, 3, true);
			}
			vecString.push_back(tempstring);
		}

		if (ReplaceASectionInFile(Filename, Keyword, vecString, false) == false)
		{
			cout << "Replacing a section in the file: " << Filename <<
			" didn't work!" << endl;
			system("Pause");
			return 0;
		}
	}

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;
	cout << "                    Time: " << time << " seconds." << endl;

	return 1;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: RunEclipse
   Task: Preprocessing, running Eclipse, Postprocessing
   Return: nothing
   Programming: 09/2009 BG / SB
   Modification:
   -------------------------------------------------------------------------*/
int CECLIPSEData::RunEclipse(long Timestep, CRFProcess* m_pcs)
{
	std::string tempstring;
	std::string projectname;
	std::string Filename;
	std::string Filename_Wells;
	ostringstream temp;
	std::string DOScommand;
	std::string folder;
	int position;
	clock_t start,finish;
	clock_t start_execute, finish_execute;
	double time;

	start = clock();

	this->Molweight_CO2 = 44.009; // [g/mol]
	this->Molweight_H2O = 18.0148; // [g/mol]
	this->Molweight_NaCl = 58.443; // [g/mol]

	this->existWells = false;

	// get actual time
	if (m_pcs->Tim->step_current == 1)
		this->actual_time = 0;
	else
		this->actual_time = this->actual_time +
		                    m_pcs->Tim->time_step_vector[m_pcs->Tim->step_current - 2];

	Filename = m_pcs->simulator_model_path; // path to eclipse input data
	if (Timestep == 1)
	{
		this->UsePrecalculatedFiles = m_pcs->PrecalculatedFiles;
		position = int(Filename.find_last_of("\\"));
		if (position >= 0)
			this->Windows_System = true;
		else
			this->Windows_System = false;
	}

	if (this->Windows_System == true)
		position = int(Filename.find_last_of("\\"));
	else
		position = int(Filename.find_last_of("/"));
	folder = Filename.substr(0, position + 1);

	if (UsePrecalculatedFiles == false)
		cout << "      Delete old result files " << endl;

	////Delete Resultfiles
	//projectname = Filename.substr(0, Filename.length() - 5);
	//DOScommand = "del " + projectname +  ".F* /Q";
	//if (system(DOScommand.c_str())){
	//	DisplayMsgLn("Could not delete result files! ");
	//	return 0;
	//}
	//projectname = Filename.substr(0, Filename.length() - 5);
	//DOScommand = "del " + projectname + "_RESTART.F* /Q";
	//if (system(DOScommand.c_str())){
	//	DisplayMsgLn("Could not delete result files! ");
	//	return 0;
	//}
	projectname = "";

	Filename_Wells = m_pcs->simulator_well_path; // path to eclipse input data

	// Read external .well-File
	if (Filename_Wells != "")
	{
		this->ReadWellData(Filename_Wells);
		this->existWells = true;
	}

	// Read Eclipse Data-File
	this->ReadDataFromInputFile(Filename);
	if (this->SurfaceCO2Density <= 0)
	{
		cout <<
		"The CO2 Density at surface conditions was not read properly (rho <= 0!)." << endl;
		system("Pause");
		return 0;
	}

	//write dissolved amount of CO2 and and later maybe also water density back to Eclipse; ToDo: Write new Water density back to Eclipse
	if ((Timestep > 1) && (this->GetVariableIndex("RS") >= 0))
	{
		if (this->UsePrecalculatedFiles == true)
			cout <<
			"Attention! If you run a simulation with precalculated Files you have no interaction between Multiphase Flow and Reactions thus you might get wrong Results for dissolved CO2 in liquid!"
			     << endl;
		else if (this->WriteDataBackToEclipse(m_pcs, folder) == 0)
		{
			cout <<
			"The WriteDataBackToEclipse funktion was not finished properly! The simulation is canceled!"
			     << endl;
			system("Pause");
			return 0;
		}
	}

	cout << "      RunEclipse() called " << endl;

	//Execute Eclipse, try several times in case of problems finding the license
	start_execute = clock();
	std::string DataFilename = "";
	int number_loops = 0;
	int maximum_loops = 10;
	do {
		cout << number_loops + 1 << ". trial" << endl;
		projectname = ExecuteEclipse(Timestep, m_pcs, folder);
		//check if Eclipse has run properly -> errors in the projectname.prt file must be zero
		DataFilename = projectname + ".F" + AddZero(Timestep, 4, true);
		number_loops += 1;
	} while ((CheckIfFileExists(DataFilename) == false) && (number_loops <= maximum_loops));

	finish_execute = clock();
	time = (double(finish_execute) - double(start_execute)) / CLOCKS_PER_SEC;
	cout << endl;
	cout << "  Timestep: " << Timestep << endl;
	cout << "        RunEclipse() called                   Time: " << time << " seconds." <<
	endl;

	if (number_loops > maximum_loops)
	{
		cout << "The Eclipse execution does not work after " << number_loops <<
		" trials!" << endl;
		system("Pause");
		return 0;
	}

	if (Timestep == 1)
	{
		// Read ECLIPSE model grid, check consistency with GeoSys mesh and construct faces
		//Read Eclipse model grid
		std::string GridFilename = "";
		GridFilename = projectname + ".FGRID";
		this->ReadEclipseGrid(GridFilename);

		//Check if Eclipse grid is identical with Geosys grid
		if (CompareElementsGeosysEclipse() == 0)
		{
			cout <<
			"The Eclipse grid is not identical to the Geosys grid! Abort Calculations!"
			     << endl;
			system("Pause");
			return 0;
		}

		this->DetermineNeighbourElements(projectname);

		//Read boundary conditions
		tempstring = projectname + ".data";

		//this->ReadPositionBoundaryCondition(tempstring);
		//Create Faces
		this->CreateFaces();
		// Connect Faces to each block
		this->ConnectFacesToElements();
		// check for radial model and one column
		int count = 0, found = 0, iindex = 0;
		double* nvec;
		if(this->Radialmodell == true)
		{
			if(this->rows > 1) // check if more than one column is active in radial model
			{
				for(long ii = 0; ii < long(this->eclgrid.size()); ii++)
				{
					if((this->eclgrid[ii]->layer == 1) &&
					   (this->eclgrid[ii]->column == 1))
					{
						found++;
						if(this->eclgrid[ii]->active == 1)
						{
							count++;
							iindex = ii;
						}
					}
					if(found == this->rows)
						break;  // foun all cells in layer ==1 and column == 1
				}
				if( count > 1)
					cout << endl << endl <<
					" Error in definition of radial flow model; Use not more than one active column "
					     << endl;
			}
			// test if faces are perpenducular to I coordinate axis
			for(int ii = 0; ii < int(this->eclgrid[iindex]->connected_faces.size());
			    ii++)
				if(this->faces[ii]->model_axis.find("I") == 0)
				{
					nvec = this->faces[ii]->PlaneEquation->GetNormalVector();
					if(nvec[0] == 1.0)
						this->Radial_I = true;
				}
			if(Radial_I != true)
				cout << endl <<
				" Error: I  Face is not perpendicular to coordinate axis " << endl;
		}
	}

	//Read Eclipse model data
	DataFilename = projectname + ".F" + AddZero(Timestep, 4, true);
	long EclipseTimeStep = 1;
	this->ReadEclipseData(DataFilename, EclipseTimeStep - 1);

	//calculate phase flux if E300 is used
	if (this->E100 == false)
		this->CalculateRSfromMassFraction_E300();

	//this->MakeNodeVector(m_pcs, folder, EclipseTimeStep - 1, 0);
	this->MakeNodeVector();

	for (int i = 0; i < int(this->Phases.size()); i++)
	{
		//Get the flow for each face
		this->GetFlowForFaces(i);
		//this->MakeNodeVector(m_pcs, folder, EclipseTimeStep - 1, 0);
		//this->CalcBlockBudget(i);
		//this->GetVelForFaces();
		//interpolate face flows to nodes
		this->InterpolateDataFromBlocksToNodes(m_pcs, folder, i);
	}
	WriteDataToGeoSys(m_pcs, folder);

	if (UsePrecalculatedFiles == false)
		CleanUpEclipseFiles(folder, projectname);

	finish = clock();
	time = (double(finish) - double(start)) / CLOCKS_PER_SEC;

	cout << "        Time for this timestep: " << time << " seconds." << endl;

	return 1;
}

//TODO:

// interpolation of pressure to nodes-> volume is just devided by the number of nodes -> not exact if shape is not a quader
