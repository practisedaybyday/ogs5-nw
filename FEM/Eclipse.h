// Classes for interface GeoSys - ECLIPSE
#include <vector>
//#define MAX_ZEILEN 512
//using namespace std;
class CRFProcess;
#include "rf_pcs.h"
#include "msh_faces.h"

typedef struct
{
	double Number;
	int Exponent;
}typeExponentialNumber;
// KB: neue Struktur zum Einlesen der Brunnendatei, siehe CECLIPSEData::ReadWellData()
struct structWell { /* deklariert den Strukturtyp well, wird beim Einlesen von .well-File neu aufgerufen */
	std::string name;
    vector <double> time;
    vector <string> rate;
	std::string phase;
	std::string open_flag;
	std::string control_mode;
	//std::string dummy_rate;
	//std::string dummy_zeile;
	//std::string rest_zeile;
};

class CECLIPSEBlock {

public:
	long index;
	long row;
	long column;
	long layer;
	vector <double> x_coordinates;
	vector <double> y_coordinates;
	vector <double> z_coordinates;
	double x_barycentre;
	double y_barycentre;
	double z_barycentre;
	int active;
	double volume;
	vector <long> connected_faces;
	
	vector <long> NeighbourElement;
	vector <long> ConnectedBoundaryCondition;

	CECLIPSEBlock(long Nodelength, long Facelength);
	~CECLIPSEBlock();

	//void CalcBarycentre(void);
	void CalculateFaceCentres(void);
};

class CReadTextfiles_ECL {

public:
    vector <string> Data;
	vector <vector <string> > Data_separated;
    long NumberOfRows;
	vector <string> SplittedString;
	vector <string> Header;

	CReadTextfiles_ECL();		//Konstruktor
	~CReadTextfiles_ECL();		//Desturktor

	bool Read_Text(string Filename);

	void SplitStrings(const string str, string delimiter);

	bool Read_SeparatedText(string Filename, string delimiter);

	//vector <string> Header;
 //   vector <string> Data;
 //   long NumberOfRows;

	//CReadTextfiles();		//Konstruktor
	//~CReadTextfiles();		//Desturktor

	//bool Read_Text(std::string Filename);
};

class CWriteTextfiles_ECL {

public:
	CWriteTextfiles_ECL();		//Konstruktor
	~CWriteTextfiles_ECL();		//Desturktor

	void Write_Text(std::string Filename, vector<std::string> Text);
};

class CPointData_ECL {
public:
	double x;
	double y;
	double z;
	//double Flow[3];
	double pressure;
	double temperature;
	//double Gas_dissolved;
	double CO2inLiquid;
	double NaClinLiquid;
	double deltaDIC;
	vector <double> phase_pressure;
	vector <double> phase_saturation;
	vector <double> phase_density;
	vector <vector <double> > q;

	CPointData_ECL() {x = 0; y = 0; z = 0;}
	~CPointData_ECL() {}
};

class CBoundaryConditions {
public:
	int index;
	int number;
	long connected_element;
	std::string boundary_position;
	double value[4];
	CBoundaryConditions() {index = -1; number = -1; connected_element = -1; boundary_position = "";}
	~CBoundaryConditions() {}
};

class CECLIPSEData {
public:
	long elements;
	long rows;
	long columns;
	long layers;
	long times;
	int numberOutputParameters;
	int activeCells;
	bool Radialmodell;
	bool Radial_J;
	bool Radial_I;
	bool RadialModellIpos;
	bool RadialModellJpos;
	double Molweight_CO2;		// [g/mol]
	double Molweight_H2O;		// [g/mol]
	double Molweight_NaCl;		// [g/mol]
	double SurfaceCO2Density;
	bool E100;
	double sumCO2removed;
	int ProcessIndex_CO2inLiquid;
	double actual_time;
	bool Windows_System;
	bool existWells;
	bool UsePrecalculatedFiles;

	vector <CECLIPSEBlock*> eclgrid;
	vector <string> SplittedString;
	vector <string> Variables;
	vector <CFaces*> faces;
	vector <CPointData_ECL*> NodeData;
	vector <CBoundaryConditions*> BC;
	vector <structWell*> ecl_well;

	double ** Data;	//array of points, times and variables
	vector <string> WellRates; // KB, abspeichern der neuen Raten für einen bestimmten Zeitschritt

	vector <long> output_x;
	vector <long> output_y;
	vector <long> output_z;
	vector <long> output_time;
	vector <long> CorrespondingEclipseElement;
    vector <long> CorrespondingGeosysElement;
	vector <string> Phases;
	vector <string> Components;
	long a [8][2];//2D Array um Keywords abzuspeichern

	CECLIPSEData();
	~CECLIPSEData();

	int GetVariableIndex(std::string Variablename);

	void SplitStrings(const std::string str,
                      std::string delimiter);

	double Round(double Number, int Decimalplaces);

	typeExponentialNumber RoundEXP(double Number, int Decimalplaces);

	std::string AddZero(double Number, int Places, bool before);

	bool CheckIfFileExists(std::string strFilename);

	bool ReplaceASectionInFile(std::string Filename, std::string Keyword, vector <std::string> Data, bool CheckLengthOfSection);

	bool ReplaceWellRate(std::string Filename, std::string Keyword_well);

	int WriteDataBackToEclipse(CRFProcess *m_pcs, std::string projectname);

	std::string ExecuteEclipse(long Timestep, CRFProcess *m_pcs, std::string folder);

	void ReadEclipseGrid(std::string Filename);

	void DetermineNeighbourElements(std::string Filename);

	bool ReadBoundaryData(int index_boundary, vector <string> Data);

	int ReadDataFromInputFile(std::string Filename);

	bool ReadPositionBoundaryCondition(std::string Filename);

	bool CorrespondingElements(void);
	
	bool CompareElementsGeosysEclipse(void);

	double CalculateDistanceBetween2Points(double Point1[3], double Point2[3]);

	bool CreateFaces(void);

	bool ConnectFacesToElements(void);

	//bool MakeNodeVector(CRFProcess *m_pcs, std::string path, int timestep, int phase_index);
	bool MakeNodeVector(void);

	void ReadEclipseData(std::string Pathname, long timestep);

	void CalculateRSfromMassFraction_E300();

	bool GetFlowForFaces(int phase_index);

	bool GetVelForFaces(void);

	bool CalcBlockBudget(int phase_index);

	void InterpolateDataFromFacesToNodes(long ele_nr, double* n_vel_x, double* n_vel_y, double* n_vel_z, int phase_index);

	void InterpolateDataFromBlocksToNodes(CRFProcess *m_pcs, std::string path, int phase_index);

	void InterpolateGeosysVelocitiesToNodes(CRFProcess *m_pcs, double *vel_nod, long node);

	void WriteDataToGeoSys(CRFProcess *m_pcs, std::string folder);

	bool CleanUpEclipseFiles(std::string folder, std::string projectname);
	
	int RunEclipse(long Timestep, CRFProcess *m_pcs);

	void ReadWellData(std::string Filename_Wells);
};
