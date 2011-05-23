// Classes for interface GeoSys - DUMUX
#include <vector>

//class CRFProcess;
//#include "rf_pcs.h"
#include "fem_ele_std.h"
#include "fem_ele.h"

class CReadTextfiles_DuMux {

public:
    std::vector <std::string> Data;
    std::vector <std::vector <std::string> > Data_separated;
    long NumberOfRows;
    std::vector <std::string> SplittedString;
    std::vector <std::string> Header;

	CReadTextfiles_DuMux();		//Konstruktor
	~CReadTextfiles_DuMux();		//Desturktor

	bool Read_Text(std::string Filename);

	void SplitStrings(const std::string str, std::string delimiter);

	bool Read_SeparatedText(std::string Filename, std::string delimiter);
};

class CWriteTextfiles_DuMux {

public:
	CWriteTextfiles_DuMux();		//Konstruktor
	~CWriteTextfiles_DuMux();		//Desturktor

	void Write_Text(std::string Filename, std::vector<std::string> Text);
};


class CPointData_DuMux {
public:
	double x;
	double y;
	double z;
	double temperature;
	double CO2inLiquid;
	double NaClinLiquid;
	std::vector <double> phase_pressure;
	std::vector <double> phase_saturation;
	std::vector <double> phase_density;
	std::vector <std::vector <double> > q;

	CPointData_DuMux() {x = 0; y = 0; z = 0;}
	~CPointData_DuMux() {}
};


class CDUMUXData {

public:
	CDUMUXData();
	~CDUMUXData();

	std::vector <CPointData_DuMux*> NodeData;
	std::vector <std::string> Phases;
	int dim;
	int ProcessIndex_CO2inLiquid;
	int ProcessIndex_NaClinLiquid;
	bool Windows_System;
	bool UsePrecalculatedFiles;
	double Molweight_CO2;		// [g/mol]
	double TotalSimulationTime;

	//CFiniteElementStd* GetAssembler() {return fem; }

	bool CheckIfFileExists(std::string strFilename);

	std::string AddZero(double Number, int Places, bool before);

	bool MakeNodeVector(void);

	void ExecuteDuMux(CRFProcess *m_pcs, std::string folder);

	int WriteInputForDuMux(CRFProcess *m_pcs, std::string Pathname, long Timestep);

	void ReadDuMuxData(CRFProcess *m_pcs, std::string Pathname, long Timestep);

	void WriteDataToGeoSys(CRFProcess *m_pcs);

	int RunDuMux(long Timestep, CRFProcess *m_pcs);

};
