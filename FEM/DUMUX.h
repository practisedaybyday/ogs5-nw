// Classes for interface GeoSys - DUMUX
#include <vector>
using namespace std;
//class CRFProcess;
//#include "rf_pcs.h"
#include "fem_ele_std.h"
#include "fem_ele.h"
//#include "Eclipse.h"
//namespace FiniteElement {class CFiniteElementStd; class CFiniteElementVec;
//                         class ElementMatrix; class ElementValue;}

class CReadTextfiles_DuMux {

public:
    vector <string> Data;
	vector <vector <string> > Data_separated;
    long NumberOfRows;
	vector <string> SplittedString;
	vector <string> Header;

	CReadTextfiles_DuMux();		//Konstruktor
	~CReadTextfiles_DuMux();		//Desturktor

	bool Read_Text(string Filename);

	void SplitStrings(const string str, string delimiter);

	bool Read_SeparatedText(string Filename, string delimiter);
};

class CWriteTextfiles_DuMux {

public:
	CWriteTextfiles_DuMux();		//Konstruktor
	~CWriteTextfiles_DuMux();		//Desturktor

	void Write_Text(string Filename, vector<string> Text);
};


class CPointData_DuMux {
public:
	double x;
	double y;
	double z;
	double temperature;
	double CO2inLiquid;
	double NaClinLiquid;
	vector <double> phase_pressure;
	vector <double> phase_saturation;
	vector <double> phase_density;
	vector <vector <double> > q;

	CPointData_DuMux() {x = 0; y = 0; z = 0;}
	~CPointData_DuMux() {}
};


class CDUMUXData {

public:
	CDUMUXData();
	~CDUMUXData();

	vector <CPointData_DuMux*> NodeData;
	vector <string> Phases;
	int dim;
	int ProcessIndex_CO2inLiquid;
	int ProcessIndex_NaClinLiquid;
	bool Windows_System;
	bool UsePrecalculatedFiles;
	double Molweight_CO2;		// [g/mol]
	double TotalSimulationTime;

	//CFiniteElementStd* GetAssembler() {return fem; }

	bool CheckIfFileExists(string strFilename);

	string AddZero(double Number, int Places, bool before);

	bool MakeNodeVector(void);

	void ExecuteDuMux(CRFProcess *m_pcs, string folder);

	int WriteInputForDuMux(CRFProcess *m_pcs, string Pathname, long Timestep);

	void ReadDuMuxData(CRFProcess *m_pcs, string Pathname, long Timestep);

	void WriteDataToGeoSys(CRFProcess *m_pcs);

	int RunDuMux(long Timestep, CRFProcess *m_pcs);

};