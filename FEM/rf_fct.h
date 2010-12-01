/**************************************************************************
FEMLib-Object: FCT
Task: Functions
Programing:
01/2005 OK Implementation
**************************************************************************/
#ifndef rf_fct_INC
#define rf_fct_INC
// C++ STL
#include <list>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

class CFunction
{
  private:
  public:
    CFunction(void);
    ~CFunction(void);
    void Set(string,string,double);
    ios::pos_type Read(ifstream*);
    void Write(fstream*);
    double GetValue(double point, bool*valid, int method =0);// add int method. WW
  public:
    string type_name;
    string geo_name;
    string geo_type_name;
    string dis_type_name;//CMCD
    vector<string>variable_names_vector;
    vector<double*>variable_data_vector;
//	int matrix_dimension_x; //NB4703
//	int matrix_dimension_y; //NB4703
	vector<int>matrix_dimension; //NB 4.8.01
    bool selected;
};

extern vector<CFunction*>fct_vector;
extern void FCTWrite(string);
extern void FCTRead(string);
extern void FCTReadTIMData(string);
extern CFunction* FCTGet(string);
extern CFunction* FCTGet(long); //YD

#define FCT_FILE_EXTENSION ".fct"

#endif
