#include "msh_core.h"
#include <string>
#include <vector>
#include "msh_node.h"

class CPlaneEquation {
private:
	double Point[3];
	double vector1[3];
	double vector2[3];
	double Lambda_NormalEquation;
	double normal_vector[3];

public:
	CPlaneEquation();
	~CPlaneEquation();

	void CalculatePlaneEquationFrom3Points(double Point1[3], double Point2[3], double Point3[3]);

	double* GetNormalVector() {return normal_vector;}

	bool CheckIfPointInPlane(double Point[3]);
};

class CFlowData {
public:
	double q[3];
	double q_norm;
	CFlowData() {q_norm = 0; q[0] = 0; q[1] = 0; q[2] = 0;}
	~CFlowData() {}
};

// Class definition
class CFaces // necessary for coupling with Eclipse
{
private:
	double gravity_centre[3];
	long nnodes;
	std::vector <MeshLib::CNode*> connected_nodes;

public:
	CPlaneEquation *PlaneEquation;
	long index;
	std::vector <long> connected_elements;
	std::string model_axis;

	//double q[3];
	//double q_norm;
	double v_norm;
	double vel[3];

	std::vector <CFlowData*> phases;

	double face_area;

	CFaces(int number_phases);
	~CFaces();

	bool Calculate_FaceGravityCentre(double Point1[3], double Point2[3], double Point3[3], double Point4[3]);

	double* GetFaceGravityCentre() {return gravity_centre;}

	void SetNodes(MeshLib::CNode* Point1, MeshLib::CNode* Point2, MeshLib::CNode* Point3, MeshLib::CNode* Point4);

	void SetElements(std::vector <long> element_indices);

	bool CreateFace(MeshLib::CNode* Point1, MeshLib::CNode* Point2, MeshLib::CNode* Point3, MeshLib::CNode* Point4);

	void Calculate_components_of_a_vector(int flag, int phase_index, bool Radialmodell);

};
