/*
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef MESHQUALITYVOLUMES_H_
#define MESHQUALITYVOLUMES_H_

#include "MeshQualityChecker.h"

namespace MeshLib
{
class MeshQualityVolumes : public MeshQualityChecker
{
public:
	MeshQualityVolumes(CFEMesh const* const mesh);
	virtual ~MeshQualityVolumes() {}

	virtual void check ();
protected:
	double _minimum, _maximum;
};
}

#endif /* MESHQUALITYVOLUMES_H_ */
