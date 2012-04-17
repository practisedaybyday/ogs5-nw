/*
 * MeshQualityNormalisedVolumes.h
 *
 *  Created on: Mar 3, 2011
 *      Author: TF
 */

#ifndef MESHQUALITYNORMALISEDVOLUMES_H_
#define MESHQUALITYNORMALISEDVOLUMES_H_

#include "MeshQualityChecker.h"

namespace MeshLib
{
class MeshQualityNormalisedVolumes : public MeshQualityChecker
{
public:
	MeshQualityNormalisedVolumes(CFEMesh const* const mesh);
	virtual ~MeshQualityNormalisedVolumes() {}

	virtual void check ();
};
}

#endif /* MESHQUALITYNORMALISEDVOLUMES_H_ */
