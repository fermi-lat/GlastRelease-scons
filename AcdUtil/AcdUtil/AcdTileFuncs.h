#ifndef ACDTILEFUNCS_H
#define ACDTILEFUNCS_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/Matrix.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

//#include "idents/AcdId.h"
#include "AcdUtil/AcdTileDim.h"

namespace AcdTileUtil {


  void planeErrorProjection(const double& activeX, const double& activeY, const double& covXX, const double& covYY,
			    double& planeError);
 
  void tileScrewHoleDoca(const AcdTileDim& tile, const double& activeX, const double& activeY, 
			 const double& covXX, const double& covYY, const double& covXY,
			 double& doca, double& docaErr, int& iHole );
 
}

#endif
