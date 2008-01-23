#ifndef ACDTILEFUNCS_H
#define ACDTILEFUNCS_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Matrix/Matrix.h"

#include "geometry/Point.h"
#include "geometry/Vector.h"

//#include "idents/AcdId.h"
#include "AcdUtil/AcdTileDim.h"

/**
 * @brief Some functions used in projecting tracks to ACD tiles
 *
 **/

namespace AcdTileUtil {

  /**
   * @brief Project the track error onto a plane
   *
   * FIMXE.  This is broken
   *
   * @param activeX
   * @param activeY
   * @param covXX
   * @param covYY
   * @param planeError
   */
  void planeErrorProjection(const double& activeX, const double& activeY, const double& covXX, const double& covYY,
			    double& planeError);
  /**
   * @brief Determine the Z value of a particular charge deposit
   *
   * FIMXE.  This is broken
   *
   * @param tile
   * @param activeX
   * @param activeY
   * @param covXX
   * @param covXY
   * @param covYY
   * @param doca
   * @param docaErr
   * @param iHole
   */ 
  void tileScrewHoleDoca(const AcdTileDim& tile, const double& activeX, const double& activeY, 
			 const double& covXX, const double& covYY, const double& covXY,
			 double& doca, double& docaErr, int& iHole );
 
}

#endif
