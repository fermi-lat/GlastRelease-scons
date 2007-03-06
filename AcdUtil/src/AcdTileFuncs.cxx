#include "../AcdUtil/AcdTileFuncs.h"

#include "CLHEP/Matrix/Matrix.h"

namespace AcdTileUtil {

  void planeErrorProjection(const double& activeX, const double& activeY, const double& covXX, const double& covYY,
			    double& planeError){
    
    double sigX2 = (activeX*activeX*covXX);
    double sigY2 = (activeY*activeY*covYY);
    if ( sigX2 > sigY2) {
      planeError = sqrt(covYY);
    } else {
      planeError = sqrt(covYY);
    }
  }

  void tileScrewHoleDoca(const idents::AcdId& tileId, const double& activeX, const double& activeY, 
			 const double& /*covXX*/, const double& /*covYY*/, const double& /*covXY*/,
			 double& doca, double& docaErr, int& iHole) {
    
    // This is the way we really should do it.
    
    //double localX(0.), localY(0.);
    //AcdTileDim::toLocalCoords(tileId,activeX,tileHit.activeY,localX,localY)
    //const std::vector< HepPoint3D >& holes = tile.screwHoles();
    //for ( int i(0); i < holes.size(); i++ ) {
    // }

    // for now, use the fact that screw holes are 5cm from edge of tiles
    double deltaX = activeX - 5.;
    double deltaY = activeY - 5.;
    doca = sqrt(deltaX*deltaX + deltaY*deltaY);
    docaErr = -1.;
    iHole = -1;
  }
 
}
