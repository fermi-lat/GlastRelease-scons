#include "../AcdRecon/AcdPatRecTools.h"

#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"

#include "idents/AcdId.h"


namespace AcdRecon {

  /// Turn this on to dump the table
  bool s_dumpRoiLookupTable = false;

  /**
   * @brief Project a Global point into the flattened ACD coords
   *
   * @param point is the point in global coords
   * @param acdUV are the projected ACD U&V coordinate
   * @return true for success, false for failure
   **/
  bool projectGlobalToUV(const HepPoint3D& point,
			 AcdUVCoord& acdUV) {
    double dist2 = (point.x() * point.x()) + (point.y() * point.y());
    double dist = sqrt(dist2);
    double zVal = point.z() + s_projZRef;
    double rad = atan2( dist , zVal );
    double aU = rad * ( point.x() / dist );
    double aV = rad * ( point.y() / dist );
    acdUV.set(aU,aV);
    return true;
  }

  /**
   * @brief Convert U or V to a bin
   *
   * @param val is the U or V value
   * @return bin number, -1 for value out of range
   **/
  int convertToUVBin(double val) {
    val += s_UVMax;
    val /= ( 2. * s_UVMax );
    val *= ((float)(s_nBins));
    int retVal = (int)(floor(val));
    if ( retVal < 0 ) return -1;
    if ( retVal >= s_nBins ) return -1;
    return retVal;
  }

  /**
   * @brief Convert UV coordiantes to bins numbers
   *
   * @param acdUV are the projected ACD U&V coordinate
   * @param binUV are the bin numbers in U&V
   * @return true for success, false for failure
   **/
  bool convertUVtoBins(const AcdUVCoord& acdUV,
		       AcdUVBin& binUV) {
    double aU = acdUV.u();
    double aV = acdUV.v();
    int binU = convertToUVBin(aU);
    int binV = convertToUVBin(aV);
    if ( binU < 0 || binV < 0 ) return false;
    binUV.set(binU,binV);
    return true;
  }

  /**
   * @brief increase a vector by toleracne
   *
   * @param inVect is the input vector
   * @param outVect is the lengthend vector
   * @param tolerence is the tolerence added to the frame vectors
   **/
  void increaseFrameVector(const HepVector3D& inVect,
			   HepVector3D& outVect,
			   double tolerance ) {
    double initSize = inVect.mag();
    double finalSize = initSize + tolerance;
    double scale = finalSize / initSize;
    outVect = inVect;
    outVect *= scale;
  }
  
  
  /**
   * @brief build the set of all bins touched by a given tile
   *
   * @param center is the tile center, in Global Coords
   * @param frameX is the localX vector, in Global Coords
   * @param frameY is the localY vector, in Global Coords
   * @param binsSet is the set of all the bins near this volume, not cleared
   * @param tolerence is the tolerence added to the frame vectors
   * @return true for success, false for failure
   **/
  bool buildBinSet(const HepPoint3D &center,
		   const HepVector3D& frameX,
		   const HepVector3D& frameY,
		   std::set<AcdUVBin>& binsSet,
		   double tolerance, int acdId ) {

    // first, increase the frame vectors by the tolerance
    HepVector3D vX;
    increaseFrameVector(frameX,vX,tolerance);
    HepVector3D vY;
    increaseFrameVector(frameY,vY,tolerance);
    AcdUVCoord acdUV;
    AcdUVBin binUV;
    float xStep(0.1);
    float yStep(0.1);
    switch ( acdId ) {
    case 500:
    case 501:
    case 502:
    case 503:
    case 600:
    case 601:
    case 602:
    case 603:
      // ribbons, more steps x, in fewer in y
      xStep = 0.002; yStep = 0.1;
      break;
    case 130:
    case 230:
    case 330:
    case 430:
      // large tiles at bottom of sides, more steps in X
      xStep = 0.002; yStep = 0.005;
      break;
    default:
      // all others use defaults
      xStep = 0.005; yStep = 0.005;
    }
    if ( s_dumpRoiLookupTable ) {
      std::cout << "ACD ROI " << acdId 
		<< center.x() << ' ' << center.y() << ' ' << center.z() << ' '
		<< vX.x() << ' ' << vX.y() << ' ' << vX.z() << ' ' 
		<< vY.x() << ' ' << vY.y() << ' ' << vY.z() << std::endl;
    }
    for ( float dx(-1.); dx <= 1.001; dx += xStep ) {
      for ( float dy(-1.); dy <= 1.001; dy += yStep ) {
	HepPoint3D p = center;
	HepVector3D delX = vX;  delX *= dx;       
	HepVector3D delY = vY;  delY *= dy;
	p += delX;
	p += delY;
	if ( ! projectGlobalToUV( p, acdUV ) ) return false;
	if ( ! convertUVtoBins( acdUV, binUV ) ) return false;
	if ( s_dumpRoiLookupTable ) {
	  std::cout << "ACD ROI TABLE " << acdId << ' ' 
		    << binUV.binU() << ' ' << binUV.binV() << ' ' 
		    << dx << ' ' << dy << ' ' 
		    << acdUV.u() << ' ' << acdUV.v() << ' ' 
		    << p.x() << ' ' << p.y() << ' ' << p.z() << std::endl;
	}	
	binsSet.insert( binUV );
      }
    }    
    return true;
  }

  /**
   * @brief Add a Tile to the pat rec maps
   *
   * @param tile the geom. info about the tile in question
   * @param theMap is the Pat. Rec. hash map
   * @param revMap is the reverse Pat. Rec. hash map (mainly for debbugging)
   * @param tolerence is the tolerence added to the frame vectors
   * @return true for success, false for failure
   **/
  bool addTileToPatRecMaps(const AcdTileDim& tile,
			   AcdPatRecMap& theMap,
			   AcdPatRecMap_Rev& revMap,
			   double tolerance ) {

    static const HepVector3D xHat(1.,0.,0.);
    static const HepVector3D yHat(0.,1.,0.);

    std::set<AcdUVBin> binsSet;
    const idents::AcdId& acdId = tile.acdId();
    for ( int i(0); i < tile.nVol(); i++ ) {
      const AcdTileSection* sect = tile.getSection(i);
      const HepPoint3D& center = sect->m_center;

      HepVector3D frameX = sect->m_invTrans * xHat;
      HepVector3D frameY = sect->m_invTrans * yHat;
      frameX *= (sect->m_dim[0] / 2.);
      frameY *= (sect->m_dim[1] / 2.);
      if ( ! buildBinSet( center, frameX, frameY, binsSet, tolerance, acdId.id() ) ) {
	std::cerr << "buildBinSet failed at " << acdId.id() << ' ' << i << ' ' << frameX << ' ' << frameY << std::endl;
	return false; 
      }
    }
    for ( std::set<AcdUVBin>::const_iterator itr = binsSet.begin();
	  itr != binsSet.end(); itr++ ) {
      if ( s_dumpRoiLookupTable ) {
	std::cout << "Adding " << itr->binU() << ':' << itr->binV() << " <-> " << acdId.id() << std::endl;
      }
      theMap[*itr].push_back(acdId);
      revMap[acdId].push_back(*itr);
    }
    return true;
  }

  /**
   * @brief Add a Ribbon to the pat rec maps
   *
   * @param ribbon the geom. info about the ribbon in question
   * @param theMap is the Pat. Rec. hash map
   * @param revMap is the reverse Pat. Rec. hash map (mainly for debbugging)
   * @param tolerence is the tolerence added to the frame vectors
   * @return true for success, false for failure
   **/
  bool addRibbonToPatRecMaps(const AcdRibbonDim& ribbon,
			     AcdPatRecMap& theMap,
			     AcdPatRecMap_Rev& revMap,
			     double tolerance ) {
    std::set<AcdUVBin> binsSet;
    const idents::AcdId& acdId = ribbon.acdId();
    static const HepVector3D unitX(1.,0.,0.);
    static const HepVector3D unitY(0.,1.,0.);
    for ( int i(AcdRibbonDim::MinusSide); i <= AcdRibbonDim::PlusSide; i++ ) {
      HepPoint3D center;
      HepVector3D frameX;
      if ( ! ribbon.setEdgeRay( i, center, frameX ) ) {
	std::cerr << "setEdgeRay failed for " << acdId.id() << std::endl;
	return false;
      }
      frameX *= 0.5;
      center += frameX;
      const HepVector3D& frameY = acdId.faceLike() == 5 ? unitY : unitX;
      if ( ! buildBinSet( center, frameX, frameY, binsSet, tolerance, acdId.id() ) ) {
	std::cerr << "buildBinSet failed at " << acdId.id() << ' ' << frameX << ' ' << frameY << std::endl;
	return false;      
      }
    }
    for ( std::set<AcdUVBin>::const_iterator itr = binsSet.begin();
	  itr != binsSet.end(); itr++ ) {
      if ( s_dumpRoiLookupTable ) {
	std::cout << "Adding " << itr->binU() << ':' << itr->binV() << " <-> " << acdId.id() << std::endl;
      }
      theMap[*itr].push_back(acdId);
      revMap[acdId].push_back(*itr);
    }
    return true;
  }

  /**
   * @brief Get the list of detector elements that with have been hit by a given track
   *
   * @param point, where the track exited the ACD
   * @param theMap is the Pat. Rec. hash map
   * @return point to the list of AcdId, null for failure
   **/
  const std::list<idents::AcdId>* findNearbyElements(const HepPoint3D& p,
						     AcdPatRecMap& theMap) {
    AcdUVCoord acdUV;
    AcdUVBin binUV;
    if ( ! projectGlobalToUV(p,acdUV) ) return 0;
    if ( ! convertUVtoBins(acdUV,binUV) ) return 0;
    AcdPatRecMap::iterator itr = theMap.find(binUV);
    if ( itr == theMap.end() ) return 0;
    const std::list<idents::AcdId>& theList = itr->second;
    return &theList;
  }


  bool buildElementSet(const HepPoint3D& p,
		       const HepVector3D& v,
		       const double& arcTol,
		       std::set<idents::AcdId>& idSet,
		       AcdPatRecMap& theMap) {
    
    std::set< const std::list<idents::AcdId>* > setList;      
    HepVector3D vScaled = arcTol*v;
    HepPoint3D p0 = p - vScaled;
    HepPoint3D p2 = p + vScaled;
    const std::list<idents::AcdId>* l0 = findNearbyElements(p0,theMap);
    if ( l0 != 0 ) setList.insert(l0);
    const std::list<idents::AcdId>* l1 = findNearbyElements(p,theMap);
    if ( l1 != 0 ) setList.insert(l1);
    const std::list<idents::AcdId>* l2 = findNearbyElements(p2,theMap);
    if ( l2 != 0 ) setList.insert(l2);
    
    for ( std::set< const std::list<idents::AcdId>* >::const_iterator itr = setList.begin();
	  itr != setList.end(); itr++ ) {
      const std::list<idents::AcdId>* aList = *itr;
      for ( std::list<idents::AcdId>::const_iterator itrId = aList->begin();
	    itrId != aList->end(); itrId++ ) {
	idSet.insert(*itrId);
      }
    }
    return true;
  }
  

}
