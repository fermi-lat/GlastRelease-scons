#ifndef ACDPATRECTOOLS_H
#define ACDPATRECTOOLS_H

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Vector3D.h"
#include <map>
#include <list>
#include <set>

namespace idents {
  class AcdId;
}

class AcdTileDim;
class AcdRibbonDim;
class Point;

namespace AcdRecon {
  

  // A couple of constants

  /// Reference point for 2-D projection
  const double s_projZRef = 600.;
  /// Range of UV coordinates
  const double s_UVMax = 1.2;
  /// Number of bins for hash set
  const int s_nBins = 36;

  struct AcdUVCoord {
  public:
    inline AcdUVCoord()
      :m_u(0.),m_v(0.){}
    inline AcdUVCoord(double u, double v)
      :m_u(u),m_v(v){}
    inline ~AcdUVCoord(){}
    inline const double& u() const { return m_u; }
    inline const double& v() const { return m_v; }
    inline const AcdUVCoord& operator=(const AcdUVCoord& other) {
      m_u = other.u();
      m_v = other.v();
      return *this;
    }
    inline void set(const double& aU, const double& aV){
      m_u = aU;
      m_v = aV;
    }
  private:
    double m_u;
    double m_v;
  };

  struct AcdUVBin {
  public:
    inline static int hashVal(int binU, int binV){
      return binU * s_nBins + binV;
    }
  public:
    inline AcdUVBin()
      :m_binU(-1),m_binV(-1),m_hash(-1){;}
    inline AcdUVBin(int binU, int binV)
      :m_binU(binU),m_binV(binV),m_hash(hashVal(binU,binV)){}
    inline AcdUVBin(const AcdUVBin& other)
      :m_binU(other.binU()),m_binV(other.binV()),m_hash(other.hashVal()){}
    inline ~AcdUVBin(){}
    inline const int& binU() const { return m_binU; }
    inline const int& binV() const { return m_binV; }
    inline const int& hashVal() const { return m_hash; }
    inline const AcdUVBin& operator=(const AcdUVBin& other) {
      m_binU = other.binU();
      m_binV = other.binV();
      m_hash = hashVal(m_binU,m_binV);
      return *this;
    }
    inline void set(int bU, int bV) {
      m_binU = bU;
      m_binV = bV;
      m_hash = hashVal(m_binU,m_binV);
    }
    inline bool operator==(const AcdUVBin& other) const {
      return m_hash == other.hashVal();
    }
    inline bool operator!=(const AcdUVBin& other) const {
      return m_hash != other.hashVal();
    }   
    inline bool operator<(const AcdUVBin& other) const {
      return m_hash < other.hashVal();
    }   
    inline bool operator<=(const AcdUVBin& other) const {
      return m_hash <= other.hashVal();
    }   
    inline bool operator>(const AcdUVBin& other) const {
      return m_hash > other.hashVal();
    }
    inline bool operator>=(const AcdUVBin& other) const {
      return m_hash > other.hashVal();
    }   
  private:
    int m_binU;
    int m_binV; 
    int m_hash;
  };

  typedef std::map<AcdUVBin, std::list<idents::AcdId> > AcdPatRecMap;
  typedef std::map<idents::AcdId, std::list<AcdUVBin> > AcdPatRecMap_Rev;


  /**
   * @brief Project a Global point into the flattened ACD coords
   *
   * @param point is the point in global coords
   * @param acdUV are the projected ACD U&V coordinate
   * @return true for success, false for failure
   **/
  bool projectGlobalToUV(const HepPoint3D& point,
			 AcdUVCoord& acdUV);

  
  /**
   * @brief Convert U or V to a bin
   *
   * @param val is the U or V value
   * @return bin number, -1 for value out of range
   **/
  int convertToUVBin(double val);


  /**
   * @brief Convert UV coordiantes to bins numbers
   *
   * @param acdUV are the projected ACD U&V coordinate
   * @param binUV are the bin numbers in U&V
   * @return true for success, false for failure
   **/
  bool convertUVtoBins(const AcdUVCoord& acdUV,
		       AcdUVBin& binUV);

  /**
   * @brief increase a vector by toleracne
   *
   * @param inVect is the input vector
   * @param outVect is the lengthend vector
   * @param tolerence is the tolerence added to the frame vectors
   **/
  void increaseFrameVector(const HepVector3D& inVect,
			   HepVector3D& outVect,
			   double tolerance = 200.);
  
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
		   double tolerance = 200.,
		   int acdId = -1);
 
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
			   double tolerance = 200.);

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
			     double tolerance = 200.);  

  /**
   * @brief Get the list of detector elements that might have been hit by a given track
   *
   * @param point, where the track exited the ACD
   * @param theMap is the Pat. Rec. hash map
   * @return point to the list of AcdId, null for failure
   **/
  const std::list<idents::AcdId>* findNearbyElements(const HepPoint3D& point,
						     AcdPatRecMap& theMap);

  
  /**
   * @brief build the set of detector elements that might have been hit by a given track
   * @param point, where the track exited the ACD
   * @param arcTol, error on the arclenght of the point
   * @return point to the list of AcdId, null for failure
   **/
  bool buildElementSet(const HepPoint3D& p,
		       const HepVector3D& v,
		       const double& arcTol,
		       std::set<idents::AcdId>& idSet,
		       AcdPatRecMap& theMap);

}

#endif
