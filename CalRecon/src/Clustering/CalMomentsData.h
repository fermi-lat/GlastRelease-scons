/**
 * @class CalMomentsData
 *
 * @brief Utility data object for the moments analysis which  attempts to make the
 * class independent of the actual Cal data objects used.
 *
 * @author Tracy Usher, Luca Baldini.
 *
 */

/// Maximum distance to declare the hit position on the xtal edge.
#define DIST_LONG_POS_INVALID  0.01
/// Maximum distance to declare the fit position near the xtal edge. 
#define DIST_FIT_POS_NEAR_EDGE 25.0

#include "geometry/Ray.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalFitParams.h"
#include "CalRecon/ICalReconSvc.h"

#include <vector>

class CalMomentsData
{
 public:

  /// Status word bits organized like:
  /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
  ///         <   Not  used   > <Analysis status> <Fit  correction>  <Xtal properties>
  /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
  ///         <                               Not used                               >
  enum StatusBits {
    /// All empty!
    ZERO              = 0x0,
    /// Longitudinal position right on the xtal edge (i.e. it was outside the xtal).
    LONG_POS_INVALID  = 0x00000001,
    /// Saturated.
    SATURATED         = 0x00000002,
    /// Longitudinal position from the fit available.
    FIT_POS_AVAILABLE = 0x00000010,
    /// The fit axis is passing close to one of the xtal edges.
    FIT_POS_NEAR_EDGE = 0x00000020,
    /// Longitudinal position from the fit out of the xtal
    /// (set right on one of the xtal edges when applying the correction).
    FIT_POS_INVALID   = 0x00000040,
    /// Use the corrected longitudinal position (i.e. the value from the fit).
    USE_FIT_POS       = 0x00000080,
    /// Do not use it in the moments analysis.
    MASKED            = 0x00000100,
  };

  /// Constructor from all the parameters.
  CalMomentsData(const Point& position, const double weight, int tower, int layer, int column);

  /// Convenience constructor from a Event::CalXtalRecData* recData object.
  CalMomentsData(Event::CalXtalRecData* recData);
    
  /// Destructor.
  ~CalMomentsData() {}
    
  /// Provides access to data.
  inline unsigned int getStatusBits()        const { return m_statusBits; }
  inline const Point& getBasePosition()      const { return m_basePosition; }
  inline const Point& getCorrPosition()      const { return m_corrPosition; }
  inline double getWeight()                  const { return m_weight; }
  inline int getTower()                      const { return m_tower; }
  inline int getLayer()                      const { return m_layer; }
  inline int getColumn()                     const { return m_column; }
  inline double getDistToAxis()              const { return m_distToAxis; }
  inline double getCoordAlongAxis()          const { return m_coordAlongAxis; }

  /// Provides set functions.
  inline void setBasePosition(const Point& pos)    { m_basePosition = pos; }
  inline void setCorrPosition(const Point& pos)    { m_corrPosition = pos; }
  inline void setWeight(double weight)             { m_weight = weight; }
  inline void setTower(int tower)                  { m_tower = tower; }
  inline void setLayer(int layer)                  { m_layer = layer; }
  inline void setColumn(int column)                { m_column = column; }
  inline void setDistToAxis(double dist)           { m_distToAxis = dist; }
  inline void setCoordAlongAxis(double coord)      { m_coordAlongAxis = coord; }

  /// Methods to check wheter the underlying xtals lies along the x or y axis.
  bool isx()                                 const { return (m_layer % 2 == 0); }
  bool isy()                                 const { return (m_layer % 2 != 0); }

  /// Manipulate status bits.
  inline void setStatusBit(StatusBits bit)         { m_statusBits |=  bit; }
  inline void clearStatusBit(StatusBits bit)       { m_statusBits &= ~bit; }
  inline bool checkStatusBit(StatusBits bit) const { return ( (m_statusBits & bit) != ZERO ); }
  inline void mask()                               { setStatusBit(MASKED); }
  inline void unmask()                             { clearStatusBit(MASKED); }
  inline bool masked()                       const { return checkStatusBit(MASKED); }
 
  /// Apply the fit longitudinal correction, i.e. set the longitudinal position of the crystal
  /// to the longitudinal position of the extrapolation of the cal direction using only the
  /// transverse information from the fit.
  void applyFitCorrection(const Event::CalFitParams fitParams, const ICalReconSvc* calReconSvc);

  /// If the fit correction is available and the corresponding position falls within
  /// the xtal, set the USE_FIT_POS status bit in order to use it.
  void enableFitCorrection();

  /// Return the "best position" (i.e. use the fit longitudinal position if the
  /// USE_FIT_POS status bit is set).
  const Point& getPosition() const;

  /// Get the value of the fit longitudinal correction (in mm), i.e. the
  /// difference between the fit position and the xtal position along the xtal.
  double getFitCorrAmount() const; 

  /// Determine distance to given axis.
  /// Note that this sets a class member, as well.
  double calcDistToAxis(const Point& centroid, const Vector& axis);

  /// Determine the coordinate along a given axis.
  /// Note that this sets a class member, as well.
  double calcCoordAlongAxis(const Point& centroid, const Vector& axis);

  /// Set the reference centroid/axis for the CAL hit.
  /// For historical reasons the previous two functions calls (calcDistToAxis() and
  /// calcCoordAlongAxis()) are used throughout the code, though it would be probably
  /// neater to use this method, followed by the proper get methods.
  void setReferenceAxis(const Point& centroid, const Vector& axis);

  /// Define how to sort.
  const bool operator<(const CalMomentsData& right)  const
  {return m_distToAxis < right.getDistToAxis();}

  /// Std output facility.
  std::ostream& fillStream(std::ostream& s) const;
  friend std::ostream& operator<< (std::ostream& s, const CalMomentsData& obj)
  {
    return obj.fillStream(s);
  }

  
 private:
  /// The status bits.
  unsigned int m_statusBits;
  /// The position of this data point.
  Point m_basePosition;
  /// The position after the fit correction of the coordinate along the xtal.
  Point m_corrPosition;
  /// A weight to assign to the point in the moments calculation.
  double m_weight;
  /// The CAL tower this hit belongs to.
  int m_tower;
  /// The CAL layer this hit belongs to.
  int m_layer;
  /// The CAL column this hit belongs to.
  int m_column;
  /// The distance from the "axis" of this point.
  double m_distToAxis;
  /// The position along the "axis" of this point (with sign, used to calculate the skewness).
  double m_coordAlongAxis;
};


typedef std::vector<CalMomentsData> CalMomentsDataVec;

