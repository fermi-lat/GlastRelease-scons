#ifndef CalibItemMgr_H
#define CalibItemMgr_H
// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL
#include "CalCalibShared.h"

// GLAST
#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"

// EXTLIB
#include "GaudiKernel/IService.h"

#include "TSpline.h"

// STD
#include <string>
#include <vector>

class CalCalibSvc;

/** @class CalibItemMgr
    @author Z.Fewtrell
    \brief abstract class for handling a single Cal calibration data type.

    provides the following services: 
    - TDS calibration data retrieval & indexing
    - validation period checking
    - configurable 'flavor'
    - support for optional spline functions & other local data store.
*/
class CalibItemMgr {
public:
  /// \param calibItem TDS path to calibration data to be managed
  /// \param ccsShared CalCalibSvc shared resources
  /// \param indexSize number of channels to manage (array size)
  /// \param nSplineTypes number of splines per channel
  CalibItemMgr(const ICalibPathSvc::CalibItem calibItem,
               const CalCalibShared &ccsShared,
               const size_t indexSize,
               const size_t nSplineTypes=0) : 
    m_calibItem(calibItem),
    m_calibBase(0),
    m_ccsShared(ccsShared),
    m_idealMode(false),
    /// size of inner array is known only by base class
    m_splineLists(nSplineTypes,CalUtil::CalVec<CalUtil::LATWideIndex, TSpline3*>(indexSize)), 
    m_splineXMin(nSplineTypes,CalUtil::CalVec<CalUtil::LATWideIndex, float>(indexSize)),
    m_splineXMax(nSplineTypes,CalUtil::CalVec<CalUtil::LATWideIndex, float>(indexSize)),
    m_splineYMin(nSplineTypes,CalUtil::CalVec<CalUtil::LATWideIndex, float>(indexSize)),
    m_splineYMax(nSplineTypes,CalUtil::CalVec<CalUtil::LATWideIndex, float>(indexSize)),
    m_rngBases(indexSize),
    m_isValid(false),
    m_serNo(SERNO_NODATA)
  {}

  virtual ~CalibItemMgr() {};

  /// \brief initialization code which must be done during Gaudi init() period
  /// 
  /// some data is not available at construction time
  virtual StatusCode initialize(const std::string &flavor);

  /// data should be invalidated at beginning of each event.
  /// just in case there is a change in validity period
  void invalidate() {m_isValid = false;} 

  /// return serial # for current calib data
  /// return calibration data serial number from CalibSvc, or SERNO_IDEAL, SERNO_NODATA
  int getSerNo(); // {return m_serNo;}

protected:
  /** \brief check calib validity period, (re)build local store if necessary

  needs to be called once per event (_before_ processing calibration data ;).
  Subsequent calls in same event will return immediately.

  */
  virtual StatusCode updateCalib();          

  /// load ideal (fake) calibration vals for my calib_type if db is down
  virtual StatusCode loadIdealVals() = 0;

  /// generate full set of spline f()'s or other local data (if
  /// applicable for calib_type)
  virtual StatusCode genLocalStore() = 0;

  /// generic spline evaluation f() works for any calib_type
  StatusCode evalSpline(const int calibType, 
                        const CalUtil::LATWideIndex idx, 
                        const float x, 
                        float &y);

  /// generate new spline from 2 STL vecs & insert into appropriate splineList
  StatusCode genSpline(const int calibType, 
                       const CalUtil::LATWideIndex idx, 
                       const std::string &name,
                       const std::vector<float> &x, 
                       const std::vector<float> &y);
  
  /// used for generate calib path
  const ICalibPathSvc::CalibItem m_calibItem;

  /// TDS path to calib data for my calib_type and path
  std::string m_calibPath;

  /// TDS location for root of my calib_type and path
  CalibData::CalCalibBase *m_calibBase;

  /// ref to data shared by all classes used by CalibDataSvc
  const CalCalibShared &m_ccsShared;

  /// boolean if we're in ideal 'fake' mode
  bool m_idealMode;

  /// 2d vector of all (optional) splines in local data store
  std::vector<CalUtil::CalVec<CalUtil::LATWideIndex, TSpline3* > >    m_splineLists;
  /// min X val for each (optional) spline in local data store
  std::vector<CalUtil::CalVec<CalUtil::LATWideIndex, float> >         m_splineXMin;
  /// max X val for each (optional) spline in local data store
  std::vector<CalUtil::CalVec<CalUtil::LATWideIndex, float> >         m_splineXMax;
  /// min Y val for each (optional) spline in local data store
  std::vector<CalUtil::CalVec<CalUtil::LATWideIndex, float> >         m_splineYMin;
  /// max Y val for each (optional) spline in local data store
  std::vector<CalUtil::CalVec<CalUtil::LATWideIndex, float> >         m_splineYMax;

  /// pointers to each data member for my calib_type                                                                                                                                                       
  CalUtil::CalVec<CalUtil::LATWideIndex, const CalibData::RangeBase* > m_rngBases;

  /** retrieve spec'd rangeBase object, update if necessary
      \return NULL if there is no data 
  */
  const CalibData::RangeBase *getRangeBase(const idents::CalXtalId xtalId) {
    return m_calibBase->getRange(xtalId);
  }
    
  /// m_serNo has this value when data has not yet been loaded
  static const int SERNO_NODATA = -1;
  /// m_serNo has this value when ideal mode data has been loaded
  static const int SERNO_IDEAL  = 1;
  
private:
  /// wipe out locally stored data (e.g. splines)
  void clearLocalStore();      

  /// calib flavor
  std::string            m_flavor;     
  /// validity state of CalibItemMgr data
  bool              m_isValid;    

  /// serial # for current calibration source
  /// may be set to SERNO_IDEAL, or SERNO_NODATA
  int               m_serNo;      
  

};
#endif
