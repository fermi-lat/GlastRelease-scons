#ifndef CalibItemMgr_H
#define CalibItemMgr_H

// LOCAL
#include "IdealCalCalib.h"

// GLAST
#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "CalUtil/CalVec.h"

// EXTLIB
#include "GaudiKernel/IService.h"

#include "TSpline.h"

// STD
#include <string>
#include <algorithm>

using namespace std;
using namespace CalUtil;
using namespace idents;

class CalCalibSvc;

/** @class CalibItemMgr
    @author Zach Fewtrell
    \brief abstract class for handling a single Cal calibration data type.

    provides the following services: 
    - TDS calibration data retrieval & indexing
    - validation period checking
    - configurable 'flavor'
    - support for optional spline functions & other local data store.
*/
class CalibItemMgr {
 public:
  CalibItemMgr(const string &calibTypePath,
               int nSplineTypes=0) : 
    m_calibTypePath(calibTypePath),
    owner(0),
    m_idealMode(false),
    m_splineLists(nSplineTypes),
    m_splineXMin(nSplineTypes),
    m_splineXMax(nSplineTypes),
    m_isValid(false),
    m_serNo(-1)
    {}

  virtual ~CalibItemMgr() {};

  StatusCode initialize(const string &flavor,
                        const CalCalibSvc &ccs);

  /// data should be invalidated at beginning of each event.
  /// just in case there is a change in validity period
  void invalidate() {m_isValid = false;} 

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
  StatusCode evalSpline(int calibType, LATWideIndex idx, float x, float &y);

  /// generate new spline from 2 STL vecs & insert into appropriate splineList
  StatusCode genSpline(int calibType, LATWideIndex idx, const string &name,
                       const vector<float> &x, const vector<float> &y);

  /// TDS path to calib data for my calib_type
  string                       m_calibTypePath;

  /// TDS path to calib data for my calib_type and path
  string                       m_calibPath;

  /// TDS location for root of my calib_type and path
  CalibData::CalCalibBase     *m_calibBase;

  /// ref to owner->CalCalibSvc object
  const CalCalibSvc           *owner;

  /// boolean if we're in ideal 'fake' mode
  bool m_idealMode;

  /// 2d vector of all (optional) splines in local data store
  vector<CalVec<LATWideIndex, TSpline3* > >    m_splineLists;
  /// min X val for each (optional) spline in local data store
  vector<CalVec<LATWideIndex, float> >         m_splineXMin;
  /// max X val for each (optional) spline in local data store
  vector<CalVec<LATWideIndex, float> >         m_splineXMax;

  /// pointers to each data member for my calib_type                                                                                                                                                       
  CalVec<LATWideIndex, CalibData::RangeBase* > m_rngBases;

  /** retrieve spec'd rangeBase object, update if necessary
      \return NULL if there is no data 
  */
  CalibData::RangeBase *getRangeBase(CalXtalId xtalId) {
    return m_calibBase->getRange(xtalId);
  }
  
 private:
  /// wipe out locally stored data (e.g. splines)
  void clearLocalStore();      
  /// calib flavor
  string            m_flavor;     
  /// validity state of CalibItemMgr data
  bool              m_isValid;    
  /// serial # for current calibration source
  int               m_serNo;      

  ///////////////////////// UTILS ///////////////////////////////////////////
  /**
     functional class deletes a pointer
     fun to use w/ for_each template
     
     I got it from here - Z.F.
     Newsgroups: comp.lang.c++.moderated
     From: Didier Trosset <didier-dot-tros...@acqiris.com>
     Date: 21 Oct 2004 15:39:18 -0400
     Subject: Re: Standard way to delete container of pointers
  */
  struct delete_ptr
  {
    template <class P>
    void operator() (P p)
    {
      delete p;
      p = 0;
    }
  };

  /// template function calls delete on all pointers stored in a STL container
  template<class T> void del_all_ptrs(T &container) {
    for_each(container.begin(),container.end(),delete_ptr());
  }
};
#endif
