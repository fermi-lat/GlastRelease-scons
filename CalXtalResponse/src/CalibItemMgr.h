#ifndef CalibItemMgr_H
#define CalibItemMgr_H 1

// LOCAL

// GLAST
#include "idents/CalXtalId.h"
#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IService.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "TSpline.h"

// STD
#include <string>
#include <algorithm>

using namespace std;
using namespace CalDefs;
using namespace idents;

////////////////////////// GENERIC UTILITIES ///////////////////////////////////////////////

/// functional class deletes a pointer
/// fun to use w/ for_each template
///
/// I got it from here - Z.F.
/// Newsgroups: comp.lang.c++.moderated
/// From: Didier Trosset <didier-dot-tros...@acqiris.com> - Find messages by this author
/// Date: 21 Oct 2004 15:39:18 -0400
/// Subject: Re: Standard way to delete container of pointers
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

class CalibItemMgr {
 public:
  CalibItemMgr(const string &calibRootPath, int nSplineTypes=0) : 
    m_calibBase(0),
    m_calibDataSvc(0),
    m_msgSvc(0),
    m_logName(0),
    m_calibRootPath(calibRootPath),
    m_splineLists(nSplineTypes),
    m_dataProviderSvc(0),
    m_isValid(false),
    m_inUpdate(false) {}

  virtual ~CalibItemMgr() {};

  StatusCode initialize(IService &calibDataSvc, 
                        const string &flavor,
                        IMessageSvc *msgSvc,
                        const string &logName
                        );

  void invalidate() {m_isValid = false;} ///< data should be invalidated at beginning of each event.

 protected:
  virtual StatusCode updateCache();          ///< check & update serno, delete & repopulate cache if necessary.
  virtual StatusCode fillRangeBases() = 0;
  virtual StatusCode genSplines() {return StatusCode::SUCCESS;}

  /// make sure that data for one rng is sane/complete
  virtual bool validateRangeBase(const CalXtalId &xtalId, CalibData::RangeBase *rngBase) = 0; 

  virtual bool checkXtalId(const CalXtalId &xtalId) = 0;

  /// generic spline evaluation f() works for any calib_type
  StatusCode evalSpline(int calibType, LATWideIndex  idx, double x, double &y);

  /// generate new spline from 2 STL vectors & insert it into appropriate splineList
  StatusCode genSpline(int calibType, LATWideIndex idx, const string &name,
                       const vector<double> &x, const vector<double> &y);

  CalibData::CalCalibBase     *m_calibBase;
  IService                    *m_calibDataSvc;

  string                       m_calibPath;
  IMessageSvc                 *m_msgSvc;
  const string                *m_logName;
  const string                &m_calibRootPath;

  CalVec<LATWideIndex, CalibData::RangeBase* > m_rngBases;
  vector<CalVec<LATWideIndex,TSpline3* > >     m_splineLists;

 private:
  virtual void flushCache();   ///< wipe out all cached data.
  IDataProviderSvc *m_dataProviderSvc;
  string            m_flavor;
  bool              m_isValid;
  int               m_serNo;

  ///< set to true during update process to avoid recursive update checking
  bool              m_inUpdate;

};
#endif
