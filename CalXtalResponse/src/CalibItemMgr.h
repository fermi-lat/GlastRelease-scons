#ifndef CalibItemMgr_H
#define CalibItemMgr_H 1

// LOCAL
#include "IdealCalCalib.h"

// GLAST
#include "idents/CalXtalId.h"
#include "CalibData/RangeBase.h"
#include "CalibData/Cal/CalCalibBase.h"
#include "CalUtil/CalDefs.h"

// EXTLIB
#include "GaudiKernel/IService.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include "TSpline.h"

// STD
#include <string>
#include <algorithm>

using namespace std;
using namespace CalDefs;
using namespace idents;

////////////////////////// GENERIC UTILITIES ///////////////////////////////////

class CalibItemMgr {
 public:
  CalibItemMgr(const string &calibRootPath, 
               const IdealCalCalib &idealCalib,
               int nSplineTypes=0
               ) : 
    m_calibBase(0),
    m_calibDataSvc(0),
    m_msgSvc(0),
    m_logName(0),
    m_calibRootPath(calibRootPath),
    m_splineLists(nSplineTypes),
    m_dataProviderSvc(0),
    m_isValid(false),
    m_idealMode(false),
    m_idealCalib(idealCalib),
    m_inUpdate(false) {}

  virtual ~CalibItemMgr() {};

  StatusCode initialize(IService &calibDataSvc, 
                        const string &flavor,
                        IMessageSvc *msgSvc,
                        const string &logName
                        );

  /// data should be invalidated at beginning of each event.
  void invalidate() {m_isValid = false;} 

 protected:
  /// check & update serno, delete & repopulate cache if necessary.
  virtual StatusCode updateCache();          
  /// load ideal (fake) calibration vals for my calib_type if db is down
  virtual StatusCode loadIdealVals() = 0;
  /// retireve pointers to full set of calib_data in TDS
  virtual StatusCode fillRangeBases() = 0;
  /// generate full set of splines (if applicable for calib_type)
  virtual StatusCode genSplines() {return StatusCode::SUCCESS;}   

  /// make sure that data for one rng is sane/complete
  virtual bool validateRangeBase(const CalXtalId &xtalId, 
                                 CalibData::RangeBase *rngBase) = 0; 

  /// check that xtalId has optional range & face information (if applicable for calib_type)
  virtual bool checkXtalId(const CalXtalId &xtalId) = 0;

  /// generic spline evaluation f() works for any calib_type
  StatusCode evalSpline(int calibType, LATWideIndex  idx, double x, double &y);

  /// generate new spline from 2 STL vecs & insert into appropriate splineList
  StatusCode genSpline(int calibType, LATWideIndex idx, const string &name,
                       const vector<double> &x, const vector<double> &y);

  /// TDS location for root of all used calib data
  CalibData::CalCalibBase     *m_calibBase;
  /// service used to retrieve calib data
  IService                    *m_calibDataSvc;

  /// TDS path to calib data for my calib_type
  string                       m_calibPath;
  /// Pointer to message log service
  IMessageSvc                 *m_msgSvc;
  /// name to use in message log
  const string                *m_logName;
  /// TDS path to root for all calib_types
  const string                &m_calibRootPath;

  /// pointers to each data member for my calib_type
  CalVec<LATWideIndex, CalibData::RangeBase* > m_rngBases;
  /// 2d vector of all splines used for my calib_type
  vector<CalVec<LATWideIndex,TSpline3* > >     m_splineLists;

  /// set to 'true' when we are using 'ideal' mode data
  bool                         m_idealMode;
  
  /// reference to 'ideal' mode calibration values
  const IdealCalCalib         &m_idealCalib;

 private:
  virtual void flushCache();      ///< wipe out all cached data.
  string            m_flavor;     ///< calib flavor
  bool              m_isValid;    ///< validity state of CalibItemMgr data
  int               m_serNo;      ///< serial # for current calibration source

  /// pointer to IDataProviderSvc interface of CalibDataSvc
  IDataProviderSvc *m_dataProviderSvc;

  /// set to true during update process to avoid recursive checking
  bool              m_inUpdate;


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
