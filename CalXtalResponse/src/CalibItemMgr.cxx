// LOCAL
#include "CalibItemMgr.h"
#include "CalCalibSvc.h"

// GLAST
#include "CalibData/DacCol.h"
#include "CalibData/Cal/Xpos.h"

// EXTLIB

// STD
#include <sstream>
#include <algorithm>

using namespace CalDefs;

using namespace std;
///////////////// GENERIC UTILITIES //////////////////////////////////
template<typename T> const T& max_val(const vector<T> &vec) {
  return *(max_element(vec.begin(),vec.end()));
}

template<typename T> const T& min_val(const vector<T> &vec) {
  return *(min_element(vec.begin(),vec.end()));
}
//////////////////////////////////////////////////////////////////////

StatusCode CalibItemMgr::initialize(const string &flavor, const CalCalibSvc &ccs) {
  StatusCode sc;

  owner = &ccs;
  
  m_flavor = flavor;

  m_calibPath = m_calibTypePath + "/" + flavor;
    
  sc = loadIdealVals();
  if (sc.isFailure()) return sc;

  //-- IDEAL MODE --//
  if (m_flavor == "ideal") {
     m_idealMode = true;

     // splines will not be automatically 
     // generated b/c we won't actually look
     // up any data.
     // it's a hack, but screw it.
     sc = genSplines();
     if (sc.isFailure()) return sc;
  }

  return StatusCode::SUCCESS;
}

StatusCode CalibItemMgr::updateCache() {
  StatusCode sc;

  // if event is already validated return quickly
  // also we don't need validation if we're in ideal mode
  if (m_isValid || m_idealMode) return StatusCode::SUCCESS;

  // enable 'inUpdate' flag for duration of this function
  m_inUpdate = true;

  /////////////////////////////////
  //-- CHECK TDS DATA VALIDITY --//
  /////////////////////////////////

  // Retrieve pointer to Gain tree from TDS
  // usually this f() should return immediately
  // if it fails then we have no valid calib data
  // for the current event.
  DataObject *pObject;
  sc = owner->m_dataProviderSvc->retrieveObject(m_calibPath, pObject);
  if (!sc.isFailure())
    m_calibBase = (CalibData::CalCalibBase *)(pObject);
  else {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    
    // else return error (can't find calib)
    msglog << MSG::ERROR << "Unable to retrieve " 
           << m_calibPath << " from calib db" << endreq;
    return sc;  
  }

  ///////////////////////////////////////
  //-- CHECK IF TDS DATA HAS CHANGED --//
  ///////////////////////////////////////

  // check serial # to see if we're still valid.
  int curSerNo = m_calibBase->getSerNo();
  if (curSerNo != m_serNo) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::INFO << "Updating " << m_calibPath << endreq;
    m_serNo = curSerNo;
    flushCache();
        
    // load up all RangeBases for new calib set
    sc = fillRangeBases();
    if (sc.isFailure()) return sc;
    
    // build all splines (if virtual f() is supplied)
    sc = genSplines();
    if (sc.isFailure()) return sc;
  }

  // reset 'in update' flag
  m_inUpdate = false;
  m_isValid = true;

  return StatusCode::SUCCESS;
}

StatusCode CalibItemMgr::evalSpline(int calibType, LATWideIndex idx, 
                                    double x, double &y) {
  // make sure we have valid calib data for this event.
  if (!m_inUpdate) {
    StatusCode sc;
    if ((sc = updateCache()).isFailure()) return sc;
  }

  // bounds check input to spline function to avoid
  // weird behavior
  x = max<double>(m_splineXMin[calibType][idx],x);
  x = min<double>(m_splineXMax[calibType][idx],x);

  y = m_splineLists[calibType][idx]->Eval(x);
  
  if (owner->m_superVerbose) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::VERBOSE << "Evaluating spline: "
           << m_splineLists[calibType][idx]->GetName()
           << " X=" << x
           << " Y=" << y
           << endreq;
  }
  
  return StatusCode::SUCCESS;
}

StatusCode CalibItemMgr::genSpline(int calibType, LATWideIndex idx, const string &name, 
                                   const vector<double> &x, const vector<double> &y) {
  int n = min(x.size(),y.size());

  // create tmp arrays for TSpline ctor
  double *xp = new double[n];
  double *yp = new double[n];

  // copy vector data into temp arrays
  copy(x.begin(),x.begin()+n,xp);
  copy(y.begin(),y.begin()+n,yp);

  TSpline3 *mySpline = new TSpline3(name.c_str(),
                                    xp,yp,n);
  mySpline->SetName(name.c_str());

  // put spline in list
  m_splineLists[calibType][idx] = mySpline;
  // populate x-axis boundaries
  m_splineXMin[calibType][idx] = xp[0];
  m_splineXMax[calibType][idx] = xp[n-1];

  if (owner->m_superVerbose) {
    // create MsgStream only when needed for performance
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::VERBOSE << "Generated spline " << name.c_str() 
           << " t="  << calibType 
           << " i="  << idx.getInt()
           << " nx=" << n
           << " ("   << m_splineXMin[calibType][idx] 
           << "->"   << m_splineXMax[calibType][idx] << ")"
           << " ny=" << n
           << " ("   << yp[0]
           << "->"   << yp[n-1]
           << ")"    << endreq;

    // clear heap variables
    delete xp;
    delete yp;
  }

  return StatusCode::SUCCESS;
}

/// Template function fills any STL type container with zero values
template <class T> static void fill_zero(T &container) {
  fill(container.begin(), container.end(), 0);
}

void CalibItemMgr::flushCache() {   
  // Gaudikeeps track of RangeBase objects, 
  // so i only need to delete the pointers
  m_rngBases.clear();

  // m_splineLists is the 'owner' of the splines, so i need to delete
  // the objects themselves as well as the pointers.
  for (unsigned i = 0; i < m_splineLists.size(); i++) {
    del_all_ptrs(m_splineLists[i]);
    fill_zero(m_splineXMin[i]);
    fill_zero(m_splineXMax[i]);
  }
}
