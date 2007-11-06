#ifndef XtalRecTool_h
#define XtalRecTool_h
// $Header$

// LOCAL
#include "CalXtalResponse/IXtalRecTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalArray.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Digi/CalDigi.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD
#include <cstring>

/*! @class XtalRecTool
  \author Zachary Fewtrell
  \brief Simple implementation of IXtalRecTool.  

  Calculate crystal energy level and centroid from single CalDigi object.

  jobOptions:
  - CalCalibSvc (default="CalCalibSvc") - Cal Calibration source

*/

// forward declarations
class Event::CalDigi;
class ICalCalibSvc;
class IGlastDetSvc;
class Point;

class XtalRecTool : public AlgTool, 
                    virtual public IXtalRecTool 
{
public:

  /// default ctor, declares jobOptions.
  XtalRecTool( const std::string& type, 
               const std::string& name, 
               const IInterface* parent);

  /// gets needed parameters and pointers to required services
  StatusCode initialize();

  StatusCode finalize() {return StatusCode::SUCCESS;}

  StatusCode calculate(const Event::CalDigi &digi,
                       Event::CalXtalRecData &xtalRec,
                       CalUtil::CalArray<CalUtil::FaceNum, bool> &belowNoise,
                       CalUtil::CalArray<CalUtil::FaceNum, bool> &saturated,
                       INeighborXtalkTool const*const xtalkTool=0);


private:
  /// reconstruct for single readout (one adc range on each crytsal face)
  /// \param belowNoise, set true on either face if signal is below noise threshold
  /// \param saturated, set true on either face if HEX1 adc range is saturated
  /// \return 0 on failure. pointer to new object on success (caller is responsible for deallocation of returned CalRangeRecData object.
  Event::CalXtalRecData::CalRangeRecData *createRangeRecon(const CalUtil::XtalIdx xtalIdx,
                                                           const Event::CalDigi::CalXtalReadout &ro,
                                                           CalUtil::CalArray<CalUtil::FaceNum, bool> &belowNoise,
                                                           CalUtil::CalArray<CalUtil::FaceNum, bool> &saturated,
                                                           INeighborXtalkTool const*const xtalkTool) const;

  /** \brief convert large diode CIDAC scale to small diode CIDAC scale
      for given xtal face & pos 

      
      \param face which xtal face to process (uses current xtal from
      m_dat struct)

      \param pos longitudinal position along xtal (mm)
      \param largeCIDAC input large diode CIDAC value
      \param smallCIDAC output small diode CIDAC value
  */
  StatusCode largeCIDAC2Small(const CalUtil::XtalIdx xtalIdx,
                              CalUtil::FaceNum face, 
                              const float pos, 
                              const float largeCIDAC, 
                              float &smallCIDAC) const;


  /** \brief convert scalar position in mm from xtal center (longitudinal
      to 3d vector position
  */
  void pos2Point(const CalUtil::XtalIdx xtalIdx,
                 const float pos, 
                 Point &pXtal) const;

  /// name of CalCalibSvc to use for calib constants.
  StringProperty m_calCalibSvcName;                         
  /// pointer to CalCalibSvc object.
  ICalCalibSvc *m_calCalibSvc;  

  /// pointer to the Glast Detector Service
  IGlastDetSvc* m_detSvc; 

  /// length of CsI xtal in mm
  float m_CsILength;
  /// the value of fLATObjects field, defining LAT towers 
  int m_eLATTowers; 
  
  /// the value of fTowerObject field, defining calorimeter module 
  int m_eTowerCAL;
  /// the value of fCellCmp field defining CsI crystal
  int m_eXtal;      
  /// number of geometric segments per Xtal
  int m_nCsISeg; 

};


#endif
