#ifndef XtalRecTool_h
#define XtalRecTool_h
// $Header$
/** @file
    @author Z.Fewtrell
*/

// LOCAL
#include "CalXtalResponse/IXtalRecTool.h"

// GLAST
#include "CalUtil/CalDefs.h"
#include "CalUtil/CalVec.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Digi/CalDigi.h"

// EXTLIB
#include "GaudiKernel/AlgTool.h"

// STD

/*! @class XtalRecTool
  \author Z.Fewtrell
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

  /// generate hit reconstruction for single Cal Crystal
  StatusCode calculate(const Event::CalDigi &digi,
                       Event::CalXtalRecData &xtalRec,
                       CalUtil::CalVec<CalUtil::FaceNum, bool> &belowNoise,
                       CalUtil::CalVec<CalUtil::FaceNum, bool> &saturated,
                       CalUtil::CalVec<CalUtil::FaceNum, bool> &acdSaturated,
                       INeighborXtalkTool const*const xtalkTool=0);

  /// Code to resolve ambiguity region of light asymmetry response
  StatusCode ambiguity(Event::CalXtalRecCol* pxtalrecs);

private:
  /// reconstruct for single readout (one adc range on each crytsal face)
  /// \param belowNoise, set true on either face if signal is below noise threshold
  /// \param saturated, set true on either face if HEX1 adc range is saturated
  /// \return 0 on failure. pointer to new object on success (caller is responsible for deallocation of returned CalRangeRecData object.
  Event::CalXtalRecData::CalRangeRecData *createRangeRecon(const CalUtil::XtalIdx xtalIdx,
                                                           const Event::CalDigi::CalXtalReadout &ro,
                                                           CalUtil::CalVec<CalUtil::FaceNum, bool> &belowNoise,
                                                           CalUtil::CalVec<CalUtil::FaceNum, bool> &saturated,
                                                           CalUtil::CalVec<CalUtil::FaceNum, bool> &acdSaturated,
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

  /** \brief convert 3d vector position to scalar long. position in mm from xtal center 
  */
  void point2Pos(const idents::CalXtalId xtalId,
                 float &pos, 
                 const Point pXtal) const;

  /** \brief convert 3d vector position to scalar transverse position in mm from Tower center 
  */
  void point2PosTrans(const idents::CalXtalId xtalId,
                 float &pos, 
                 const Point pXtal) const;

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
