#ifndef CALIBDATASVC_ICALIBPATHSVC_H
#define CALIBDATASVC_ICALIBPATHSVC_H

#include "GaudiKernel/IInterface.h"
#include <string>

/**  @class ICalibPathSvc
     @brief Provide service to translate enums to strings used in
            tds path.  "Calib" is a bit of a misnomer since it's 
            anticipated that this service will also handle config items
*/
static const InterfaceID IID_ICalibPathSvc("ICalibPathSvc", 1, 0);

class ICalibPathSvc : virtual public IInterface {
public:

  static const InterfaceID& interfaceID() { return IID_ICalibPathSvc; }

  // List all the possibilities in the following enum
  typedef enum {
    CALIB = 0,
    LATC_SRC = 1
  } CalibPathType;

  typedef enum {
    Calib_TKR_HotChan       = 0,
    Calib_TKR_DeadChan,
    Calib_TKR_BadChan,
    Calib_TKR_TOTSignal,
    Calib_TKR_TOTDist,
    Calib_TKR_MIPEff,
    Calib_TKR_Splits,
    Calib_TKR_ChargeScale,
    Calib_TKR_TrgThresh,
    Calib_TKR_DataThresh,
    Calib_TKR_TowerAlign,
    Calib_TKR_InternalAlign,

    Calib_CAL_LightAtt,
    Calib_CAL_LightAsym,
    Calib_CAL_LightYield,
    Calib_CAL_Ped,
    Calib_CAL_ElecGain,
    Calib_CAL_IntNonlin,
    Calib_CAL_DiffNonlin,
    Calib_CAL_HotChan,
    Calib_CAL_DeadChan,
    Calib_CAL_MuSlope,

    Calib_CAL_MevPerDac,
    Calib_CAL_TholdCI,
    Calib_CAL_TholdMuon,
    Calib_CAL_Asym,

    Calib_ACD_Eff,
    Calib_ACD_ThreshHigh,
    Calib_ACD_ThreshVeto,
    Calib_ACD_Ped,
    Calib_ACD_ElecGain,
    Calib_ACD_Range,
    Calib_ACD_HighRange,
    Calib_ACD_CoherentNoise,
    Calib_ACD_Ribbon,
    Calib_ACD_HighPed,
    Calib_ACD_Carbon,
    Calib_ACD_VetoFit,
    Calib_ACD_CnoFit,
    Calib_ACD_PE,
               
    Calib_CalibTest1,
    
    Calib_NAS_TowerCfg,
    Calib_NAS_SAABoundary,
    Calib_NAS_LATAlignment,

    Calib_ANC_TaggerPed,
    Calib_ANC_TaggerGain,
    Calib_ANC_QdcPed,
    Calib_COUNT
  } CalibItem;

  /// Return std::string TDS path from inputs @a item (element from
  /// enumeration defined above) and optional @a flavor.  
  virtual const std::string
  getCalibPath(const CalibItem item, const std::string& flavor="") const = 0;
};

#endif
