#ifndef IAcdCalibSvc_H
#define IAcdCalibSvc_H
//  $Header$

// LOCAL INCLUDES
#include "AcdUtil/AcdCalib.h"

// GLAST INCLUDES
#include "idents/AcdId.h" 

// EXTLIB INCLUDES
#include "GaudiKernel/IInterface.h"

// STD INCLUDES
//#include <vector>

// Forward declares
namespace CalibData {
  class AcdPed;
  class AcdGain;
  class AcdVeto;
  class AcdCno;
  class AcdRange;
  class AcdHighRange;
  class AcdCoherentNoise;
  class AcdRibbon;
}

class AcdCalibMgr;


// Declaration of the interface ID ( interface id, major version,
// minor version)
static const InterfaceID IID_IAcdCalibSvc("IAcdCalibSvc", 2, 0);

/** 
 * @class AcdUtil::IAcdCalibSvc
 * 
 * @brief Abstract interface for provision of GLAST LAT ACD calibration constants
 *
 * Functions are provided to access every type of ACD calibration 
 * by Channel and PMT:
 *   - CalibData::AcdPed           Pedestal (in PHA counts)
 *   - CalibData::AcdGain          Gain (aka MIP peak in PHA above pedestal)
 *   - CalibData::AcdVeto          Veto threshold (in PHA counts)
 *   - CalibData::AcdCno           CNO threshold (in PHA counts)
 *   - CalibData::AcdRange         Range crossover (in PHA counts)
 *   - CalibData::AcdHighRange     High-Range PHA -> MIPs calibration (pedestal, slope, saturation)
 *   - CalibData::AcdCoherentNoise Coherent Noise effect ( [0] * exp([1]*x) * sin( [2]*x + [3] )
 *   - CalibData::AcdRibbon        Ribbon light attenuation
 *   
 * All functions return Success or Failure and fill the provided parameter.
 * Type checking is included.
 *
 * @author Eric Charles (from Zach Fewtrell's CalCalibSvc)
 *
 */

namespace AcdUtil {
  
  class IAcdCalibSvc : virtual public IInterface {
  public:
    static const InterfaceID& interfaceID() { return IID_IAcdCalibSvc; }

     /** \brief get a calibration for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param calib a pointer to the relevent calibration data
    */
    virtual StatusCode getCalibMgr(AcdCalibData::CALTYPE type,
				   AcdCalibMgr*& calib) = 0;
    

    /** \brief get pedestals for given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param pedestal a pointer to the relevent pedestal data
    */
    virtual StatusCode getPedestal(idents::AcdId id, unsigned pmt,
				   CalibData::AcdPed*& pedestal);
			      
    /** \brief get mip peak for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param pedestal a pointer to the relevent gain data
    */
    virtual StatusCode getMipPeak(idents::AcdId id, unsigned pmt,
				  CalibData::AcdGain*& mipPeak);

    /** \brief get veto threshold for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param veto a pointer to the relevent threshold data
    */
    virtual StatusCode getVeto(idents::AcdId id, unsigned pmt,
		     CalibData::AcdVeto*& veto);

    /** \brief get cno threshold for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param cno a pointer to the relevent threshold data
    */
    virtual StatusCode getCno(idents::AcdId id, unsigned pmt,
			      CalibData::AcdCno*& cno);
    
    /** \brief get range crossover for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param range a pointer to the relevent data
    */
    virtual StatusCode getRange(idents::AcdId id, unsigned pmt,
				CalibData::AcdRange*& range);

    /** \brief get high range calibration for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param highRange a pointer to the relevent high range data
    */
    virtual StatusCode getHighRange(idents::AcdId id, unsigned pmt,
				   CalibData::AcdHighRange*& highRange);
    
    /** \brief get coherent noise calibration for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param calib a pointer to the relevent high range data
    */
    virtual StatusCode getCoherentNoise(idents::AcdId id, unsigned pmt,
					CalibData::AcdCoherentNoise*& noiseCalib);

    /** \brief get ribbon calibration for a given channel
	\param id  the tile or ribbon id
	\param pmt A(0) or B(1) pmt
	\param calib a pointer to the relevent high range data
    */
    virtual StatusCode getRibbon(idents::AcdId id, unsigned pmt,
				 CalibData::AcdRibbon*& ribbon);
  };
};

#endif // IAcdCalibSvc_H
