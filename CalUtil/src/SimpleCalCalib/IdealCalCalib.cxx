// $Header$

/** @file
    @author Z.Fewtrell
*/

// LOCAL INCLUDES
#include "CalUtil/SimpleCalCalib/IdealCalCalib.h"
// LOCAL

// GLAST

// EXTLIB

// STD
#include <vector>
#include <string>

using namespace std;

// GLAST INCLUDES
#include "xmlBase/IFile.h"

// EXTLIB 

// STD

using namespace std;

namespace CalUtil {
  const string IdealCalCalib::PEDS("PEDESTALS");
  const string IdealCalCalib::ASYM("ASYMMETRY");
  const string IdealCalCalib::THOLD_CI("THOLD_CI");
  const string IdealCalCalib::INL("INT_NONLIN");
  const string IdealCalCalib::MPD("MEV_PER_DAC");

  IdealCalCalib::IdealCalCalib() :
    asymLrgNeg(0),
    asymLrgPos(0),
    asymSmNeg(0),
    asymSmPos(0),
    asymPSNBNeg(0),
    asymNSPBNeg(0),
    asymSigPct(0),
    mpdLrg(0),
    mpdSm(0),
    mpdSigPct(0),
    inlSigPct(0) {}

  StatusCode IdealCalCalib::readCfgFile(const string &path) {
    xmlBase::IFile ifile(path.c_str());

    pedVals       = ifile.getDoubleVector(PEDS.c_str(), "VALS");
    pedCos        = ifile.getDoubleVector(PEDS.c_str(), "COS");
    pedSigs       = ifile.getDoubleVector(PEDS.c_str(), "SIG");

    asymLrgNeg    = ifile.getDouble(ASYM.c_str(), "LARGE_NEG");
    asymLrgPos    = ifile.getDouble(ASYM.c_str(), "LARGE_POS");
    asymSmNeg     = ifile.getDouble(ASYM.c_str(), "SMALL_NEG");
    asymSmPos     = ifile.getDouble(ASYM.c_str(), "SMALL_POS");
    asymPSNBNeg   = ifile.getDouble(ASYM.c_str(), "PSNB_NEG");
    asymPSNBPos   = ifile.getDouble(ASYM.c_str(), "PSNB_POS");
    asymNSPBNeg   = ifile.getDouble(ASYM.c_str(), "NSPB_NEG");
    asymNSPBPos   = ifile.getDouble(ASYM.c_str(), "NSPB_POS");
    asymSigPct    = ifile.getDouble(ASYM.c_str(), "SIG_PCT");

    mpdLrg        = ifile.getDouble(MPD.c_str(), "LARGE");
    mpdSm         = ifile.getDouble(MPD.c_str(), "SMALL");
    mpdSigPct     = ifile.getDouble(MPD.c_str(), "SIG_PCT");

    ciFLE         = ifile.getDouble(THOLD_CI.c_str(), "FLE");
    ciFHE         = ifile.getDouble(THOLD_CI.c_str(), "FHE");
    ciLAC         = ifile.getDouble(THOLD_CI.c_str(), "LAC");
    ciULD         = ifile.getDoubleVector(THOLD_CI.c_str(), "ULD");
    ciSigPct      = ifile.getDouble(THOLD_CI.c_str(), "SIG_PCT");
    ciPeds        = ifile.getDoubleVector(THOLD_CI.c_str(), "PEDS");

    inlADCPerCIDAC  = ifile.getDoubleVector(INL.c_str(), "ADC_PER_CIDAC");
    inlSigPct     = ifile.getDouble(INL.c_str(), "SIG_PCT");

    return StatusCode::SUCCESS;
  }

}
