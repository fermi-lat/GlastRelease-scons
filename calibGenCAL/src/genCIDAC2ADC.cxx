// $Header$

/** @file Gen CIDAC2ADC calibrations from singlex16 charge injection event files
    @author Zachary Fewtrell
*/

// LOCAL INCLUDES
#include "lib/Algs/IntNonlinAlg.h"
#include "lib/CalibDataTypes/CIDAC2ADC.h"
#include "lib/Util/CGCUtil.h"
#include "lib/Util/CfgMgr.h"

// GLAST INCLUDES
#include "CalUtil/CalDefs.h"

// EXTLIB INCLUDES

// STD INCLUDES
#include <iostream>
#include <string>
#include <fstream>

using namespace std;
using namespace CalUtil;
using namespace CGCUtil;
using namespace CfgMgr;

class AppCfg {
public:
  AppCfg(const int argc,
         const char **argv) :
    cmdParser(path_remove_ext(__FILE__)),
    rootFileHE("rootFileHE",
               'h',
               "input HE DIODE singlex16 digi root file",
               ""),
    rootFileLE("rootFileLE",
               'h',
               "input LE DIODE singlex16 digi root file",
               ""),
    columnMode("bcastMode",
               'b',
               "singlex16 pulses 12 columns individually"),
    outputBasename("outputBasename",
                   "all output files will use this basename + some_ext",
                   "")
  {
    cmdParser.registerArg(outputBasename);

    cmdParser.registerVar(rootFileHE);
    cmdParser.registerVar(rootFileLE);


    try {
      cmdParser.parseCmdLine(argc, argv);
    } catch (exception &e) {
      cout << e.what() << endl;
      cmdParser.printUsage();
      throw e;
    }

  }
  /// construct new parser
  CmdLineParser cmdParser;
  
  CmdOptVar<string> rootFileHE;
  CmdOptVar<string> rootFileLE;

  CmdSwitch columnMode;
  
  CmdArg<string> outputBasename;

};

int main(int argc,
         const char **argv) {
  // libCalibGenCAL will throw runtime_error
  try {
    AppCfg cfg(argc,argv);

    // i can process 1 or 2 files, but not none
    if (cfg.rootFileLE.getVal().length() == 0 && cfg.rootFileHE.getVal().length() == 0) {
      cout << __FILE__ << ": no input files specified." << endl;
      return -1;
    }

    //-- SETUP LOG FILE --//
    /// multiplexing output streams
    /// simultaneously to cout and to logfile
    LogStream::addStream(cout);
    // generate logfile name
    string logfile(cfg.outputBasename.getVal() + ".log.txt");
    ofstream tmpStrm(logfile.c_str());

    LogStream::addStream(tmpStrm);

    //-- LOG SOFTWARE VERSION INFO --//
    output_env_banner(LogStream::get());

    // txt output filename
    string       outputTXTFile(cfg.outputBasename.getVal() + ".txt");

    CIDAC2ADC    adcMeans;
    CIDAC2ADC    cidac2adc;
    IntNonlinAlg inlAlg;

    string       adcMeanFile(cfg.outputBasename.getVal() + ".adcmean.txt");

    if (cfg.rootFileLE.getVal().length()) {
      LogStream::get() << __FILE__ << ": reading LE calibGen event file: " << cfg.rootFileLE.getVal() << endl;
      inlAlg.readRootData(cfg.rootFileLE.getVal(), adcMeans, LRG_DIODE, !cfg.columnMode.getVal());
    }

    if (cfg.rootFileHE.getVal().length()) {
      LogStream::get() << __FILE__ << ": reading HE calibGen event file: " << cfg.rootFileHE.getVal() << endl;
      inlAlg.readRootData(cfg.rootFileHE.getVal(), adcMeans, SM_DIODE,  !cfg.columnMode.getVal());
    }

    LogStream::get() << __FILE__ << ": saving adc means to txt file: "
                     << adcMeanFile << endl;
    adcMeans.writeTXT(adcMeanFile);

    LogStream::get() << __FILE__ << ": generating smoothed spline points: " << endl;
    inlAlg.genSplinePts(adcMeans, cidac2adc);

    LogStream::get() << __FILE__ << ": writing smoothed spline points: " << outputTXTFile << endl;
    cidac2adc.writeTXT(outputTXTFile);
  } catch (exception &e) {
    cout << __FILE__ << ": exception thrown: " << e.what() << endl;
  }

  return 0;
}

