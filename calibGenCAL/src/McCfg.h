#ifndef McCfg_H
#define McCfg_H 1

// LOCAL INCLUDES
#include "ICfg.h"
#include "CGCUtil.h"

// GLAST INCLUDES
#include "xmlBase/IFile.h"

// EXTLIB INCLUDES

// STD INCLUDES
#include <fstream>

using namespace std;

class McCfg : ICfg {
 public:
  /// basic ctor
  McCfg() {valid = false;}
  virtual ~McCfg() {};

  /// clear all values, delete all pointers
  void clear();

  /// read in config data from config file, calculate dependent vars
  void readCfgFile(const string& cfgPath);

  /// return data valid flag.
  bool isValid() {return valid;}

  /// print summary to ostream
  void summarize();
 public:  // i know I'm not supposed to make data members public, but it's just easier this way!
  // CONFIGURABLE PARAMETERS //

  // SECTION: TEST_INFO //
  string timestamp;     ///< time of observation/measurement

  string instrument;    ///< instrument name "EM", "FM101", etc...
  vector<int> towerList; ///< list of installed towers in lat.

  // SECTION: PATHS //
  string rootFileListStr; ///< list of input root files

  string intNonlinFile; ///< input txt filename for integral non-linearity
  string dtdPath;       ///< Data descriptoin file for .xml output
  string dtdFilename;   ///< 

  string outputDir;     ///< folder for autonamed output files

  string pedFileXML;    ///< output xml filename for pedestals
  string asymFileXML;   ///< output xml filename for asymmetry calibrations
  string mpdFileXML;    ///< output xml filename for MevPerDAC calibrations
  string adc2nrgFileXML; ///< output xml filename for adc2nrg table


  string pedHistFile;   ///< output ROOT histogram file - pedestal phase
  string asymHistFile;  ///< output ROOT histogram file - asymmetry phase
  string mpdHistFile;   ///< output ROOT histogram file - MevPerDAC phase

  string pedFileTXT;    ///< output txt filename for pedestals
  string asymFileLLTXT; ///< output txt filename for asymmetry Large Diode Pos face 2 Large Diode Neg face 
  string asymFileLSTXT; ///< output txt filename for asymmetry Large Diode Pos face 2 Small Diode Neg face 
  string asymFileSLTXT; ///< output txt filename for asymmetry Small Diode Pos face 2 Large Diode Neg face 
  string asymFileSSTXT; ///< output txt filename for asymmetry Small Diode Pos face 2 Small Diode Neg face 
  string largeMPDFileTXT; ///< output txt filename for Mev per DAC Large Diode
  string smallMPDFileTXT; ///< output txt filename for Mev per DAC Small Diode

  string logfile;       ///< duplicate of stdout log
    
  // SECTION: CONSTANTS //
  double hitThresh;     ///< threshold to count a hit 
  
  double cellHorPitch;  ///< horizontal pitch between 2 cal xtals
  double cellVertPitch; ///< vertical pitch between 2 cal xtals
  double csiLength;     ///< length of one cal CsI crystal

  double maxAsymLL;     ///< used in omission of events w/ bad asymmetry logratio
  double maxAsymLS;     ///< used in omission of events w/ bad asymmetry logratio
  double maxAsymSL;     ///< used in omission of events w/ bad asymmetry logratio
  double maxAsymSS;     ///< used in omission of events w/ bad asymmetry logratio
  double minAsymLL;     ///< used in omission of events w/ bad asymmetry logratio
  double minAsymLS;     ///< used in omission of events w/ bad asymmetry logratio
  double minAsymSL;     ///< used in omission of events w/ bad asymmetry logratio
  double minAsymSS;     ///< used in omission of events w/ bad asymmetry logratio
  
  // SECTION: GENERAL //
  int nEvtRoughPed;     ///< number of events for rough pedestal calibration
  int nEvtPed;          ///< number of events for Pedestal calibration
  int nEvtAsym;         ///< number of events for Asymmetry calibration
  int nEvtMPD;          ///< number of events for MevPerDAC calibration

  bool readInPeds;      ///< skip pedestal calibration and read in previous results from .txt file
  bool readInAsym;      ///< skip Asymmetry calibration and read in previous results from .txt file
  bool skipMPD;         ///< smip MevPerDAC calibration

  bool genXML;          ///< generate xml output
  bool genTXT;          ///< generate text output
  bool genHistfiles;    ///< generate histogram output
  bool genLogfile;      ///< clone stdout stream to a logfile

  bool genOptAsymHists; ///< generate optional asymmetry histograms

  // DERIVED FROM CFG PARAMES //
  vector<string> rootFileList;
  
  /// multiplexing output stream will contain at least cout, but
  /// may also contain a logfile stream if the user requests it.
  CGCUtil::multiplexor_ostream ostrm;
  ofstream logStrm;

  /// CVS Tag string cleaned up for proper format in output XML files
  string creator;  

 private:
  string baseFilename;  ///< string shared by all autogenerated output filenames 

  // Section decription strings
  static const string TEST_INFO; ///< TEST_INFO xml IFile section name
  static const string PATHS;     ///< PATHS xml IFile section name
  static const string CONSTANTS; ///< CONSTANTS xml IFile section name
  static const string GENERAL;   ///< GENERAL xml IFile section name

  bool valid;   // set to false member data is incomplete/invalid.
};

#endif // McCfg_H
