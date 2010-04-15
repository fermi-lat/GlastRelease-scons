//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 25 15:10:52 2007 by ROOT version 5.10/00
// from TTree ConfigReport/ConfigReporturation Tree
// found on file: e2e_LAT-20xGammaThruMIPCNONoPer_0.01hr.root
//////////////////////////////////////////////////////////

#ifndef ConfigReport_h
#define ConfigReport_h

#include <Rtypes.h>
#include <iostream>

class ConfigBase;
class TList;
class TH1;

class ConfigReport {
public :

  ConfigReport() {;}
  virtual ~ConfigReport(){;}

  virtual void report(ConfigBase& c, TList& l, Option_t* options, std::ostream& os);
  
  void reportSingleton(const char* name, ULong64_t val, std::ostream& os = std::cout);   
  void reportSingleton(const char* name, UInt_t val, std::ostream& os = std::cout);
  void reportSingleton(const char* name, Char_t val, std::ostream& os = std::cout);
  
  void reportV1(const char* name, UInt_t bcast, UInt_t* val, int size, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportV1(const char* name, UShort_t bcast, UShort_t* val, int size, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportV2(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportV2(const char* name, UShort_t bcast, UShort_t* val, int size1, int size2, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportVeto(const char* name, UShort_t* veto, UShort_t* vernier,
		  std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportV3(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, int size3, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportV3(const char* name, ULong64_t bcast, ULong64_t* val, int size1, int size2, int size3, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
  void reportMask(const char* name, ULong64_t bcast, ULong64_t* val, int size1, int size2, int size3, 
		  std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);  
  void reportV4(const char* name, UInt_t bcast, UInt_t* val, int size1, int size2, int size3, int size4, 
		std::ostream& os = std::cout, TH1* h = 0, TH1* hv = 0);
 
};

#endif

