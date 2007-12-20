//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 25 15:10:52 2007 by ROOT version 5.10/00
// from TTree ConfigCompare/ConfigCompareuration Tree
// found on file: e2e_LAT-20xGammaThruMIPCNONoPer_0.01hr.root
//////////////////////////////////////////////////////////

#ifndef ConfigCompare_h
#define ConfigCompare_h

#include <Rtypes.h>
#include <iostream>

class ConfigBase;

class ConfigCompare {
public :

  ConfigCompare(){;}
  virtual ~ConfigCompare(){;}
  
  virtual Bool_t   compare(ConfigBase& one, ConfigBase& other, 
			   Bool_t shortCompare=kFALSE, Bool_t onlyOne=kTRUE);
  
protected:
  
  Bool_t compareSingleton(const char* name, ULong64_t val1, ULong64_t val2);
  Bool_t compareSingleton(const char* name, UInt_t val1, UInt_t val2);
  Bool_t compareSingleton(const char* name, Char_t val1, Char_t val2);
  Bool_t compareV1(const char* name, UInt_t* val1, UInt_t* val2, int size, Bool_t onlyOne=kTRUE);
  Bool_t compareV1(const char* name, UShort_t* val1, UShort_t* val2, int size, Bool_t onlyOne=kTRUE);
  Bool_t compareV2(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, Bool_t onlyOne=kTRUE);
  Bool_t compareV2(const char* name, UShort_t* val1, UShort_t* val2, int size1, int size2, Bool_t onlyOne=kTRUE);
  Bool_t compareV3(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, int size3, Bool_t onlyOne=kTRUE);
  Bool_t compareV3(const char* name, ULong64_t* val1, ULong64_t* val2, int size1, int size2, int size3, Bool_t onlyOne=kTRUE);
  Bool_t compareV4(const char* name, UInt_t* val1, UInt_t* val2, int size1, int size2, int size3, int size4, Bool_t onlyOne=kTRUE);
  
};

#endif
