// $Header$

#ifndef CalibData_CalibModel_h
#define

/**
       @file CalibModel.h
 Definition of strings holding paths to calibration data objects
 in the TDDS (transient detector data store)
   @author : adapted from similar structures for TDS
   @author   J. Bogart
 */ 

#include <string>
#include <vector>

// For everybody except the CalibModel class implementation file,
// the variables are extern.  CalibModel actually defines them.

#if defined(_CalibData_CalibModel_cxx)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

    namespace CalibData {
      _EXTERN_ std::string   TkrBadStrips;
      //   ... more TKR calib types

      //   ... more for CAL calib types

      //   ... more for ACD calib types 

      //   ... simple type whose "data" just come from MySQL metadata row
      //       for testing
      _EXTERN_ std::string   TestMetadataInfo;

      typedef  std::vector<std::pair <std::string, CLID>> CalibPairCol;
      typedef  CalibPairCol::const_iterator PairIt;
      _EXTERN_    CalibPairCol pairs;
    }
    

#undef _EXTERN_
#endif 
