// $Header$
/** @file test_meta.cxx
    Sample program to exercise calibration metadata database services
*/

#include <string>
#include <iostream>
#include "calibUtil/Metadata.h"

void gotIt(unsigned int ser, calibUtil::Metadata::eDataFmt dataFmt,
           const std::string& fmtVersion, const std::string& dataIdent);

calibUtil::Metadata::eRet lookup(calibUtil::Metadata::eCalibType ctype,
                                 const calibUtil::Timestamp& ts,
                                 unsigned int levelMask,
                                 calibUtil::Metadata::eInstrument inst);

int main(int argc, char* argv[]) {
  using calibUtil::Metadata;
  using calibUtil::Timestamp;

  Metadata  meta;
  Timestamp t_ok("2001-11-10 08:00");
  Timestamp t_none("2000-09-08 10:00");

  Metadata::eRet ret = lookup(Metadata::CTYPE_TKRBadChan, t_ok, 
                              Metadata::LEVELDev, Metadata::INSTBtem);

  ret = lookup(Metadata::CTYPE_TKRBadChan, t_ok, 
               Metadata::LEVELProd | Metadata::LEVELDev, Metadata::INSTBtem);

  ret = lookup(Metadata::CTYPE_ACDEff, t_ok, 
               Metadata::LEVELProd | Metadata::LEVELDev, 
               Metadata::INSTBtem);
  ret = lookup(Metadata::CTYPE_TKRBadChan, t_none, Metadata::LEVELDev,
               Metadata::INSTBtem);

  // Try to insert a record
  ret = meta.openRecord(Metadata::INSTEm, Metadata::CTYPE_TKRBadChan,
                   Metadata::FMTXml, "1.0",
                   "$CALIBUTILROOT/xml/test/testHot-2002-05-02.xml",
                   Metadata::CMPLOk);
  if (ret) {
    std::cerr << "openRecord failed with return value " << ret << std::endl;
    return ret;
  }
  ret = 
    meta.addInputDesc("This is the standard invented hot strips file");
  if (ret) {
    std::cerr << "Bad return from addInputDesc: " << ret << std::endl;
    return ret;
  }
  ret = meta.addNotes("Fake record, added from test_meta");
  if (ret) {
    std::cerr << "Bad return from addNotes: " << ret << std::endl;
    return ret;
  }
  ret = meta.addValidInterval(Timestamp(2000, 8, 2), Timestamp());
  if (ret) {
    std::cerr << "Bad return from addValidInterval: " << ret << std::endl;
    return ret;
  }
  unsigned int newSerial;
  ret = meta.insertRecord(&newSerial);
  if (ret) {
    std::cerr << "Bad return from insertRecord: " << ret << std::endl;
  }
  else {
    std::cout << "Successfully inserted new record, serial number " 
              << newSerial << std::endl;
  }
  return(ret);

}


void gotIt(unsigned int ser, calibUtil::Metadata::eDataFmt dataFmt,
           const std::string& fmtVersion, const std::string& dataIdent) {

  std::cout << "Success reading info for record #" << ser << std::endl;
  
  std::cout << "Data format = " << dataFmt << std::endl;
  std::cout << "Format version = " << fmtVersion << std::endl;
  std::cout << "Data ident = " << dataIdent << std::endl;
}

calibUtil::Metadata::eRet lookup(calibUtil::Metadata::eCalibType ctype,
                                 const calibUtil::Timestamp& ts,
                                 unsigned int levelMask,
                                 calibUtil::Metadata::eInstrument inst) {
  using calibUtil::Metadata;
  unsigned int ser;

  
  std::cout << std::endl;
  std::cout << "lookup called with input " << std::endl;
  std::cout << "   calibType = " << ctype <<std::endl;
  std::cout << "   timestamp = " << ts.timeString() << std::endl;
  std::cout << "   levelMask = " << levelMask << std::endl;
  std::cout << "   instrument = " << inst << std::endl;

  Metadata       meta;
  Metadata::eRet ret = meta.findBest(&ser, ctype, ts, levelMask, inst);

  if (ret != Metadata::RETOk) {
    std::cout << "findBest failed with status" << ret << std::endl;
  }
  else if (!ser) {
    std::cout << "Query succeeded; no rows found." << std::endl;
  }
  else {
    std::string  fmtVersion;
    Metadata::eDataFmt dataFmt;
    std::string dataIdent;

    ret = meta.getReadInfo(ser, dataFmt, fmtVersion, dataIdent);
    
    if (ret == Metadata::RETOk) { 
      gotIt(ser, dataFmt, fmtVersion, dataIdent);
    }

    else {
      std::cout << "Failed reading info for record #" << ser;
      std::cout << " with code " << ret << std::endl;
    }
  }

  return ret;
}


