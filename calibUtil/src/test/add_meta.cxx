// $Header$
/** @file test_meta.cxx
    Sample program to add record(s) to the MySQL database. For now user
    only gets to specify
           instrument type
           location of file
           validity start and end times

    Ultimately need something more flexible and easier to use.

*/

#include <string>
#include <iostream>
#include "calibUtil/Metadata.h"


void gotIt(unsigned int ser, calibUtil::Metadata::eDataFmt dataFmt,
           const std::string& fmtVersion, const std::string& dataIdent);

// bundle up all the fields which must be specified by caller to
// make a new record
typedef struct s_MetaInput {
  calibUtil::Metadata::eInstrument instr;
  calibUtil::Metadata::eCalibType  calibType;
  calibUtil::Metadata::eDataFmt    fmt;
  calibUtil::Metadata::eCompletion cmpl;
  calibUtil::Metadata::eLevel      procLevel;
  std::string           path;
  std::string           inputDesc;
  std::string           notes;
  facilities::Timestamp vstart;
  facilities::Timestamp vend;
  std::string           flavor;
}   MetaInput;

bool getUserInput(MetaInput* inp);

int main(int, char**) {
  using calibUtil::Metadata;
  using facilities::Timestamp;

  Metadata  meta;

  MetaInput userInput;

  while (getUserInput(&userInput)) {
    Metadata::eRet ret;
    // Try to insert a record
    ret = meta.openRecord(userInput.instr, userInput.calibType,
                          userInput.fmt, "1.0",
                          userInput.path,
                          userInput.cmpl,
                          userInput.procLevel);
    if (ret) {
      std::cerr << "openRecord failed with return value " << (int) ret 
                << std::endl;
      return ret;
    }  
    ret = 
      meta.addInputDesc(userInput.inputDesc);
    if (ret) {
      std::cerr << "Bad return from addInputDesc: " << (int) ret << std::endl;
      return ret;
    }
    ret = meta.addNotes(userInput.notes);
    if (ret) {
      std::cerr << "Bad return from addNotes: " << (int) ret << std::endl;
      return ret;
    }
    ret = meta.addValidInterval(userInput.vstart, userInput.vend);
    if (ret) {
      std::cerr << "Bad return from addValidInterval: " << (int) ret 
                << std::endl;
      return ret;
    }
    
    ret = meta.addFlavor(userInput.flavor);
    
    unsigned int newSerial;
    ret = meta.insertRecord(&newSerial);
    if (ret) {
      std::cerr << "Bad return from insertRecord: " << (int) ret << std::endl;
    }
    else {
      std::cout << "Successfully inserted new record, serial number " 
                << newSerial << std::endl;
    }

  }
}

// Get info needed to make a new record in the MySQL database interactively.
// Some day might want to allow input from a file as well
bool getUserInput(MetaInput* inp) {
  using calibUtil::Metadata;

  // Put in some defaults which user may override
  inp->instr=Metadata::INSTLat;
  inp->calibType=Metadata::CTYPE_TestGen;
  inp->fmt = Metadata::FMTXml;
  inp->cmpl = Metadata::CMPLOk;
  inp->procLevel = Metadata::LEVELProd;
  inp->flavor = "vanilla";
  inp->inputDesc = "No real input description";
  inp->notes = "No real notes";

  std::cout << "Instrument: LAT (" << (int) Metadata::INSTLat << ") BTEM (";
  std::cout << (int) Metadata::INSTBtem <<  ") BFEM (" 
            << (int) Metadata::INSTBfem;
  std::cout << ") EM (" << (int) Metadata::INSTEm << ")";
  std::cout << "CU (" << (int) Metadata::INSTCu << ")"  << std::endl;
  int inst;
  std::cout << "Your type:  (or -1 to quit):  ";
  std::cin >> inst;
  if (inst < 0) return false;

  inp->instr = (Metadata::eInstrument) inst;

  std::cout << "Calib type is TestGen (" 
            << (int) Metadata::CTYPE_TestGen << ")";
  std::cout << std::endl;

  std::cout << 
    "Path to XML file, e.g. $(CALIBUTILROOT)/xml/test/gen1_1.xml (required): ";
  std::cin >> inp->path;
  std::string ascTime;
  std::cout << 
    "Validity interval start timestamp ending with *, e.g. 2003-01-12 08:00*";
  std::cout << std::endl;

  // std::cin always seems to have a \n at the front, even after
  // std::cin.sync(), causing getline to come back immediately with no 
  // data unless we use a non-standard character to signal end of input.
  // std::cout << "Validity interval start timestamp, e.g. 2003-01-12 08:00" 
  //     << std::endl;
  //  std::cin.sync();
  // std::getline(std::cin, ascTime);

  std::getline(std::cin, ascTime, '*');
  inp->vstart = facilities::Timestamp(ascTime);
  std::cout << 
    "Validity interval end timestamp ending with *, e.g. 2003-01-12 17:00*  :";
  std::cout << std::endl;
  std::getline(std::cin, ascTime, '*');
  inp->vend = facilities::Timestamp(ascTime);
  return true;
}


void gotIt(unsigned int ser, calibUtil::Metadata::eDataFmt dataFmt,
           const std::string& fmtVersion, const std::string& dataIdent) {

  std::cout << "Success reading info for record #" << ser << std::endl;
  
  std::cout << "Data format = " << (int) dataFmt << std::endl;
  std::cout << "Format version = " << fmtVersion << std::endl;
  std::cout << "Data ident = " << dataIdent << std::endl;
}


