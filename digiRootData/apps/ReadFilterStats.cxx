#include <string>
#include <map>
#include <list>

#include "enums/Lsf.h"
#include "digiRootData/DatagramInfo.h"
#include "digiRootData/RunInfo.h"

#include <TFile.h>
#include <TTree.h>
#include <TObjString.h>

#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
#endif

struct readData;
struct fillData;

void usage(std::ostream& os){
  os << std::endl;
  os << "Usage:  ReadFilterStats.exe -i <inputFile> -o <outputFile>" << std::endl;
  return;
};

/**
 * Go from string to enum for the Platfom ( lat | testbed | host )
 *
 **/
enums::Lsf::Platform        run_platform_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Platform> theMap;
  if ( theMap.size() == 0 ) {
    theMap["lat"]     = enums::Lsf::Lat;
    theMap["testbed"] = enums::Lsf::Testbed;
    theMap["host"]    = enums::Lsf::Host;    
  }
  std::map<std::string,enums::Lsf::Platform>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::NoPlatform;
};


/**
 * Go from string to enum for the Data Origin ( orbit | montecarlo | ground  )
 *
 **/  
enums::Lsf::DataOrigin      run_origin_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::DataOrigin> theMap;
  if ( theMap.size() == 0 ) {
    theMap["orbit"]      = enums::Lsf::Orbit;
    theMap["montecarlo"] = enums::Lsf::MonteCarlo;
    theMap["ground"]     = enums::Lsf::Ground;
  }
  std::map<std::string,enums::Lsf::DataOrigin>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::NoOrigin;  
};

/**
 * Go from string to enum for the Datagram Open action ( start | resume | continue  )
 *
 **/  
enums::Lsf::Open::Action    open_action_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Open::Action> theMap;
  if ( theMap.size() == 0 ) {
    theMap["start"] = enums::Lsf::Open::Start;
    theMap["resume"] = enums::Lsf::Open::Resume;
    theMap["continue"] = enums::Lsf::Open::Continue;
  }
  std::map<std::string,enums::Lsf::Open::Action>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::Open::Unspecified;  
};

/**
 * Go from string to enum for the Datagram Open reason 
 *     ( operator | modechange | timedout | countedout | full )
 *
 **/  
enums::Lsf::Open::Reason    open_reason_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Open::Reason> theMap;
  if ( theMap.size() == 0 ) {
    theMap["operator"]   = enums::Lsf::Open::Operator;
    theMap["modechange"] = enums::Lsf::Open::ModeChange;
    theMap["timedout"]   = enums::Lsf::Open::TimedOut;
    theMap["countedout"] = enums::Lsf::Open::CountedOut;
    theMap["full"]       = enums::Lsf::Open::Full;
  }
  std::map<std::string,enums::Lsf::Open::Reason>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::Open::Unknown;  
};

/**
 * Go from string to enum for the Datagram Close action ( stop | abort | pause | continue  )
 *
 **/  
enums::Lsf::Close::Action   close_action_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Close::Action> theMap;
  if ( theMap.size() == 0 ) {
    theMap["stop"]     = enums::Lsf::Close::Stop;
    theMap["abort"]    = enums::Lsf::Close::Abort;
    theMap["pause"]    = enums::Lsf::Close::Pause;
    theMap["continue"] = enums::Lsf::Close::Continue;
  }
  std::map<std::string,enums::Lsf::Close::Action>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::Close::Unspecified;  
};

/**
 * Go from string to enum for the Datagram Close reason 
 *     ( operator | modechange | timedout | countedout | full )
 *
 **/  
enums::Lsf::Close::Reason   close_reason_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Close::Reason> theMap;
  if ( theMap.size() == 0 ) {
    theMap["operator"]   = enums::Lsf::Close::Operator;
    theMap["modechange"] = enums::Lsf::Close::ModeChange;
    theMap["timedout"]   = enums::Lsf::Close::TimedOut;
    theMap["countedout"] = enums::Lsf::Close::CountedOut;
    theMap["full"]       = enums::Lsf::Close::Full;
  }
  std::map<std::string,enums::Lsf::Close::Reason>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::Close::Unknown;  
};  

/**
 * Go from string to enum for the proccesing crate
 *     ( epu[0-2] | siu[0-1] | aux | mixed )
 *
 **/  
enums::Lsf::Crate           open_crate_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Crate> theMap;
  if ( theMap.size() == 0 ) {
    theMap["epu0"]  = enums::Lsf::Epu0;
    theMap["epu1"]  = enums::Lsf::Epu1; 
    theMap["epu2"]  = enums::Lsf::Epu2; 
    theMap["siu0"]  = enums::Lsf::Siu0; 
    theMap["siu1"]  = enums::Lsf::Siu1; 
    theMap["aux"]   = enums::Lsf::Aux;
    theMap["mixed"] = enums::Lsf::Mixed; 
  }
  std::map<std::string,enums::Lsf::Crate>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::NoCrate;  
};

/**
 * Go from string to enum for the datagram mode
 *     ( normal | too | grb[0-2] | solar | calibration | diagnostic )
 *
 **/  
enums::Lsf::Mode            open_mode_enum(const std::string& str){
  static std::map<std::string,enums::Lsf::Mode> theMap;
  if ( theMap.size() == 0 ) {
    theMap["normal"]      = enums::Lsf::Normal;
    theMap["too"]         = enums::Lsf::TOO;
    theMap["grb0"]        = enums::Lsf::GRB0;
    theMap["grb1"]        = enums::Lsf::GRB1;
    theMap["grb2"]        = enums::Lsf::GRB2;
    theMap["solar"]       = enums::Lsf::Solar;
    theMap["calibration"] = enums::Lsf::Calibration;
    theMap["diagnostic"]  = enums::Lsf::Diagnostic;
  }
  std::map<std::string,enums::Lsf::Mode>::const_iterator itr = theMap.find(str);
  return itr != theMap.end() ? itr->second : enums::Lsf::NoMode;  
};


/**
 * fillData
 *     Data structure used to fill ntpule
 *
 **/  

struct fillData {
public:
  /// Simple c'tor, need name for TTree
  fillData(const char* treeName)
    :m_tree(0),
     m_runInfo(new RunInfo),
     m_datagramInfo(new DatagramInfo),
     m_handlerName(new TObjString){
    TString tN(treeName);
    tN.ReplaceAll("-","_");
    m_tree = new TTree(tN,"Filter statistics");
  }
  /// Simple d'tor, clean up
  ~fillData(){
    delete m_runInfo;
    delete m_datagramInfo;
    delete m_handlerName;
    m_tree = 0;
  }
  /// return the TTree this data is filling
  inline TTree* getTree() { return m_tree; }

  /// Declare branches and return the TTree this data is filling
  TTree* fillTree( ) {    
    m_tree->Bronch("RunInfo",RunInfo::Class()->GetName(),&m_runInfo);
    m_tree->Bronch("DatagramInfo",DatagramInfo::Class()->GetName(),&m_datagramInfo);
    m_tree->Bronch("HandlerName",TObjString::Class()->GetName(),&m_handlerName);

    m_tree->Branch("HandlerVersion",&m_handlerVersion,"HandlerVersion/i");
    m_tree->Branch("handlerId",&m_handlerId,"handlerId/i");
    m_tree->Branch("iRecord",&m_iRecord,"iRecord/i");
    m_tree->Branch("nRecords",&m_nRecords,"nRecords/i");
    m_tree->Branch("nMergedRecords",&m_nMergedRecords,"nMergedRecords/i");
    
    m_tree->Branch("beginSequence",&m_beginSequence,"beginSequence/i");
    m_tree->Branch("endSequence",&m_endSequence,"endSequence/i");
    
    m_tree->Branch("beginTime",&m_beginTime,"beginTime/l");
    m_tree->Branch("endTime",&m_endTime,"endTime/l");
    
    m_tree->Branch("Passed",&m_nPassed,"Passed[32]/i");
    m_tree->Branch("Vetoed",&m_nVetoed,"Vetoed[32]/i");
    return m_tree;
  }

  /// Fill this object with data read from csv file, return true for success
  bool convert(const readData& rData, bool first);

  /// The tree this stuff goes on
  TTree*                      m_tree;

  // Data going into ROOT file  

  RunInfo*                    m_runInfo;
  DatagramInfo*               m_datagramInfo;
  TObjString*                 m_handlerName;
  UInt_t                      m_handlerVersion;
  UInt_t                      m_handlerId;  

  UInt_t                      m_iRecord;
  UInt_t                      m_nRecords;
  UInt_t                      m_nMergedRecords;

  ULong64_t                   m_beginTime;
  UInt_t                      m_beginSequence;  
  ULong64_t                   m_endTime;
  UInt_t                      m_endSequence;

  UInt_t                      m_nPassed[32];
  UInt_t                      m_nVetoed[32];
};


struct readData {
public:

  /// Simple c'tor
  readData(){;}
  /// Simple c'tor
  ~readData(){;}

  /// Read a single line for a csv file, return true for success
  bool readline(std::istream& is, bool first) {

    if ( is.eof() ) return true;
    if ( !is.good() ) {
      std::cerr << "Failed before read" << std::endl;
      return false;
    }

    is >> m_run_startedAt >> m_run_softwareKey >> m_run_platform >> m_run_origin;
    is >> m_open_crate >> m_open_mode 
       >> m_handlerId >> m_handlerVersion 
       >> m_open_masterKey >> m_open_cfgKey >> m_open_cfgId >> m_open_datagrams 
       >> m_open_action >> m_open_reason >> m_close_action >> m_close_reason;
    is >> m_handlerName;
    is >> m_iRecord >> m_nRecords;    
    is >> m_cur_beginSequence >> m_cur_beginTime 
       >> m_cur_endSequence >> m_cur_endTime 
       >> m_cur_dSequence >> m_cur_dTime;
    is >> m_cum_beginSequence >> m_cum_beginTime 
       >> m_cum_endSequence >> m_cum_endTime 
       >> m_cum_dSequence >> m_cum_dTime;

    UInt_t i(0);
    for ( i = 0 ; i < 32; i++ ) {
      is >> m_curPassed[i];
    }
    for ( i = 0 ; i < 32; i++ ) {
      is >> m_curVetoed[i];
    }
    for ( i = 0 ; i < 32; i++ ) {
      m_prevPassed[i] = first ? 0 : m_cumPassed[i];
      is >> m_cumPassed[i];    
      m_deltaPassed[i] = m_cumPassed[i] - m_prevPassed[i];
    }
    for ( i = 0 ; i < 32; i++ ) {
      m_prevVetoed[i] = first ? 0 : m_cumVetoed[i];
      is >> m_cumVetoed[i];
      m_deltaVetoed[i] = m_cumVetoed[i] - m_prevVetoed[i];
    }
    if ( is.eof() ) return true;
    if ( ! is.good() )  {
      std::cerr << "Bad line " << std::endl;
      return false;
    }
    return true;
  }

  /// read an entire csv file, filling the TTree in the fData struct
  bool readFile( const char* fileName, fillData& fData, bool firstFile) {
    std::ifstream is(fileName);
    if ( ! is.good() ) {
      std::cerr << "File Open of " << fileName << " failed." << std::endl;
      return false;
    }
    TTree* tout = firstFile ? fData.fillTree() : fData.getTree();
    
    // burn the first line
    char lineBuf[8192];  
    is.getline(lineBuf,8192);
    int i(0);
    bool first = firstFile;
    while ( ! is.eof() ) {
      i++;
      if ( ! readline(is,first) ) {
	std::cerr << "readline failed at line " << i << " in file " << fileName << std::endl;
	return false;
      }
      if ( ! fData.convert(*this,first) ) {
	std::cerr << "convert failed at line " << i << " in file " << fileName << std::endl;
	return false;
      }
      first = false;
      tout->Fill();
    }
    is.close();
    return true;
  }

  // Data read from csv file

  ULong64_t                   m_run_startedAt;
  UInt_t                      m_run_softwareKey;  
  std::string                 m_run_platform;
  std::string                 m_run_origin;
  std::string                 m_open_crate;
  std::string                 m_open_mode;
  UInt_t                      m_open_masterKey;
  UInt_t                      m_open_cfgKey;
  UInt_t                      m_open_cfgId;
  UInt_t                      m_open_datagrams;
  std::string                 m_open_action;
  std::string                 m_open_reason;
  std::string                 m_close_action;
  std::string                 m_close_reason;  
  UInt_t                      m_handlerVersion;
  UInt_t                      m_handlerId;  
  UInt_t                      m_iRecord;
  UInt_t                      m_nRecords;
  ULong64_t                   m_cur_beginTime;
  UInt_t                      m_cur_beginSequence;
  ULong64_t                   m_cur_endTime;
  UInt_t                      m_cur_endSequence;
  ULong64_t                   m_cur_dTime;
  UInt_t                      m_cur_dSequence;
  ULong64_t                   m_cum_beginTime;
  UInt_t                      m_cum_beginSequence;
  ULong64_t                   m_cum_endTime;
  UInt_t                      m_cum_endSequence;
  ULong64_t                   m_cum_dTime;
  UInt_t                      m_cum_dSequence;
  std::string                 m_handlerName;  
  UInt_t                      m_curPassed[32];
  UInt_t                      m_curVetoed[32];
  UInt_t                      m_cumPassed[32];
  UInt_t                      m_cumVetoed[32];

  // Caches

  UInt_t                      m_prevPassed[32];
  UInt_t                      m_prevVetoed[32];  
  UInt_t                      m_deltaPassed[32];
  UInt_t                      m_deltaVetoed[32];
};

/**
 * Use the data read from csv file to fill the structure in the output TTree
 *
 **/

bool fillData::convert(const readData& rData, bool first){  
  m_runInfo->initialize( run_platform_enum(rData.m_run_platform),
			 run_origin_enum(rData.m_run_origin),
			 0,
			 rData.m_run_startedAt );
  m_datagramInfo->initialize(  open_action_enum(rData.m_open_action),
			       open_reason_enum(rData.m_open_reason),
			       open_crate_enum(rData.m_open_crate),				  
			       open_mode_enum(rData.m_open_mode),
			       close_action_enum(rData.m_close_action),
			       close_reason_enum(rData.m_close_reason),
			       rData.m_open_datagrams,
			       0 );
  
  m_handlerVersion = rData.m_handlerVersion;
  m_handlerId = rData.m_handlerId;
  m_nMergedRecords = first ? m_iRecord+1 : rData.m_iRecord - m_iRecord;
  m_iRecord = rData.m_iRecord;
  m_nRecords = rData.m_nRecords;
  m_beginTime = first ? rData.m_cur_beginTime : m_endTime;
  m_beginSequence = first ? rData.m_cur_beginSequence : m_endSequence;
  m_endTime = rData.m_cur_endTime;
  m_endSequence = rData.m_cur_endSequence;  
  m_handlerName->SetString( rData.m_handlerName.c_str() );
  for ( UInt_t i(0); i < 32; i++ ) {
    m_nPassed[i] = rData.m_deltaPassed[i];
    m_nVetoed[i] = rData.m_deltaVetoed[i];
  }
  return true;
}

/**
 * Read the input file, which is just a list of csv files.
 *
 * return true for success
 **/
bool readInputFile(const char* inFileName, std::list<std::string>& csvFiles) {

  std::ifstream is(inFileName);
  if ( ! is.good() ) {
    std::cerr << "File Open of " << inFileName << " failed." << std::endl;
    return false;
  }

  char lineBuf[200];    
  while ( ! is.eof() ) {  
    is.getline(lineBuf,8192);
    csvFiles.push_back(lineBuf);
  }

  if ( is.eof() ) return true;
  if ( ! is.good() ) return false;
  is.close();
  return true;
}

/**
 *  main() 
 *
 *  return codes:
 *     0 -> success
 *     1 -> printed help message and exited
 *     2 -> bad arguments
 *     3 -> failed to read input file
 *     4 -> failed to open output file
 *     5 -> falied to read one of the csv data files
 **/
int main(int argn, char** argc) {

  // files names
  std::string inFileName;  // input file (with the list of csv files)
  std::string outFileName; // output file (ROOT file)
  std::list<std::string> csvFileNames; // list of csv files with FilterStatistics data

  // Option parsing
  int opt;
#ifdef WIN32
  while ( (opt = facilities::getopt(argn, argc, "hi:o:")) != EOF ) {
#else
  while ( (opt = getopt(argn, argc, "hi:o:")) != EOF ) {
#endif
    switch (opt) {
    case 'h':   // help      
      usage(std::cout);
      return 1;
    case 'i':   //  input
      inFileName = optarg;
      break;
    case 'o':   //  output
      outFileName  = optarg;
      break;
    case '?':
      usage(std::cout);
      return 2;
    default:
      std::cerr << opt << " not parsable..." << std::endl;
      usage(std::cerr);
      return 2;  
    }
  }

  if ( optind != argn  ) {
    usage(std::cerr);
    return 2;    
  }

  if ( outFileName.size() == 0 ||
       inFileName.size() == 0 ) {
    usage(std::cerr);
    return 2; 
  }


  // get the list of csv files from the inputFile
  if ( ! readInputFile(inFileName.c_str(),csvFileNames) ) {
    std::cerr << "Failed to read input file " << inFileName << std::endl;
    return 3;
  }

  // parse the csv files by type
  std::list<std::string> sourceList;
  sourceList.push_back("epu0-stats-h01");
  sourceList.push_back("epu0-stats-h04");
  sourceList.push_back("epu0-stats-h05");
  sourceList.push_back("epu1-stats-h01");
  sourceList.push_back("epu1-stats-h04");
  sourceList.push_back("epu1-stats-h05");
  
  std::map<std::string,std::list<std::string> > fileMap;
  for ( std::list<std::string>::const_iterator itr = csvFileNames.begin(); 
	itr != csvFileNames.end(); itr++ ) {
    for ( std::list<std::string>::const_iterator itrFind = sourceList.begin();
	  itrFind != sourceList.end(); itrFind++ ){
      if ( itr->find(*itrFind) != std::string::npos ) {
	fileMap[*itrFind].push_back(*itr);
	std::cout << *itrFind << ' ' << *itr << std::endl;
      }
    } 
  }
  

  // Open the output file
  TFile* fout = TFile::Open(outFileName.c_str(),"RECREATE");
  if ( fout == 0 ) {
    std::cerr << "Open output " << outFileName << " failed" << std::endl;
    return 4;
  }

  // Read the csv files and fill the output file
  readData rData; 
  for ( std::map<std::string,std::list<std::string> >::const_iterator itrF = fileMap.begin(); 
	itrF != fileMap.end(); itrF++ ) {
    // get the list of files of this type
    const std::list<std::string>& sList = itrF->second;
    fillData* nF = new fillData( itrF->first.c_str() );
    bool first(true);
    // loop on files of this type
    for ( std::list<std::string>::const_iterator itrFF = sList.begin(); itrFF != sList.end(); itrFF++ ) {
      if ( ! rData.readFile( itrFF->c_str(), *nF, first ) ) {
	return 5;
      }      
    }
    nF->getTree()->Write();
  }

  // all done, clean up and exit
  fout->Close();
  return 0;
}







