//##############################################################
//
// Minimal job options file for reading root events + calib services
//

// List of Dlls
ApplicationMgr.DLLs   = { "GlastSvc", "RootIo" };
ApplicationMgr.DLLs   += {"CalibSvc"};

// List of required services
 ApplicationMgr.ExtSvc   = {"GlastDetSvc", "GlastEventSelector/EventSelector" , "EventCnvSvc" };

// Use the new RootIoSvc, implementing IRunable
ApplicationMgr.ExtSvc += { "RootIoSvc" };
ApplicationMgr.Runable= "RootIoSvc"; 

// There are several other possibilities for test algroithm, including
//  UseBadStrips, UseGains,...
ApplicationMgr.TopAlg = {"mcRootReaderAlg", "UsePeds" };
// digiRootReaderAlg.digiRootFile = "$CALIBSVCROOT/src/test/digi_1409.root";
mcRootReaderAlg.mcRootFileList = {"$CALIBSVCROOT/src/test/mc_1409.root"};


// Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
MessageSvc.OutputLevel      = 2;

//=========================================================================
// Persistency service setup:
EventPersistencySvc.CnvServices = {"EventCnvSvc"};

// Next few lines from AcdDigi test jobOptions
EventSelector.Input = "NONE";

GlastDetSvc.xmlfile="$(XMLGEODBSROOT)/xml/flight/flightSegVols.xml";
GlastDetSvc.topVolume="LAT";

//ApplicationMgr.EvtMax = 10000;
ApplicationMgr.EvtMax = 20;

// Model for immediately following lines is LHCb example.  See
// http://lhcbsoft.web.cern.ch/LHCbSoft/Ex/DetConExamples/v2r0/options/
//   in particular, files contDB.opts and common.opts

// Start up a CalibDataSvc 
ApplicationMgr.ExtSvc += {"CalibDataSvc"};

ApplicationMgr.ExtSvc += {"CalibMySQLCnvSvc", "CalibXmlCnvSvc"};
ApplicationMgr.ExtSvc += { "CalibRootCnvSvc" };

DetectorPersistencySvc.CnvServices += {"CalibMySQLCnvSvc"};
DetectorPersistencySvc.CnvServices += {"CalibXmlCnvSvc"};
DetectorPersistencySvc.CnvServices += {"CalibRootCnvSvc"};

CalibDataSvc.CalibFlavorList = {"ideal"};
// not used yet
CalibDataSvc.CalibNameList = {"TKR_HotChan/chocolate",
			   "Test_1/mocha"};

// CalibDataSvc properties below are explicitly set to the
// default values as established in CalibDataSvc.cxx.  
// They're listed below in order of likelihood of need to override, 
// from most to least likely.
// Storage type of 14 corresponds to MYSQL
CalibDataSvc.CalibInstrumentName = "LAT";
CalibDataSvc.CalibStorageType = 14;
CalibDataSvc.CalibRootName = "Calib";

// Value of "*" means 'use default MySQL host', so currently (May, 2003)
// equivalent to value of "centaurusa.slac.stanford.edu".
// For local MySQL server, use value "localhost"
CalibMySQLCnvSvc.Host = "*";

// By default, use production database "calib", but user may override.
//  calib_user is a development database for use of, e.g., calibration
//  algorithm developers
// CalibMySQLCnvSvc.DbName = "calib_user";

CalibDataSvc.CalibTimeSource = "mc";
/////CalibDataSvc.CalibTimeSource = "clock";
/////CalibDataSvc.startTime = "2003-2-25 00:20";
/////CalibDataSvc.delayTime = 900;

// In order to circumvent normal calibration-fetching mechanism, which
// relies on event timestamp, you *must* include the first two lines
// below.  The third and fourth lines have defaults of 
//   "2003-09-01" (calibration must have been entered sept. 1 2003 or later)
//          and 
//   forever
// CalibDataSvc.UseEventTime = false;
// CalibMySQLCnvSvc.UseEventTime = false;

// CalibMySQLCnvSvc.EnterTimeStart = "2003-01-01";
// CalibMySQLCnvSvc.EnterTimeEnd = "2004-03-01";


//==============================================================
//
// End of job options file
//
//##############################################################

