// $Header$
#ifndef RootBaseCnv_h
#define RootBaseCnv_h

/** @class RootBaseCnv 

  Base class for calibration converters from XML files to TCDS.
  All such converters need to do certain things, which are
  handled here.

  @author J. Bogart
*/
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/Time.h"
#include <string>

#include "CalibSvc/ICalibRootSvc.h"
#include "TObject.h"

class ISvcLocator;
class GenericAddress;
class ICalibRootSvc;
class ICalibMetaCnvSvc;
class IInstrumentName;
//class ITime;

class TFile;
class TTree;
class TDirectory;
class TObject;

namespace CalibData {
  class CalibTime;
  class CalibBase;
  class DacCol;
}

class  RootBaseCnv : public Converter {


public:

  virtual ~RootBaseCnv();

  // Standard public converter stuff
  virtual StatusCode initialize();
  virtual StatusCode finalize();

  /**
   Create the transient representation of an object, given an opaque
   address.  This and the following update method comprise the core 
   functionality of calibration converters.
  */
  virtual StatusCode createObj(IOpaqueAddress* addr,
                               DataObject*& refpObject);

  ICalibRootSvc* getCalibRootSvc() {
    return m_rootSvc;
  }

  static const unsigned char storageType() {return CALIBROOT_StorageType;}

  /**
     Constructor for this converter
     @param svc a ISvcLocator interface to find services
     @param clid the type of object the converter is able to convert
   */
  RootBaseCnv(ISvcLocator* svc, const CLID& clid);

  // End standard public converter stuff

  /**
     Create ROOT file corresponding to TDS object input.  
     Default implementation is to return an error.  Must be separately
     implemented for each calibration type.
     

     @param fname     Filename for output file
     @param pTDSObj   Pointer to tds object to be converted
   */
  virtual StatusCode createRoot(const std::string& fname, 
                                CalibData::CalibBase* pTDSObj);

  /**
     Read in object (by default the first) from specified branch. 
  */
  virtual StatusCode readRootObj(const std::string& treename,
                                 const std::string& branch, TObject*& pCalib,
                                 unsigned index=0);

  virtual StatusCode readRootObj(TTree*  tree,
                                 const std::string& branch, TObject*& pCalib,
                                 unsigned index=0);

  /// Retrieve the class type of the data store the converter uses.
  virtual long repSvcType() const {return Converter::i_repSvcType();}

protected:
  /** This creates the transient representation of an object from the
   *  corresponding ROOT object it, then fills it and process it.
   *  This implementation actually only calls the i_* methods of the
   *  "right" converter to do the job; so the very first thing it
   *  does is get a pointer to the appropriate derived converter.
   *  Converters typically don't need to override this method
   *  but only to  override/implement some of the i_* methods.

   *  @param pRootObj pointer to the ROOT object
   *  @param refpObject the object to be built
   *  @param address the opaque address for this object
   *  @return status depending on the completion of the call
   */
  virtual StatusCode internalCreateObj (const std::string& fname,
                                        // was TObject*     pRootObj,
                                        DataObject*& refpObject,
                                        IOpaqueAddress* address);
  
  /** This creates the transient representation of an object from the
   *  corresponding ROOT object. This actually does the "new" operation
   *  and deals with the attributes of the node. This base class implementation
   *  does nothing; it should not normally be called because it doesn't
   *  correspond to any TCDS class.  Instead, 
   *  i_createObj of some derived class will be called.
   *  @param fname The ROOT file to be read in
   *  to be used to builds the object
   *  @param refpObject the object to be built
   *  @return status depending on the completion of the call
   */
  virtual StatusCode i_createObj (const std::string& fname,
                                  DataObject*& refpObject);

  /// In case there is additional work to do on the created object
  virtual StatusCode i_processObj(DataObject* pObject,
                                  IOpaqueAddress* address);

  /**
     Given a pointer to a TDS object which can be cast to "our" type, fill
     in corresponding information in the corresponding root class

     @param pTDSObj   Pointer to tds object to be converted
     @param pRootObj  Pointer to destination root object

     ...maybe don't need pRootObj argument; keep this as the (protected)
     data member m_rootObj. Or else this routine could set the protected
     member to this passed-in value
  */
  virtual StatusCode fillRoot(CalibData::CalibBase* pTDSObj, 
                              TObject* pRootObj);


  /**
     Utility used by derived converters to start writing a ROOT file
     (open TFile, make a TTree, give it a branch)
     @param fname       Name for new file
     @param className   Name of class for object specified in next parameter;
                        used to name branch as well)
     @param pCalib      pointer to object used to create the branch

  */
  virtual StatusCode openWrite(const std::string& fname, 
                               const std::string& className,
                               TObject*& pCalib);

  /**
     Finish up writing file opened with openWrite:
            fill the tree
            write the file
            close the file
            Delete TFile (causes associated Tree to be deleted)
  */
  StatusCode closeWrite();
  

  /**
     Utility for "leaf" converters to call
     @param    Root file to open for read
//     @param    Name of branch to be read in
//     @param    ref. to pCalib pointer which will be set to address of
               read-in object
  */
  //  StatusCode openRead(const std::string& fname, const std::string& branch,
  //                      TObject*& pCalib);
  StatusCode openRead(const std::string& fname);

  /** Clean up when we've finished reading in */
  StatusCode closeRead();

  /// Another convenience for derived classes: sets information belonging
  /// to the calibration base class, namely validity interval and serial
  /// number.
  void setBaseInfo(CalibData::CalibBase* pObj);

  // Might want to verify that instrument, calType are correct,
  // for example.  If so, might as well provide the service in
  // the base converter.
  //////  virtual StatusCode readHeader(const DOM_Element&);


  ICalibRootSvc* m_rootSvc;
  ICalibMetaCnvSvc* m_metaSvc;
  IInstrumentName* m_instrSvc;

  int m_serNo;
  Gaudi::Time*  m_vstart;
  Gaudi::Time*  m_vend;

  // Note by keeping this stuff here we're disallowing possibility of 
  // interleaved writes of different calibrations
  TFile*  m_outFile;
  TTree*  m_ttree;

  TFile*  m_inFile;

  TDirectory* m_saveDir;

  // Don't think we want this after all
  // store pointer to intermediate (between TDS and permanent file) form
  // of calibration data as root obj.  It will normally be the responsibility
  // of the "leaf" calibration converter to create and delete it.
  //    calibRootData::Base* m_rootObj;

  // Is it reasonable to always call our TTree "Calib" ?
  // Branches will be named after calibration type
private:
  // Return true if something was there to be cleaned
  bool doClean();
  
};

#endif
