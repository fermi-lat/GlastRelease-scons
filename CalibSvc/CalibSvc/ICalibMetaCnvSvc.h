//$Header$
#ifndef ICalibMetaCnvSvc_h
#define ICalibMetaCnvSvc_h 1

#include "GaudiKernel/IInterface.h"

// External constants
#include "GaudiKernel/ClassID.h"

// Type definition
/////#include "ConditionsDB/CondDBKey.h"

static const InterfaceID IID_ICalibMetaCnvSvc ("ICalibMetaCnvSvc", 1, 0);

// Forward declarations
class DataObject;
// class IConditionsDBGate;
class IRegistry;
class ITime;

namespace calibUtil {
  class Metadata;
}


/** @class ICalibMetaCnvSvc 

    Abstract interface of a conversion service for GLAST calibration
    persistency.

    Handles creation and updating condition data objects (i.e. DataObjects
    implementing IValidity).

    Adapted from Andrea Valassi's IConditionsDBCnvSvc interface

    @author Joanne Bogart
    @date November 2002
*/

class ICalibMetaCnvSvc : virtual public IInterface
{

 public:

  // Re-implemented from IInterface

  /// Interface ID of this class
  static const InterfaceID& interfaceID() { return IID_ICalibMetaCnvSvc; }

 public:
  
  // Create/update condition MetadataObject not necessarily registered in 
  // the TDS.
  /// Create a calib DataObject -- corresponds to info from metadata -- by 
  /// --> whatever information is needed by Metadata::findBest,
  /// namely calibType, event time, instrument, flavor, maybe level
  //  [was folder name, tag and time for IConditionsDBCnvSvc]
  /// This method does not register DataObject in the transient data store,
  /// The string storage type is discovered at runtime in the MySQL metadata
  /// database.
 virtual StatusCode createCalibMetadata(DataObject*& refpObject,
                                        const std::string& calibType,
                                        const std::string& flavor,
                                        const ITime&       time,
                                        const std::string& instrumentName,
                                        const CLID&        classID,
                                        IRegistry*         entry=0) = 0;

  /// Update a condition DataObject by type, flavor, time (& instrument?).
  /// This method does not register DataObject in the transient data store,
  /// but may register TDS addresses for its children if needed (e.g. Catalog).
  /// The string storage type is discovered at runtime in the MySQL metadata
  /// database.
  virtual StatusCode updateCalibData(DataObject*         pObject,
                                     const std::string& calibType,
                                     const std::string& flavor,
                                     const ITime&       time,
                                     const std::string& instrumentName,
                                     const CLID&        classID,
                                     IRegistry*         entry=0) = 0;

  
  // The stuff below is maybe replaced by an abstract version of
  //   Metadata::getReadInfo
  /// Decode the string storage type from the folder description string
  /*
  virtual StatusCode decodeDescription(const std::string&  description,
                                       unsigned char&      type) = 0;
  */
  /// Encode the string storage type into the folder description string
  /*
  virtual StatusCode encodeDescription(const unsigned char& type,
                                       std::string&         description ) = 0;
  */


  // Get handle for metadata access from calibUtil.
  virtual calibUtil::Metadata* getMeta() = 0;

};

#endif




