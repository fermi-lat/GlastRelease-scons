//	============================================================
//	DbEventCnvSvc.h
//	------------------------------------------------------------
//
//	Package   : Db Converters
//
//	Author    : Markus Frank
//	============================================================
#ifndef DBCNV_DBEVENTCNVSVC_H
#define DBCNV_DBEVENTCNVSVC_H 1
// $Header$

// Include files
#include "DbCnvSvc.h"

// Forward declarations
template <class TYPE> class SvcFactory;

/** DbEventCnvSvc class implementation definition.

    <P> System:  The LHCb Offline System
    <P> Package: PersistencySvc/Db

     Description:
     The DbEventCnvSvc is the specialized version of the 
     DbCnvSvc used by the event data store to load event data.

     Base class:
     DbCnvSvc

     Dependencies:
    <UL>
    <LI> DbCnvSvc definition:  DbCnv/DbSvc/DbCnvSvc.h
    </UL>

    History:
    <PRE>
    +---------+----------------------------------------------+---------+
    |    Date |                 Comment                      | Who     |
    +---------+----------------------------------------------+---------+
    | 3/11/98 | Initial version                              | M.Frank |
    +---------+----------------------------------------------+---------+
    </PRE>
   @author Markus Frank
   @version 1.0
*/
class DbEventCnvSvc : public DbCnvSvc	{
  /// Allow the factory class access to the constructor
  friend SvcFactory<DbEventCnvSvc>;
public:
  /**@name DbEventCnvSvc overrides     */
  //@{
  /// Initialize the service.
  virtual StatusCode initialize();

  /// Update state of the service
  virtual StatusCode updateServiceState(IOpaqueAddress* pAddress);
  //@}

protected:
  /**@name: Object implementation     */
  //@{
  /// Standard Constructor
  DbEventCnvSvc(const std::string& name, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~DbEventCnvSvc();
  //@}
};

#endif // DBCNV_DBEVENTCNVSVC_H
