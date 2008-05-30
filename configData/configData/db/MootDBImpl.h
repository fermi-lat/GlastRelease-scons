/** @class MootDBImpl
    @author M. Kocian
    

    Implementation to retrieve sets of files from MOOT

    $Id$
*/
#ifndef MOOTDBIMPL_H
#define MOOTDBIMPL_H

#include "configData/db/LatcDB.h"
#include <string>
#include <map>
#include "mootCore/MootQuery.h"
#include "mootCore/FileDescrip.h"
#include "enums/Lsf.h"

class FswEfcSampler;

using namespace MOOT;

class MootDBImpl:public LatcDB{
public:
  /// Default constructor
  MootDBImpl();

  /// D'tor constructor
  virtual ~MootDBImpl();

  /// function to retrieve filename of type type given key key.
  const std::string getFilename(const char* type, unsigned key);

  /// function to retrieve the list of filenames given a MOOT key
  const std::vector<std::string> getFilenameList(unsigned key);

  /// function to get the RunSummaryInfo and the MOOT master key
  /// given the run start time and the source ID
  AcqSummaryInfo* getAcqSummary(unsigned startedAt,
				unsigned scid, 
				unsigned& configKey) const;
  
  /// function to get the prescale factor
  /// given the MOOT master key
  unsigned prescaleFactor(unsigned configKey,			  			  
			  MOOT::LpaMode lpaMode,
			  unsigned handlerId,			  
			  unsigned filterStatusWord,
			  enums::Lsf::RsdState rsdState,
			  enums::Lsf::LeakedPrescaler leakPrescaler) const;

  bool readEfcFromFile( MOOT::LpaMode lpaMode, unsigned handlerId, const char* fileName );
  
protected:

  void clearCache() const;
  
  FswEfcSampler* getFswCdb(unsigned configKey, MOOT::LpaMode lpaMode, unsigned handlerId) const;
  
private:
  
  mutable std::map<unsigned,FswEfcSampler*>       m_fswEfcCache;
  mutable AcqSummaryInfo*                         m_acqSummaryCache;
  mutable unsigned m_configKeyCache;
  mutable unsigned m_startedAtCache;
  mutable unsigned m_scidCache;

  MootQuery* m_mq;
};

#endif
