/**
 * @class CalSimpleClustering
 *
 * @brief Implements a Gaudi Tool for performing very simple clustering in the Cal 
 *
 * @author Tracy Usher
 *
 * File and Version Information:
 *      $Header$
 */

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "CalClustering.h"


class CalSimpleClustering : public CalClustering
 {
  public :
  
    /// Standard Gaudi Tool interface constructor
    CalSimpleClustering(
      const std::string& type,
      const std::string& name,
      const IInterface* parent ) ;
    virtual ~CalSimpleClustering() ;
    
  protected:

    //! Distinguish sets of related xtals
    virtual void makeSets( const XtalDataVec & xtals, XtalDataVecVec & clusters ) ;

  private:
  
    XtalDataVec            getXtalsInLayer(XtalDataVec& xTalVec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInSameLayer(XtalDataVec& xTalVec, XtalDataVec& NNvec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInDiffLayer(XtalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer);

 } ;

