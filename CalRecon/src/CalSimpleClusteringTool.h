/**
 * @class CalSimpleClusteringTool
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
#include "CalClusteringTool.h"


class CalSimpleClusteringTool : public CalClusteringTool
 {
  public :
  
    /// Standard Gaudi Tool interface constructor
    CalSimpleClusteringTool(
      const std::string& type,
      const std::string& name,
      const IInterface* parent ) ;
    virtual ~CalSimpleClusteringTool() ;
    
  protected:

    /// This finds the next highest energy cluster in a vector of CalXtalRecData pointers
    xTalDataVec            nextXtalsSet(xTalDataVec& xTalVec);

  private:
  
    xTalDataVec            getXtalsInLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInSameLayer(xTalDataVec& xTalVec, xTalDataVec& NNvec, Event::CalXtalRecData* xTal);
    Event::CalXtalRecData* getNearestXtalInDiffLayer(xTalDataVec& xTalVec, Event::CalXtalRecData* xTal, int layer);

 } ;

