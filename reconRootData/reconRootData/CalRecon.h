#ifndef CALRECON_H
#define CALRECON_H

#include "TObject.h"
#include "TObjArray.h"

#include <vector>

#include "CalXtalRecData.h"
#include "CalCluster.h"

/** @class CalRecon
 * @brief Primary Root object containing CAL reconstruction data.
   portal through which one gains access to the CAL recon data. 
   Currently there are 2 lists available:
   - A collection of reconstruction CAL crystals (CalXtalRecData)
   - A collection of CAL clusters (CalCluster)
*/

class CalRecon : public TObject
{

public:

    CalRecon();

    virtual ~CalRecon();

    void Clear(Option_t *option="");

    void Print(Option_t *option="");

    /// provide access to list of CalClusters
    TObjArray* getCalClusters() { return m_clusterCol; };

    /// provide access to the list of CalXtalRecData
    std::vector<CalXtalRecData> getCalXtalRecCol() { return m_xtalRecCol; };

    void addCalCluster(CalCluster *cluster) { m_clusterCol->Add(cluster); };

    void addXtalRecData(const CalXtalRecData& xtal) { m_xtalRecCol.push_back(xtal); };

private:
    /// list of CalXtalRecData
//    TObjArray *m_xtalRecCol;      
    std::vector<CalXtalRecData> m_xtalRecCol;

    /// list of CalClusters
    TObjArray *m_clusterCol;  

  
    ClassDef(CalRecon,2)
};

#endif
