#ifndef CALRECON_H
#define CALRECON_H

#include "TObject.h"
#include "TObjArray.h"

// ========================================================
//      CalRecon
// ========================================================
//! Primary Root object containing CAL reconstruction data.
/*! Primary CAL reconstruction Data Root object
portal through which one gains access to the CAL recon data. 
Currently there are 2 lists available:
1.) A list of hit CsI Logs (CalLogEne)
2.) A list of CAL clusters (CalCluster)
*/

class CalRecon : public TObject
{

private:
    //! list of hit logs
    TObjArray *m_logs;      

    //! list of CalClusters
    TObjArray *m_clusters;  

public:

    //! constructor
    // default constructor
    // so client must call Create() to setup 
    // TObjArray properly
    CalRecon();

    //! destructor
    virtual ~CalRecon();

    //! frees all dynamically allocated memory
    void Clean();

    //! allocates memory as necessary
    void Create();

    //! provide access to list of CalClusters
    TObjArray* getCalClusters() { return m_clusters; };

    //! provide access to the list of CalLogs
    TObjArray* getCalLogs() { return m_logs; };
  
    ClassDef(CalRecon,1)
};

#endif
