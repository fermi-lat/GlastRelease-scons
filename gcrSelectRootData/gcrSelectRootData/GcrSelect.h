#ifndef ROOT_GCRSELECT_H
#define ROOT_GCRSELECT_H

#include "TObject.h"
#include "TObjArray.h"

//#include <vector>
//#ifdef WIN32
//using namespace std;
//#endif

//#include "GcrXtal.h"
#include "GcrSelectedXtal.h"
#include "GcrSelectVals.h"

/** @class CalRecon
 * @brief Primary Root object containing CAL reconstruction data.
 *
 * Portal through which one gains access to the CAL recon data. 
 * Currently there are 2 lists available:
 * - A collection of reconstructed CAL crystals (CalXtalRecData)
 * - A collection of CAL clusters (CalCluster)
 *
 * $Header$
*/

class GcrSelect : public TObject
{

public:

    GcrSelect();

    virtual ~GcrSelect();

    void initialize();

    void Clear(Option_t *option="");
    //void Fake( Int_t ievent, Float_t randNum ) ; // for tests
    //Bool_t CompareInRange( GcrSelect &, const std::string & name = "" ) ; // for tests
    void Print(Option_t *option="") const;
    
    TObjArray* getGcrSelectedXtalCol() { return m_gcrSelectedXtalCol; };
    TObject* getGcrSelectVals() { return m_gcrSelectVals; };
    
    /// C.L. 08/22/06:
    void addGcrSelectedXtal(GcrSelectedXtal* gcrSelectedXtal){ m_gcrSelectedXtalCol->Add(gcrSelectedXtal); };
    void addGcrSelectVals(GcrSelectVals* gcrSelectVals){ m_gcrSelectVals= gcrSelectVals; };
   
private:
    
    ///CL: 08/22/06 list of gcrXtals
    //TObjArray* m_gcrXtalCol;
    TObjArray* m_gcrSelectedXtalCol;
    TObject* m_gcrSelectVals;
    

    ClassDef(GcrSelect,5)
};

#endif
