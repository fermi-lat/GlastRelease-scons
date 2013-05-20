/** @file TkrFilterParams.h

* @class TkrFilterParams
*
* @brief This is the root version of the TkrFilterParams object. The TkrFilterParams structure is the top
*        level object returned by the Hough Filter
*
* last modified 02/01/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TkrFilterParams_H
#define ROOT__TkrFilterParams_H

#include "TObject.h"
#include "reconRootData/TkrVecNodes.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

class TkrFilterParams: public TObject
{
public:
    // constructors - null case
    TkrFilterParams() {Clear();}

    // constructurs - have the information case
//    TkrFilterParams(TkrFilterParams* parent, const TkrVecPointsLink* associatedLink);

    // destructor
    virtual ~TkrFilterParams() {Clear();}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrFilterParams &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // Initializer
    void initializeInfo(TVector3& position, 
                        TVector3& direction,
                        UInt_t    statusBits,
                        Double_t  energy,
                        Int_t     numBiLayers,
                        Int_t     numIterations,
                        Int_t     numHitsTotal,
                        Int_t     numDropped,
                        Double_t  chiSquare,
                        Double_t  aveDist,
                        Double_t  transRms, 
                        Double_t  longRms,
                        Double_t  longRmsAsym
                        );

    /// Status word bits organized like:
    ///        |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <               > <               > <               >  <              >
    enum StatusBits {CALPARAMS  = 0x0001,  //Set if using Calorimeter parameters
                     TKRPARAMS  = 0x0002,  //Set if using Tracker parameters
                     FIRSTPASS  = 0x1000,  //Set if first pass numbers
                     SECONDPASS = 0x2000}; //Set if second pass numbers used

    /// Access data members
    UInt_t       getStatusBits()      const {return m_statusBits;}
    Double_t     getEventEnergy()     const {return m_eventEnergy;}
    TVector3     getEventPosition()   const {return m_eventPosition;}
    TVector3     getEventAxis()       const {return m_eventAxis;}

    Int_t        getNumBiLayers()     const {return m_numBiLayers;}
    Int_t        getNumIterations()   const {return m_numIterations;}
    Int_t        getNumHitsTotal()    const {return m_numHitsTotal;}
    Int_t        getNumDropped()      const {return m_numDropped;}

    Double_t     getChiSquare()       const {return m_chiSquare;}
    Double_t     getAverageDistance() const {return m_aveDist;}
    Double_t     getTransRms()        const {return m_transRms;}
    Double_t     getLongRms()         const {return m_longRms;}
    Double_t     getLongRmsAsym()     const {return m_longRmsAsym;}

    /// Modify data members
    void  setStatusBit(const UInt_t stat)        {m_statusBits    |= stat;}
    void  setEventEnergy(const Double_t energy)  {m_eventEnergy    = energy;}
    void  setEventPosition(const TVector3& pos)  {m_eventPosition  = pos;}
    void  setEventAxis(const TVector3& axis)     {m_eventAxis      = axis;}
    void  setNumBiLayers(UInt_t numBiLayers)     {m_numBiLayers    = numBiLayers;}
    void  setNumIterations(UInt_t numIterations) {m_numIterations  = numIterations;}
    void  setNumHitsTotal(UInt_t numHits)        {m_numHitsTotal   = numHits;}
    void  setNumDropped(UInt_t numDropped)       {m_numDropped     = numDropped;}
    void  setChiSquare(Double_t chiSquare)       {m_chiSquare      = chiSquare;}
    void  setAverageDistance(Double_t aveDist)   {m_aveDist        = aveDist;}
    void  setTransRms(Double_t transRms)         {m_transRms       = transRms;}
    void  setLongRms(Double_t longRms)           {m_longRms        = longRms;}
    void  setLongRmsAsym(Double_t longAsym)      {m_longRmsAsym    = longAsym;}

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);

private:
    UInt_t     m_statusBits;
    Double32_t m_eventEnergy;
    TVector3   m_eventPosition;
    TVector3   m_eventAxis;

    // Keep also some diagnostic information
    Int_t      m_numBiLayers;
    Int_t      m_numIterations;
    Int_t      m_numHitsTotal;
    Int_t      m_numDropped;

    // As well as moments output
    Double32_t m_chiSquare;
    Double32_t m_aveDist;
    Double32_t m_transRms; 
    Double32_t m_longRms;
    Double32_t m_longRmsAsym;


    ClassDef(TkrFilterParams,2)
};


#endif // __TkrFilterParams_H
