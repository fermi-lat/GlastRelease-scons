/** @file TkrVecPointInfo.h

* @class TkrVecPointInfo
*
* @brief This is the root version of the TkrVecPointInfo object. The purpose of the TkrVecPointInfo object is
*        to provide summary information about the TkrVecPoints (see TkrVecPoint.h)
*
* last modified 01/29/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TkrVecPointInfo_H
#define ROOT__TkrVecPointInfo_H

#include "TObject.h"
#include "TkrVecPoint.h"

class TkrVecPointInfo: public TObject
{
public:
    // constructors
    TkrVecPointInfo() : m_maxNumSkippedLayers(0),
                        m_numVecPoints(0),
                        m_numBiLayersWVecPoints(0),
                        m_maxNumLinkCombinations(0)
    {}

    TkrVecPointInfo(Int_t    maxNumSkippedLayers,
                    Int_t    numVecPoints,
                    Int_t    numBiLayersWVecPoints,
                    Double_t maxNumLinkCombinations) :
                    m_maxNumSkippedLayers(maxNumSkippedLayers),
                    m_numVecPoints(numVecPoints),
                    m_numBiLayersWVecPoints(numBiLayersWVecPoints),
                    m_maxNumLinkCombinations(maxNumLinkCombinations)
                    {}

    // destructor
    virtual ~TkrVecPointInfo() {}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrVecPointInfo&, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // initializer
    void inializeInfo(Int_t    maxNumSkippedLayers,
                      Int_t    numVecPoints,
                      Int_t    numBiLayersWVecPoints,
                      Double_t maxNumLinkCombinations)
    {
        m_maxNumSkippedLayers    = maxNumSkippedLayers;
        m_numVecPoints           = numVecPoints;
        m_numBiLayersWVecPoints  = numBiLayersWVecPoints;
        m_maxNumLinkCombinations = maxNumLinkCombinations;
    }

    /// @name Data set methods
    //@{
    /// Set the maximum number of skipped layers
    void setMaxNumSkippedLayers(int maxNumSkippedLayers)          {m_maxNumSkippedLayers    = maxNumSkippedLayers;}
    void setNumTkrVecPoints(int numVecPoints)                     {m_numVecPoints           = numVecPoints;}
    void setNumBiLayersWVecPoints(int numBiLayersWVecPoints)      {m_numBiLayersWVecPoints  = numBiLayersWVecPoints;}
    void setMaxNumLinkCombinations(double maxNumLinkCombinations) {m_maxNumLinkCombinations = maxNumLinkCombinations;}
    //@}

    /// @name Data access methods
    //@{
    /// How many layers are we allowed to skip when building links? 
    int    getMaxNumSkippedLayers()    {return m_maxNumSkippedLayers;}
    /// Total number of TkrVecPoints in the collection? 
    int    getNumTkrVecPoints()        {return m_numVecPoints;}
    /// Number of bilayers with TkrVecPoints
    int    getNumBiLayersWVecPoints()  {return m_numBiLayersWVecPoints;}
    /// Get the maximum number of link combinations which are possible
    double getMaxNumLinkCombinations() {return m_maxNumLinkCombinations;}
    //@}

private:

    // data members
    /// How many layers are we allowed to skip?
    Int_t      m_maxNumSkippedLayers;

    /// Also keep track of the total number of TkrVecPoints
    Int_t      m_numVecPoints;

    /// Keep count of the number of bilayers with TkrVecPoints
    Int_t      m_numBiLayersWVecPoints;
 
    /// Finally, keep track of max possible link combinations
    Double32_t m_maxNumLinkCombinations;

    ClassDef(TkrVecPointInfo,1)
};


#endif // __TkrVecPointInfo_H
