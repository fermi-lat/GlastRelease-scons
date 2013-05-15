/** @file TkrVecPointsLink.h

* @class TkrVecPointsLink
*
* @brief This is the root version of the TkrVecPointsLink object. The purpose of the TkrVecPointsLink object is
*        to represent a 3D "point" in the Tracker definied by the crossing of two TkrClusters
*
* last modified 01/29/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TkrVecPointsLink_H
#define ROOT__TkrVecPointsLink_H

#include "TObject.h"
#include "commonRootData/idents/TkrId.h"
#include "TkrVecPoint.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

class TkrVecPointsLink: public TObject
{
public:
    // Status bit definitions    
    enum StatusBits {ASSOCIATED   = 0x00000001,
                     FIRSTLINK    = 0x00000002,
                     SKIP1LAYER   = 0x00000010,
                     SKIP2LAYER   = 0x00000020,
                     SKIP3LAYER   = 0x00000040,
                     SKIPNLAYER   = 0x00000080,
                     SAMETOWER    = 0x00000100,
                     VERIFIED     = 0x00000800,
                     INTERTOWER   = 0x00010000,
                     WAFERGAP     = 0x00020000,
                     WAFERGAPPLUS = 0x00040000,
                     GAPANDCLUS   = 0x00080000,
                     TRUNCATED    = 0x00100000,
                     BADSTRIPS    = 0x00080000
                                         };

    enum skipsBits  {SKIPSLAYERS = SKIP1LAYER | SKIP2LAYER | SKIP3LAYER | SKIPNLAYER};

    // constructors - null case
    TkrVecPointsLink() {Clear();}

    // constructurs - have the information case
    TkrVecPointsLink(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, double ang);

    // destructor
    virtual ~TkrVecPointsLink() {Clear();}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrVecPointsLink &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // Initializer
    void initializeInfo(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, UInt_t status, double ang);

    void setUnAssociated()                {m_statusBits      &= ~ASSOCIATED;}
    void setAssociated()                  {m_statusBits      |=  ASSOCIATED;}
    void setNotFirstLink()                {m_statusBits      &= ~FIRSTLINK;}
    void setFirstLink()                   {m_statusBits      |=  FIRSTLINK;}
    void setSameTower()                   {m_statusBits      |=  SAMETOWER;}
    void setSkip1Layer()                  {m_statusBits      |=  SKIP1LAYER;}
    void setSkip2Layer()                  {m_statusBits      |=  SKIP2LAYER;}
    void setSkip3Layer()                  {m_statusBits      |=  SKIP3LAYER;}
    void setSkipNLayer()                  {m_statusBits      |=  SKIPNLAYER;}
    void setVerified()                    {m_statusBits      |=  VERIFIED;}
    void updateStatusBits(unsigned bits)  {m_statusBits      |=  bits;}
    void clearStatusBits(unsigned bits)   {m_statusBits      &= ~bits;}
    void setMaxScatAngle(double ang)      {m_maxScatAngle     =  ang;}

    const bool         associated()          const {return (m_statusBits & ASSOCIATED)  == ASSOCIATED;}
    const bool         firstLink()           const {return (m_statusBits & FIRSTLINK)   == FIRSTLINK;}
    const bool         sameTower()           const {return (m_statusBits & SAMETOWER)   == SAMETOWER;}
    const bool         skipsLayers()         const {return (m_statusBits & SKIPSLAYERS) != 0;}
    const bool         skip1Layer()          const {return (m_statusBits & SKIP1LAYER)  == SKIP1LAYER;}
    const bool         skip2Layer()          const {return (m_statusBits & SKIP2LAYER)  == SKIP2LAYER;}
    const bool         skip3Layer()          const {return (m_statusBits & SKIP3LAYER)  == SKIP3LAYER;}
    const bool         skipNLayer()          const {return (m_statusBits & SKIPNLAYER)  == SKIPNLAYER;}
    const bool         verified()            const {return (m_statusBits & VERIFIED)    == VERIFIED;}
    const unsigned int getStatusBits()       const {return m_statusBits;}
    const TVector3&    getPosition()         const {return m_position;}
    const TVector3     getPosition(double z) const;
    const TVector3     getBotPosition()      const {return getPosition(m_secondVecPoint->getPosition().z());}
    const TVector3&    getVector()           const {return m_vector;}
    const TkrVecPoint* getFirstVecPoint()    const {return m_firstVecPoint;}
    const TkrVecPoint* getSecondVecPoint()   const {return m_secondVecPoint;}
    const Double_t     getMaxScatAngle()     const {return m_maxScatAngle;}

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);
    

private:
    const TkrVecPoint* m_firstVecPoint;
    const TkrVecPoint* m_secondVecPoint;
    
    UInt_t             m_statusBits;
    TVector3           m_position;
    TVector3           m_vector;

    // Calculated expected maximum scattering angle
    Double32_t         m_maxScatAngle;

    ClassDef(TkrVecPointsLink,1)
};


#endif // __TkrVecPointsLink_H
