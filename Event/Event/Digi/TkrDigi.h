#ifndef GlastEvent_TkrDigi_H
#define GlastEvent_TkrDigi_H 1

//include files

#include <iostream>
#include <vector>

#include "Gaudi/Kernel/Kernel.h"
#include "Gaudi/Kernel/ContainedObject.h"
#include "Gaudi/Kernel/SmartRefVector.h"

#include "GlastEvent/TopLevel/Definitions.h"
#include "GlastEvent/TopLevel/ObjectVector.h"
#include "GlastEvent/TopLevel/ObjectList.h"

/*!
* \class TkrDigi
* \author Leon Rochester
* \brief TDS version of the TKR digi
*
* Contains layer and tower identification, ToT and a list of hit strips
* The digis are produced either from MC hit output or from the actual data
*/

extern const CLID& CLID_TkrDigi;

class TkrDigi : virtual public ContainedObject {
    
public:
    //! typedefs
    typedef std::vector<int> HitList;
    typedef HitList::const_iterator const_iterator;
    
    //! Constructors
    //! Null constructor
    TkrDigi() {};
    
    //! constructor with layer, tower and ToT, but with an empty list of strips
    TkrDigi(int l, int v, int t, int* tot)
        : m_layer (l), m_view (v), m_tower (t) {
        m_tot[0] = *tot;
        m_tot[1] = *(++tot);
    };
    //! Destructor
    virtual ~TkrDigi() {
    };
    
    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID() const   { return TkrDigi::classID(); }
    static const CLID& classID()       { return CLID_TkrDigi; }
    
    //! Retrieve methods 
    //! Retrieve layer info
    int layer() const {return m_layer;} 
    //! Retrieve view info
    int view () const {return m_view;}
    //! Retrieve tower info
    int tower() const {return m_tower;}
    //! Retrieve ToT info
    int ToT( int i) const {return m_tot[i];}
    //! number of hits
    int num() const{ return m_hits.size(); }
    //! Get the pointer to the ith hit strip
    int  hit( int i ) const { return m_hits[i];}
    
    //! Set Methods
    //! Set layer
    void setLayer( int layer) {m_layer = layer;}
    //! set view
    void setView( int view) {m_view = view;}
    //! Set tower
    void setTower( int tower ) {m_tower = tower;}
    //! Set the ToT array
    void setToT ( int i, int tot) {m_tot[i] = tot;}
    //! Add a hit to the hit list
    void addHit( int strip ) {m_hits.push_back(strip);}
    
    //! Serialize the object for reading
    virtual StreamBuffer& serialize( StreamBuffer& s );
    //! Serialize the object for writing
    virtual StreamBuffer& serialize( StreamBuffer& s ) const;
    //! Fill the ASCII output stream
    virtual std::ostream& fillStream( std::ostream& s ) const;
    
    //! begin iterator
    const_iterator begin()const {return m_hits.begin();}
    //! end iterator
    const_iterator end()const {return m_hits.end();}
    
private:
    
    int m_layer;
    int m_view;
    int m_tower;
    int m_tot[2];
    HitList m_hits;
};


//! Serialize the object for writing
inline StreamBuffer& TkrDigi::serialize( StreamBuffer& s ) const {
    ContainedObject::serialize(s);  
    s << m_layer
        << m_view
        << m_tower
        << m_tot[0]
        << m_tot[1]
        << m_hits.size();
    const_iterator ih;
    for (ih = m_hits.begin(); ih!=m_hits.end(); ih++) {
        s << *ih;
    }
    
    return s;
}

//! Serialize the object for reading
inline StreamBuffer& TkrDigi::serialize( StreamBuffer& s )       {
    ContainedObject::serialize(s);
    int size;
    s >> m_layer
        >> m_view
        >> m_tower
        >> m_tot[0]
        >> m_tot[1]
        >> size;
    
    m_hits.resize( size, 0);
    std::vector<int>::iterator ih;
    for (ih = m_hits.begin(); ih!=m_hits.end(); ih++) {
        s >> *ih;
    }
    
    return s;
}

//! Fill the ASCII output stream

inline std::ostream& TkrDigi::fillStream( std::ostream& s ) const {
    int j;
    int size = m_hits.size();
    s << "class TkrDigi :" << std::endl
        << "Layer: " << m_layer 
        << " view: " << m_view
        << " tower: " << m_tower
        << " ToT: " << m_tot[0] << " " << m_tot[1] << std::endl
        << "Number of hits strips: " << size << std::endl;
    
    const_iterator ih;
    for(ih = m_hits.begin(), j=0; ih != m_hits.end();ih++,j++) {
        if (j==10) {j = 0; s << std::endl;}
        s << *ih << " ";
    }
    s << std::endl;
    return s;
}

//! Definition of all container types of TkrDigi
/*! These are the containers that will go into the TDS.
*/

template <class TYPE> class ObjectVector;
template <class TYPE> class ObjectList;

typedef ObjectVector<TkrDigi> TkrDigiCol;

#endif
