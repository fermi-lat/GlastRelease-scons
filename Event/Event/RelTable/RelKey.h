
#ifndef RELKEY_H
#define RELKEY_H


#include "GaudiKernel/SmartRef.h"
#include <iostream>

/** 
 * @class RelKey
 *
 * @brief This class is used to relate events data.
 *
 * The RelKey class is used by the Relation class to relate events data. It
 * contains three pointers: one points to a particular event data, another
 * points to the next relation containing the same data pointer and the third
 * one points to the first relation non containing the same data pointer.  This
 * information is useful to efficiently search for all relations involving a
 * particular event data.
 * 
 *
 * @author Marco Frailis 
 * @author Riccardo Giannitrapani
 *    
 * $Header$
 */


namespace Event {

template <class T1, class T2>
class Relation;

template <class T1, class T2, class T3>
class RelKey {
    
public:
    
    RelKey() {}
    RelKey(T1* obj): m_data(obj) {}
    
    ~RelKey() {}
    
    void setData(T1* obj) {m_data = obj;}  
    const T1* getData() const { return m_data;}
    T1* getData() {return m_data;}

    void setSame(Relation<T2,T3>* rel) {m_same = rel;}
    Relation<T2,T3>* getSame() {return m_same;}

    void setFirst(Relation<T2,T3>* rel) {m_first = rel;}
    Relation<T2,T3>* getFirst() {return m_first;}

    /// Fill the ASCII output stream
    void toStream(std::ostream& s) const;

private:
    
    /// Pointer to the object to be related
    SmartRef<T1> m_data;
    /// Pointer to the next relation containing m_data
    SmartRef< Relation<T2,T3> > m_same;
    /// Pointer to the first relation which not contains m_data
    SmartRef< Relation<T2,T3> > m_first;
        
};


template <class T1, class T2, class T3>
inline void RelKey<T1,T2,T3>::toStream(std::ostream& s) const
{
  /// Fill the ASCII output stream
  s  << "\n        Data                    = " << m_data
     << "\n        Next Relation           = " << m_same
     << "\n        First Different Data    = " << m_first;
}

}

#endif // RELKEY_H 
