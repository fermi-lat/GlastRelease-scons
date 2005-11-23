/**@file XTtupleVars.h
@brief Contains class definitions for implementing a "local" ntuple  row
@author T. Usher
$Header$
*/

#ifndef XTtupleVars_h
#define XTtupleVars_h

#include "GlastClassify/ITupleInterface.h"

#include <iostream>
#include <map>

/** @class Exception 
    @brief hold a string
*/
class XTexception : public std::exception
{
public: 
    XTexception(std::string error):m_what(error){}
    ~XTexception() throw() {;}
    virtual const char *what( ) const  throw() { return m_what.c_str();} 
    std::string m_what;
};

// Forward declaration for compiler issues
template <class T> class XTcolumnVal;
template <class T> std::ostream& operator <<(std::ostream& stream, const XTcolumnVal<T>& tupleVal);

// Real stuff
template <class T> class XTcolumnVal
{
public:
    typedef typename std::map<std::string,XTcolumnVal<T>* > XTtupleMap;
    
    XTcolumnVal(const std::string& name) : m_name(name), m_valid(false), m_data(0.) {}
    ~XTcolumnVal() {} //delete T;}

    const bool         dataIsValid() const {return m_valid;}
    const std::string& getName()     const {return m_name;}

    const T*           operator()()  const
    {
        if (!m_valid)
        {
            throw XTexception("XTcolumnVal: Attempting to access invalid data");
        }
        return &m_data;
    }
    
    T                  value()       const 
    {
        if (!m_valid)
        {
            throw XTexception("XTcolumnVal: Attempting to access invalid data");
        }
        return m_data;
    }

    void clearValidFlag() {m_valid = false;}
    void setDataValue(const T data) {m_data = data; m_valid = true;}
private:
    std::string m_name;
    bool        m_valid;
    T           m_data;

    // This is for making a fancy output... I'm not necessarily proud of it...
    friend std::ostream& operator <<<T>(std::ostream& stream, const XTcolumnVal<T>& node);
};

template <class T> std::ostream& operator <<(std::ostream& stream, const XTcolumnVal<T>& tupleVal)
{
    stream << tupleVal.m_name << "+" << tupleVal.m_data;

    return stream;
}

#endif
