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

#include "ImPrecision.h"

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

// Define base class for storing in tuple map
class XTcolumnValBase
{
public:
    XTcolumnValBase(const std::string& name, const std::string& type) : m_name(name), m_type(type), m_valid(false) {}
    virtual ~XTcolumnValBase() {}

    const bool         dataIsValid() const {return m_valid;}
    const std::string& getType()     const {return m_type;}
    const std::string& getName()     const {return m_name;}

    void               clearValidFlag() {m_valid = false;}
    void               setValidFlag()   {m_valid = true;}
private:
    std::string m_name;
    std::string m_type;
    bool        m_valid;
};

// Typedef to define the overall container for our objects
typedef std::map<std::string,XTcolumnValBase* > XTtupleMap;

// Real stuff
template <class T> class XTcolumnVal : public XTcolumnValBase
{
public: 
    XTcolumnVal(const std::string& name, const std::string& type="continuous") : XTcolumnValBase(name,type) {}
    ~XTcolumnVal() {} //delete T;}

    const T*           operator()()  const
    {
        if (!dataIsValid())
        {
            throw XTexception("XTcolumnVal: Attempting to access invalid data: " + getName() + ", " + getType());
        }
        return &m_data;
    }
    
    T                  value()       const 
    {
        if (!dataIsValid())
        {
            throw XTexception("XTcolumnVal: Attempting to access invalid data: " + getName() + ", " + getType());
        }
        return m_data;
    }

    void setDataValue(const T& data) {m_data = data; setValidFlag();}
private:
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
