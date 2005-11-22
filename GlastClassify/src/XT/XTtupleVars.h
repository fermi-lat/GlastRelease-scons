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
    friend std::ostream& operator <<(std::ostream& stream, const XTcolumnVal<T>& node);
};

template <class T> std::ostream& operator <<(std::ostream& stream, const XTcolumnVal<T>& tupleVal)
{
    stream << tupleVal.m_name << "+" << tupleVal.m_data;

    return stream;
}

/*
// Define the class to contain the above objects
template <class T> class XTtupleVars
{
public:
    typedef typename std::map<std::string,XTcolumnVal<T>* > XTtupleMap;

    XTtupleVars(const GlastClassify::ITupleInterface& tuple) : m_iTuple(tuple) {m_tupleMap.clear();}
    ~XTtupleVars() {}

    inline void initEvent();

    inline XTcolumnVal<T>* getColumnVal(const std::string& name);

    inline XTcolumnVal<T>* addNewDataItem(const std::string& name); 
    
    inline void print(std::ostream& out=std::cout) const;

private:
    const GlastClassify::ITupleInterface& m_iTuple;
    XTtupleMap                            m_tupleMap;
};

template <class T> inline void XTtupleVars<T>::initEvent() 
{
    // Loop over entries in the map
    for(XTtupleMap::iterator dataIter = m_tupleMap.begin(); dataIter != m_tupleMap.end(); dataIter++)
    {
        XTcolumnVal<T>*    colVal = dataIter->second;
        const std::string& name   = colVal->getName();

        // Use a "try" to catch the exception if not in the tuple
        try
        {
            const Item* item = m_iTuple.getItem(name);

            if (item != 0)
            {
                double value = (*item);
                colVal->setDataValue(value);
            }
            else
            {
                // Set the "value"
                colVal->setDataValue(0.);
                colVal->clearValidFlag();
            }
        }
        catch (std::invalid_argument&)
        {
            int j = 0;
        }
    }
    return;
}

template <class T> inline XTcolumnVal<T>* XTtupleVars<T>::getColumnVal(const std::string& name)
{
    XTcolumnVal<T>* dataPtr = 0;

    XTtupleMap::iterator dataIter = m_tupleMap.find(name);

    if (dataIter != m_tupleMap.end())
    {
        dataPtr = dataIter->second;
    }

    return dataPtr;
}

template <class T> inline XTcolumnVal<T>*  XTtupleVars<T>::addNewDataItem(const std::string& name) 
{
    // Let's make sure the item doesn't already exist
    XTcolumnVal<T>* colValPtr = getColumnVal(name);

    if (colValPtr == 0)
    {
        colValPtr = new XTcolumnVal<T>(name);
    
        m_tupleMap[name] = colValPtr;
    }

    return colValPtr;
}
    
template <class T> inline void XTtupleVars<T>::print(std::ostream& out=std::cout) const
{
    int numVars = m_tupleMap.size();
    out << "Local Tuple Map size: " << numVars << std::endl;

    for(XTtupleMap::const_iterator dataIter = m_tupleMap.begin(); dataIter != m_tupleMap.end(); dataIter++)
    {
        out << dataIter->first << ",  " << dataIter->second->getName() << std::endl;
    }
    return;
}
*/
#endif
