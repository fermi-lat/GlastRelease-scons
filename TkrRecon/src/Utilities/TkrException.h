/**
 * @file TkrException.h
 *

 * $Header$
 *
*/

#ifndef __TkrException_H
#define __TkrException_H

    /** @class TkrException 
         @brief hold a string
         */
class TkrException : public std::exception
{
public: 
    TkrException(std::string error):m_what(error){}
    ~TkrException() throw() {;}
        
    virtual const char *what( ) const  throw() { return m_what.c_str();} 
        
    std::string m_what;
};

#endif
