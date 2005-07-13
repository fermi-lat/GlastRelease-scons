
#ifndef __CalException_H
#define __CalException_H

#include <exception>
#include <string>

/**
@class CalException 
@brief hold a string
*/

class CalException : public std::exception {

public: 

    CalException( const std::string & error )
      : m_what(error) {
    }
    ~CalException() throw() {
    }
    virtual const char * what( ) const  throw() {
        return m_what.c_str() ;
    }
     
private:

    std::string m_what ;

} ;

#endif
