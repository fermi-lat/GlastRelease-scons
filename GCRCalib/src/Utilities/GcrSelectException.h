
#ifndef __GcrSelectException_H
#define __GcrSelectException_H

#include <exception>
#include <string>

/**
@class GcrSelectException 
@brief hold a string
*/

class GcrSelectException : public std::exception {

public: 

    GcrSelectException( const std::string & error )
      : m_what(error) {
    }
    ~GcrSelectException() throw() {
    }
    virtual const char * what( ) const  throw() {
        return m_what.c_str() ;
    }
     
private:

    std::string m_what ;

} ;

#endif
