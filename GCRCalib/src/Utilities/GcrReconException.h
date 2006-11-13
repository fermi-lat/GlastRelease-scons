
#ifndef __GcrReconException_H
#define __GcrReconException_H

#include <exception>
#include <string>

/**
@class GcrReconException 
@brief hold a string
*/

class GcrReconException : public std::exception {

public: 

    GcrReconException( const std::string & error )
      : m_what(error) {
    }
    ~GcrReconException() throw() {
    }
    virtual const char * what( ) const  throw() {
        return m_what.c_str() ;
    }
     
private:

    std::string m_what ;

} ;

#endif
