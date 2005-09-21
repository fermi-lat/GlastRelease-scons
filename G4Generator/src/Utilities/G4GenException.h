#ifndef G4GenException_H
#define G4GenException_H

#include <exception>
#include <string>

/**
@class G4GenException 
@brief hold a string
*/

class G4GenException : public std::exception 
{
public: 
    G4GenException( const std::string & error ) : m_what(error) {}
    ~G4GenException() throw() {}
    virtual const char * what( ) const  throw() {return m_what.c_str();}
     
private:
    std::string m_what ;
};

#endif
