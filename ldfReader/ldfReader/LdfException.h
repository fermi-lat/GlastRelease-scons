#ifndef LdfException_H
#define LdfException_H 1

/** @file LdfException.h
@brief A place to keep definitions of all ebf exception classes
$Header$
*/

class LdfException : public std::exception {
public:
    LdfException(std::string error):m_what(error) {}
    ~LdfException() throw() {;}
    virtual const char *what() const throw() { return m_what.c_str(); }
    std::string m_what;
};

//class LdfBadInstrument : public LdfExceptionBase {
//public:
//    LdfBadInstrument(const std::string& instrument) : m_name(instrument) {}
//    std::string m_name;
//};

#endif
