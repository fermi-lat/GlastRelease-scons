#ifndef LdfException_H
#define LdfException_H 1

/** @file LdfException.h
@brief A place to keep definitions of all ebf exception classes
$Header$
*/

class LdfExceptionBase {
public:
    LdfExceptionBase() {}
};

class LdfBadInstrument : public LdfExceptionBase {
public:
    LdfBadInstrument(const std::string& instrument) : m_name(instrument) {}
    std::string m_name;
};

#endif
