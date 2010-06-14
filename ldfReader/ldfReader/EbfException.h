#ifndef EBFEXCEPTION_H
#define EBFEXCEPTION_H 1

/** @file EbfException.h
@brief A place to keep definitions of all ebf exception classes
$Header$
*/

class EbfExceptionBase {
public:
    EbfExceptionBase() {}
};

class EbfBadInstrument : public EbfExceptionBase {
public:
    EbfBadInstrument(const std::string& instrument) : m_name(instrument) {}
    std::string m_name;
};

#endif
