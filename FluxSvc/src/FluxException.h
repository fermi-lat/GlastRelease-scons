#ifndef FluxExecption_H_
#define FluxExecption_H_

/** 
* \class FluxException
*
* \brief Define FATAL MACRO, which throws exception with error message
* \author Toby Burnett
* 
* $Header $
*/

#include <strstream>


#define FATAL_MACRO(output)\
    std::cerr << output << std::endl;\
    do{std::ostrstream message; \
    message <<__FILE__<<":"<<__LINE__<<": "<<output<<'\0';\
throw(message.str());}while(0)

#define WARNING_MACRO(output)\
std::cerr << output << std::endl;

#endif // _ERROR_H_

