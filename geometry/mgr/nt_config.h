//This file copied from \packages\glast_settings\v0\src\nt_config.h 
// Include file for NT with Microsoft Visual Studio
// Sets some defines and some pragmas

#ifndef NT_CONFIG_H
#define NT_CONFIG_H

#define NT_MSVCPP
#define __STANDARD_CPLUSPLUS
#define HAVE_NEW_IOSTREAMS
#define HAVE_NAMESPACE_STD
#define DEFECT_NO_VIRTUAL_COVARIANCE
//#define DEFECT_NO_MEMBER_TEMPLATES
#define DEFECT_CMATH_NOT_STD
//#define NO_HEPTUPLE
//#define DEFECT_GUIDING_DECLS
#define DEFECT_MSVC
#define G__LITE


#define _WIN32_WINNT 0x0400 
#define DEFECT_STANDARD_CPLUSPLUS 
#define DEFECT_TYPE_INFO_GLOBAL 
#define DEFECT_CERNLIB_UPPERCASE 
#define DEFECT_NO_DELETE_CONST 
//#define DEFECT_CMATH_NOT_STD  
#define DATAREP_LITTLE_ENDIAN  
#define DS_LITTLE_ENDIAN
#define DSXIO 
// #define _MT
#define STLTYPE 1

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_2PI
#define M_2PI 6.28318530717958647692
#endif


// Work around for multiple definitions in for loops
#define for if(0);else for

#pragma warning( disable:4786 ) // Identifier is truncated to 255 characters in debug information
#pragma warning( disable:4244 ) // Conversion from 'unsigned short' to 'char'
#pragma warning( disable:4018 ) // Signed/unsigned mismatch
#pragma warning( disable:4800 ) // 'const int': forcing value to bool 'true' or 'false'
#pragma warning( disable:4355 ) // 'this': used in base member initializer list

#endif
