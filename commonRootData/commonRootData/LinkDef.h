/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RelTable+;
#pragma link C++ class Relation+;
#pragma link C++ class BaseException;
#pragma link C++ class IndexOutOfBoundsException;
//#pragma link C++ class BadPointerException;

#endif
