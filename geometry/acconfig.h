/* 
 * Macross defining non-supported C++ language features and compiler
 * specific pragmas.  
 *
 * These should be replaced with a autoconf generated define.
 *
 * $Id$
 *
 * Author: Paul_Kunz@slac.stanford.edu
 * */

#ifndef GLASTSIM_CONFIG_H
#define GLASTSIM_CONFIG_H

#ifdef _MSC_VER
/* define if max() and min() template functions come with STL */
#undef CLHEP_MAX_MIN_DEFINED
# pragma warning(disable:4786)  // '255' characters in the debug information
# pragma warning(disable:4305)  // conversion from double to float
#else
#define CLHEP_MAX_MIN_DEFINED
#endif //_MSC_VER

@TOP@
@BOTTOM@

#endif /* GLASTSIM_CONFIG_H */

