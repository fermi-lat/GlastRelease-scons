
#ifndef ROOT_OBJ_H
#define ROOT_OBJ_H

/*!

 @class RootObj        
 @brief Wrap any class into a kind of TObject

 This class is meant to transform any class T into a kind of TObject,
 so to benefit from some ROOT features, such as containers.

 For each kind of T, one is expected to add a pragma link in
 a LinkDef file, and it is also adviced to add an explicit
 template instanciation in some body file.

 @author David Chamont - CNRS IN2P3 LLR Ecole Polytechnique

 */

#include <TObject.h>

template <class T>
class RootObj : public TObject
 { 
  public:
    RootObj() ;
    RootObj( const T & ) ;
    T & value() ;
  private :
    T t_ ;
    ClassDef(RootObj,1)
 } ;


//===================================================================
// inline implementation
//===================================================================

template <class T>
RootObj<T>::RootObj()
 {
# ifdef DEBUG
//  rtb::out(1)<<"[RootObj::RootObj]"<<endl ;
# endif
 }
 
template <class T>
RootObj<T>::RootObj( const T & val )
 : t_(val)
 {
# ifdef DEBUG
//  rtb::out(1)<<"[RootObj::RootObj(const T &)] "<<val<<endl ;
# endif
 }
 
template <class T>
T & RootObj<T>::value()
 {
//# ifdef DEBUG
//  rtb::out(1)<<"[RootObj::operator T &] "<<endl ;
//# endif
  return t_ ;
 }

#endif
