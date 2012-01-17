/** @file Medium.cxx
     @brief  implementation of class Medium

   $Header$
*/

#include "Medium.h"
#include "geometry/Box.h"
#include <cstring>

//////////////////////////////////////////////////////////////////////////////


Medium::Medium(Medium* parent, Shape* vol, const char* matName)
   :  _volume(vol)
   , _material(matName)
   ,  _parent(parent)
   ,  _title(0)
{
   if(_parent) {
      _parent->addMedium(this);
   }
   s_count ++;
}

Medium& Medium::setTitle(const char* newTitle)
{
   if( _title) delete [] _title;
   strcpy(_title = new char[strlen(newTitle)+1], newTitle);
   return *this;
}
Medium::~Medium()
{
   if(_title) delete [] _title;
   if(_volume) delete _volume;
    s_count--;
}

const char* Medium::nameOf() const {
  return "Medium";
}

int Medium::isComposite() const {
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//                manage inner media

Medium& Medium::addMedium ( Medium* nextMedium)
{
    //FATAL("Non-composite Medium tried to add medium");
    return *nextMedium;
}

Medium& Medium::removeMedium (Medium* oldMedium)
{
    //FATAL("Non-composite Medium tried to remove medium");
    return *oldMedium;
}

///////////////////////////////////////////////////////////////////////////////
//                       coordinate transformation
GeomObject&
Medium::transform(const CoordTransform& T)
{
    if( _volume )
	_volume->transform(T); // could happen during initialzation
   return *this;
}


//-----------------------------------------------------------------
//                  gobal statics
unsigned Medium::s_count=0;


