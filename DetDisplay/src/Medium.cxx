// $Id$
//
//

#include "Medium.h"

#include "geometry/Box.h"

//////////////////////////////////////////////////////////////////////////////
//              constructors

Medium::Medium(Medium * prnt, float size)
   : _volume(new Box(size,size,size))
   , _material(0)
   , _detector(0)
   , _field(0)
   , _parent(prnt)
   , _title(0)
{
     set_defaults();
}
Medium::Medium(Medium* parent, Shape* vol, const char* matName, Detector* det)
   :  _volume(vol)
   , _material(0)
   ,  _detector(det)
   ,  _field(0)
   ,  _parent(parent)
   ,  _title(0)
{
	set_defaults();
}

void Medium::set_defaults()
{
   if(_parent) {
     // parent was specified: set appropiate attributes to be identical, then add

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

void Medium::createDetectorView(gui::DisplayRep& v) {

}

//-----------------------------------------------------------------
//                  gobal statics
const Medium * Medium::lastMedium = 0;
unsigned Medium::s_count=0;


