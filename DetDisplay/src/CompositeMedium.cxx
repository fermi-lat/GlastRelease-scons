// $Id$
// File: CompositeMedium.cxx


#include "CompositeMedium.h"

#include "gui/DisplayRep.h"

#include "geometry/Shape.h"
#include "geometry/Volume.h"
#include <algorithm>
//////////////////////////////////////////////////////////////////////////////
//              constructors

CompositeMedium::CompositeMedium(Medium * prnt, float size)
: Medium(prnt, size)
{}

CompositeMedium::CompositeMedium(Medium* prnt, Shape* vol, const char* matName, Detector* det)
: Medium(prnt, vol, matName, det)
{}

CompositeMedium::~CompositeMedium()
{
    deleteInnerMedia();
    if(_parent) {
        _parent->removeMedium(this);
    }
}
////////////////////////////////0/////////////////////////////////////////////
//                manage inner media
void CompositeMedium::deleteInnerMedia()
{
    for(reverse_iterator it=rbegin(); it !=rend(); ++it)
	delete *it;
}

Medium& CompositeMedium::addMedium ( Medium* nextMedium)
{
    push_back(nextMedium);
    nextMedium->setParent(this);
    return *nextMedium;
}

Medium& CompositeMedium::removeMedium (Medium* oldMedium)
{
    iterator it = std::find(begin(), end(), oldMedium);
    if( it != end())
	 erase(it);
    return *oldMedium;
}
const char* CompositeMedium::nameOf() const {
  return "CompositeMedium";
}

int CompositeMedium::isComposite() const {
  return 1;
}

///////////////////////////////////////////////////////////////////////////////
//                       coordinate transformation
GeomObject&
CompositeMedium::transform(const CoordTransform& T)
{
    Medium::transform(T);
    for(iterator it=begin(); it !=end(); ++it)
	(*it)->transform(T);
    return *this;
}
void CompositeMedium::createDetectorView(gui::DisplayRep& v)
{

    Medium::createDetectorView(v);
    for(iterator it=begin(); it !=end(); ++it )
       (*it)->createDetectorView(v.nested());
}

