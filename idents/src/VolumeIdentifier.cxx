// $Header$

#include "idents/VolumeIdentifier.h"

#include <sstream>
#include <cassert>

using namespace idents;

VolumeIdentifier::VolumeIdentifier():  m_value(0){}

std::string VolumeIdentifier::name(const char* delimiter) const
{
    std::stringstream s;
    s << delimiter;
    std::copy(begin(),end(), std::ostream_iterator<unsigned int>(s,delimiter));
    std::string tmp=s.str();
    return tmp.substr(0,tmp.size()-1);
}
void VolumeIdentifier::prepend( const VolumeIdentifier& id)
{
    std::copy(id.rbegin(),id.rend(), std::front_inserter(*this));
    update();
}

void VolumeIdentifier::append( const VolumeIdentifier& id)
{
    std::copy(id.begin(), id.end(), std::back_inserter(*this));
    update();
}

void VolumeIdentifier::append( unsigned int id)
{
    push_back(id);
    update();
}

/*
std::ostream& operator<<(std::ostream& out, const VolumeIdentifier& vol)
{
std::copy(vol.begin(),vol.end(), std::ostream_iterator<unsigned int>(out,"/"));
return out;
}
*/
void VolumeIdentifier::update()
{
    m_value= 0;
    assert(size()<11); 
    for(VolumeIdentifier::iterator i = begin(); i!=end(); ++i){
        assert ((*i)<32);
        m_value =(m_value >> 6) | (*i) << 54;
    }
    
}