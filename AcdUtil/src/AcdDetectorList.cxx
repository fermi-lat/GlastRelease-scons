// File and Version Information:
// $Header$
// Description:

#include "AcdUtil/AcdDetectorList.h"
#include <iostream>

namespace AcdUtil {

IGeometry::VisitorRet 
AcdDetectorList::pushShape(ShapeType s, const UintVector& idvec, 
                       std::string name, std::string material, 
                       const DoubleVector& params, 
                       VolumeType type, SenseType sense)
{
    // concatenate the id for this node to current id.
    m_idcount.push_back(idvec.size());
    for( UintVector::const_iterator u=idvec.begin(); u!=idvec.end(); ++u){
        m_idValues.push_back(static_cast<unsigned int>(*u));
    }

    // is this what we want? add to the list if so and abort
    if (name.substr(0,7)=="topTile") {
        this->push_back(getId());
        return AbortSubtree;
    } else if (name.substr(0,8)=="sideTile" ) {
        if (name.substr(0,11) == "sideTileRow") return More;
        if (name.substr(0,10) == "sideTileR3") {
            this->push_back(getId());
            return More;
        }
        this->push_back(getId());
        return AbortSubtree;
    } else if (name.substr(0,10) == "ACDScrewSq") {
        this->push_back(getId());
        return AbortSubtree; 
    } else if (name.substr(0,10) == "sideRibbon" ) {
        if (name.substr(0,11) == "sideRibbons" ) return More;
        // ignore top ribbons - we just want a count of whole ribbons not the segments
        // since there are 2 SideRibbons per ribbon - 
        // will need to check to see if we found this guy already

        // Also check to see if the ribbons are position detectors or not
        if (sense == posSensitive) this->push_back(getId());
        return AbortSubtree;
    }

    // otherwise continue
    return More;
}

void AcdDetectorList::popShape()
{
    unsigned int count = m_idcount.back(); m_idcount.pop_back();
    while(count--){
        m_idValues.pop_back();
    }
}

idents::VolumeIdentifier AcdDetectorList::getId()const
{
    idents::VolumeIdentifier id(m_prefix);
    for(UintVector::const_iterator i=m_idValues.begin(); i!=m_idValues.end(); ++i) id.append(*i);
    return id;

}

}
