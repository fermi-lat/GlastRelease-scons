// File and Version Information:
// $Header$
// Description:

#include "AcdTileList.h"

IGeometry::VisitorRet 
AcdTileList::pushShape(ShapeType s, const UintVector& idvec, 
                       std::string name, std::string material, 
                       const DoubleVector& params, 
                       VolumeType type)
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
        this->push_back(getId());
        return AbortSubtree;
    }

    // otherwise continue
    return More;
}

void AcdTileList::popShape()
{
    unsigned int count = m_idcount.back(); m_idcount.pop_back();
    while(count--){
        m_idValues.pop_back();
    }
}

idents::VolumeIdentifier AcdTileList::getId()const
{
    idents::VolumeIdentifier id(m_prefix);
    for(UintVector::const_iterator i=m_idValues.begin(); i!=m_idValues.end(); ++i) id.append(*i);
    return id;

}
