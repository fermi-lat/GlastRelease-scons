#ifndef RELTABLE_H
#define RELTABLE_H

#include "Relation.h"
#include <vector>
#include <algorithm>
#include <iterator>

/** 
* @class RelTable
*
* @brief This class is used to wrap a collection of Relations.
*
* The RelTable class wraps a list (a Gaudi ObjectList) of Relations. It
* lets the user search for all object related to a given one. The search can be
* done with respect to an object of the first or the second field of the
* relations. The user can also modify or delete relations.
* 
*
* @author Marco Frailis
* @author Riccardo Giannitrapani
*   
* $Header$
*/
namespace Event 
{

template <class T1, class T2> class RelTable 
{    
public: 
    RelTable() : m_relations(0), m_ownRelationList(true), m_firstMMap(0), m_secondMMap(0), m_ownMaps(true) {}
    RelTable(ObjectList < Relation<T1,T2> >* rels);
    RelTable(const RelTable<T1,T2>& table);

    // Destructor
    ~RelTable();
  
    /// Initialize the internal pointer to an ObjectList of relations
    void init(); // { m_relations = new std::list< Relation<T1,T2>* >;}

    /**
    * The following method add a Relation to the table if it doesn't contain
    * a relation between the same two objects, otherwise it appends the info
    * vector to the exsisting relation
    * @param rel is a pointer to a relation between two objects
    * @return true if the relation has been added and false if it is a duplicate
    * and has not been added (in this case the user has to delete it)
    */
    bool addRelation(Relation<T1,T2>* rel);
  
    /**
    * This method search for all relations having obj in the first
    * field. 
    * @param obj it's a pointer to the object given by the user
    * @return A vector of pointers to the relations involving the given object.
    */
    std::vector< Relation<T1,T2>* > getRelByFirst(const T1* pobj) const;
  
    /**
    * This method search for all relations having pobj in the second
    * field. 
    * @param pobj it's a pointer to the object given by the user
    * @return A vector of pointers to the relations involving the given object.
    */
    std::vector< Relation<T1,T2>* > getRelBySecond(const T2* pobj) const;
  
    /**
    * This method erase a particular relation from the table (keeping the 
    * integrity).
    * @param rel it's a pointer to the relation to be erased
    */
    void erase(Relation<T1,T2> *rel);
  
    /**
    * This method clears the Relation list and the multimaps
    */
    void clear();
  
    /**
    * This method change the first data pointer of a given relation contained
    * into the table.
    * @param rel it's a pointer to the relation to be modified
    * @param pobj is the new data value provided by the user
    */
    void changeFirst(Relation<T1,T2> *rel, T1 *pobj);
  
    /**
    * This method change the second data pointer of a given relation contained
    * into the table.
    * @param rel it's a pointer to the relation to be modified
    * @param pobj is the new data value provided by the user
    */
    void changeSecond(Relation<T1,T2> *rel, T2 *pobj);
  
    /// This method returns the number of relations in the table
    unsigned long size() const ;
  
    /// Returns the pointer to the collection of relations.
    RelationList<T1,T2>* getAllRelations(bool disown = true);
  
private:

    typedef typename RelationList<T1,T2>::iterator RelationListIter;

    typedef typename RelKeyMultiMap<T1,T1,T2>::iterator    mapT1RelIter;
    typedef typename std::pair<mapT1RelIter, mapT1RelIter> mapT1RelIterPair;

    typedef typename RelKeyMultiMap<T2,T1,T2>::iterator    mapT2RelIter;
    typedef typename std::pair<mapT2RelIter, mapT2RelIter> mapT2RelIterPair;
  
    /// Pointer to a collection of relations
    RelationList<T1,T2>*      m_relations;

    /// Flag to declare ownership of above list
    bool                      m_ownRelationList;

    /// Pointer to T1 multimap
    RelKeyMultiMap<T1,T1,T2>* m_firstMMap;

    /// Pointer to T2 multimap
    RelKeyMultiMap<T2,T1,T2>* m_secondMMap;

    /// Flag to declare ownership of maps
    bool                      m_ownMaps;
};
  
template <class T1, class T2> inline RelTable<T1,T2>::RelTable(ObjectList < Relation<T1,T2> >* rels) :
                              m_relations(0), m_ownRelationList(false), 
                                  m_firstMMap(0), m_secondMMap(0), m_ownMaps(true)
{
    // We are given a RelationList to work with
    m_relations  = dynamic_cast<Event::RelationList<T1,T2>*>(rels);

    // Will need to create and fill the multi-maps
    m_firstMMap  = new RelKeyMultiMap<T1,T1,T2>;
    m_secondMMap = new RelKeyMultiMap<T2,T1,T2>;

    // Set up and loop through the provided relations, resetting pointers to MMaps
    typename ObjectList<Event::Relation<T1,T2> >::iterator relIter;
    for(relIter = rels->begin(); relIter != rels->end(); relIter++)
    {
        Relation<T1,T2>* relation = *relIter;

//        relation->insertInList(m_relations);  // Not necessary to do this since we have set pointer to it
        relation->m_listIter = relIter;
        relation->insertFirst(m_firstMMap);
        relation->insertSecond(m_secondMMap);
    }
}

template <class T1, class T2> inline RelTable<T1,T2>::RelTable(const RelTable<T1,T2>& table)
{
    m_relations       = table.m_relations;
    m_ownRelationList = false;
    m_firstMMap       = table.m_firstMMap;
    m_secondMMap      = table.m_secondMMap;
    m_ownMaps         = false;
    return;
}

template <class T1, class T2> inline RelTable<T1,T2>::~RelTable()
{
    // If we don't "own" it then don't delete it (presumed in TDS)
    if (m_relations  && m_ownRelationList) delete m_relations;
    if (m_firstMMap  && m_ownMaps)         delete m_firstMMap;
    if (m_secondMMap && m_ownMaps)         delete m_secondMMap;
    return;
}

template <class T1, class T2> inline void RelTable<T1,T2>::init()
{
    if (m_relations)  delete m_relations;
    if (m_firstMMap)  delete m_firstMMap;
    if (m_secondMMap) delete m_secondMMap;

    m_ownRelationList = true;
    m_ownMaps         = true;

    m_relations  = new RelationList<T1,T2>;
    m_firstMMap  = new RelKeyMultiMap<T1,T1,T2>;
    m_secondMMap = new RelKeyMultiMap<T2,T1,T2>;
}

//
// Define a class for the sorting algorithm
// This will be used to sort a vector of pointers to TrackElements
//
template <class T1,class T2> class CompareRelations
{
public:
    CompareRelations(const Relation<T1,T2>* rel) : m_rel(rel) {}

    const bool operator()(const Relation<T1,T2>* rel) const
    {
        return *rel == *m_rel;
    }
private:
    const Relation<T1,T2>* m_rel;
};

template <class T1,class T2> bool RelTable<T1,T2>::addRelation(Relation<T1,T2>* rel) 
{
    // Purpose and Method:  This routine add a relation to the table if it doesn't 
    // contain a relation between the same two objects, otherwise it appends the info
    // vector to the exsisting relation
    // Inputs:  rel is a pointer to the relation to be added.
    // Outputs: a boolean value which is true if the realtion has been added to the
    //          table and false it it is a duplicate and thus has not been added.
    //          In the latter case the user has to delete the relation
    bool addedToList = false;

    // First search to see if this relation is already in the list
    Relation<T1,T2>* locRel = 0;

    // Use MMap to get vector of relations to test against
    mapT1RelIterPair iterPair = m_firstMMap->equal_range(rel->getFirst());

    // Loop through to find match on second
    for(mapT1RelIter mapIter = iterPair.first; mapIter != iterPair.second; mapIter++)
    {
        typename RelationList<T1,T2>::RelationListIter relIter = (*mapIter).second;
        Relation<T1,T2>* relation = *relIter;
        if (relation->getSecond() == rel->getSecond()) locRel = relation;
    }

    // If not in list then add it
    if (!locRel)
    {
        // This adds the relation to the list
        rel->insertInList(m_relations);

        // Take care of the maps
        rel->insertFirst(m_firstMMap);
        rel->insertSecond(m_secondMMap);

        addedToList = true;
    }
    // Otherwise, update the "info" vector 
    else if (!rel->getInfos().empty())
    {
        locRel->addInfo(rel->getInfos().front());
    }

    return addedToList;
}

template <class T1,class T2>
    std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelByFirst(const T1* pobj) const 
{
    // Purpose and Method: This routine finds all relations having pobj in the
    // first field.  
    // Inputs: pobj is a pointer to the object to be searched in the first field
    // Outputs: A pointer to a vector of Relation* including pobj
    
    std::vector< Relation<T1,T2>* > rels;
    if (!m_relations->size()) return rels;

    T1* pObjLocal = const_cast<T1*>(pobj);

    mapT1RelIterPair iterPair = m_firstMMap->equal_range(pObjLocal);

    for(mapT1RelIter mapIter = iterPair.first; mapIter != iterPair.second; mapIter++)
    {
        typename RelationList<T1,T2>::RelationListIter relIter = (*mapIter).second;
        Relation<T1,T2>* relation = *relIter;
        rels.push_back(relation);
    }

    return rels;
} 
  
template <class T1,class T2>
    std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelBySecond(const T2* pobj) const 
{
    // Purpose and Method: This routine finds all relations having pobj in the
    // second field.  
    // Inputs: pobj is a pointer to the object to be searched in the second field
    // Outputs: A pointer to a vector of Relation* including pobj
    std::vector< Relation<T1,T2>* > rels;
    if (!m_relations->size()) return rels;

    T2* pObjLocal = const_cast<T2*>(pobj);

    mapT2RelIterPair iterPair = m_secondMMap->equal_range(pObjLocal);

    for(mapT2RelIter mapIter = iterPair.first; mapIter != iterPair.second; mapIter++)
    {
        typename RelationList<T1,T2>::RelationListIter 
          relIter = (*mapIter).second;
        Relation<T1,T2>* relation = *relIter;
        rels.push_back(relation);
    }

    return rels;
}
  
template <class T1,class T2> void RelTable<T1,T2>::erase(Relation<T1,T2>* rel) 
{
    // Purpose: This method remove the given relation from the table

    rel->removeFirst(m_firstMMap);
    rel->removeSecond(m_secondMMap);
    rel->removeFromList(m_relations);

    delete rel;
}
  
template <class T1,class T2>
    void RelTable<T1,T2>::changeFirst(Relation<T1,T2> *rel, T1 *pobj) 
{
    // Purpose: This method change the first data pointer of a relation with the
    // one given by the user

    T2* pSecond = rel->getSecond();
    std::vector<std::string> info = rel->getInfos();

    erase(rel);

    rel = new Event::Relation<T1,T2>(pobj, pSecond, info);

    addRelation(rel);
}


template <class T1,class T2>
    void RelTable<T1,T2>::changeSecond(Relation<T1,T2> *rel, T2 *pobj) 
{
    // Purpose: This method change the second data pointer of a relation with the
    // one given by the user

    T1* pFirst = rel->getFirst();
    std::vector<std::string> info = rel->getInfos();

    erase(rel);

    rel = new Event::Relation<T1,T2>(pFirst, pobj, info);

    addRelation(rel);
}
  
template <class T1,class T2> void RelTable<T1,T2>::clear() 
{
    // Purpose: This method removes all relations from the table

    m_firstMMap->clear();
    m_secondMMap->clear();

    // Note that this calls ObjectList's clear which will delete
    // the relations... 
    m_relations->clear();

    return;
}


template <class T1,class T2>
    inline unsigned long RelTable<T1,T2>::size() const 
{
    // Purpose: This method returns the total number of relations contained in the
    // collection
    return m_relations->size();
}
  
template <class T1,class T2> inline RelationList<T1,T2>* RelTable<T1,T2>::getAllRelations(bool disown)
{
    // Does caller want ownership of the list (default mode)
    if (disown) m_ownRelationList = false;

    return m_relations;
}

}; // Namespace
#endif // RELTABLE_H 
