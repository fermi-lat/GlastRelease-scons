
#ifndef RELTABLE_H
#define RELTABLE_H


#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/SmartRef.h"
#include "Relation.h"
#include <vector>

/** 
 * @class RelTable
 *
 * @brief This class is used to wrap a collection of Relations.
 *
 * The RelTable class wraps a vector (a Gaudi ObjectList) of Relations. It
 * lets the user search for all object related to a given one. The search can be
 * done with respect to an object of the first or the second field of the
 * relations
 * 
 *
 * @author Marco Frailis
 * @author Riccardo Giannitrapani
 *   
 * $Header$
 */


namespace Event {

template <class T1, class T2>
class RelTable {
    
public:
    
    RelTable() {}
    RelTable(ObjectList < Relation<T1,T2> >* rels);
    

    /// Initialize the internal pointer to an ObjectList of relations
    void init() { m_relations = new ObjectList< Relation<T1,T2> >;}
      
    /// The following method add a new Relation to the vector of relations.
    void addRelation(Relation<T1,T2>* rel);

    /**
     * This method search for all relations having obj in the first
     * field. 
     * @param obj it's a pointer to the object given by the user
     * @return A vector of pointers to the relations involving the given object.
     */
    std::vector< Relation<T1,T2>* > getRelByFirst(T1* pobj);

    /**
     * This method search for all relations having pobj in the second
     * field. 
     * @param pobj it's a pointer to the object given by the user
     * @return A vector of pointers to the relations involving the given object.
     */
    std::vector< Relation<T1,T2>* > getRelBySecond(T2* pobj);
 
    /// This method returns the number of relations in the table
    unsigned long size();

    /// Returns the pointer to the collection of relations.
    ObjectList< Relation<T1,T2> >* getAllRelations();
       
private:
    
    /// Pointer to a collection of relations
    ObjectList < Relation<T1,T2> >* m_relations;
        
};




template <class T1,class T2>
inline RelTable<T1,T2>::RelTable(ObjectList < Relation<T1,T2> >* rels) {

  m_relations = rels;

}

template <class T1,class T2>
void RelTable<T1,T2>::addRelation(Relation<T1,T2>* rel) {
  // Purpose and Method:  This routine add a new relation to the collection.
  // Inputs:  rel is a pointer to the relation to be added.

  if (m_relations->size())
    {
      SmartRef< Relation<T1,T2> > r = m_relations->front();
      while ((r->getFirst() != rel->getFirst()))
        {
	  if (r->m_first.getFirst())
	    {
	      r = r->m_first.getFirst();
	    }
	  else
	    {
	      break;
	    }
        }
  
      if (r->getFirst() != rel->getFirst())
        {
          r->m_first.setFirst(rel);
        }
      else
        {
          rel->m_first.setSame(r->m_first.getSame());
          rel->m_first.setFirst(r->m_first.getFirst());
          r->m_first.setSame(rel);
        }

      r = m_relations->front();
      while ((r->getSecond() != rel->getSecond()))
        {
	  if (r->m_second.getFirst())
	    {
	      r = r->m_second.getFirst();
	    }
	  else
	    {
	      break;
	    }
        }
  
      if (r->getSecond() != rel->getSecond())
        {
          r->m_second.setFirst(rel);
        }
      else
        {
          rel->m_second.setSame(r->m_second.getSame());
          rel->m_second.setFirst(r->m_second.getFirst());
          r->m_second.setSame(rel);
        }
    }
  m_relations->push_back(rel);
}


template <class T1,class T2>
std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelByFirst(T1* pobj) {
  // Purpose and Method: This routine finds all relations having pobj in the
  // first field.  
  // Inputs: pobj is a pointer to the object to be searched in the first field
  // Outputs: A pointer to a vector of Relation* including pobj
  
  std::vector< Relation<T1,T2>* > rels;
  SmartRef< Relation<T1,T2> > r = m_relations->front();
  while (pobj != r->getFirst() && r->m_first.getFirst())
    {
      r = r->m_first.getFirst();
    }
  
  if (pobj == r->getFirst())
    {
      rels.push_back(r);
      while (r->m_first.getSame())
	      {
	        rels.push_back(r->m_first.getSame());
		r = r->m_first.getSame();
	      }
    }
  return rels;
}
  

template <class T1,class T2>
std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelBySecond(T2* pobj) {
  // Purpose and Method: This routine finds all relations having pobj in the
  // second field.  
  // Inputs: pobj is a pointer to the object to be searched in the second field
  // Outputs: A pointer to a vector of Relation* including pobj
  std::vector< Relation<T1,T2>* > rels;
  SmartRef< Relation<T1,T2> > r = m_relations->front();
  while (pobj != r->getSecond() && r->m_second.getFirst())
    {
      r = r->m_second.getFirst();
    }
  
  if (pobj == r->getSecond())
    {
      rels.push_back(r);
      while (r->m_second.getSame())
	      {
	        rels.push_back(r->m_second.getSame());
		r = r->m_second.getSame();
	      }
    }
  return rels;
}
 


template <class T1,class T2>
inline unsigned long RelTable<T1,T2>::size() {
  // Purposae: This method returns the total number of relations contained in the
  // collection
  return m_relations->size();
}

template <class T1,class T2>
inline ObjectList< Relation<T1,T2> >* RelTable<T1,T2>::getAllRelations() {

  return m_relations;

}

}

#endif // RELTABLE_H 
