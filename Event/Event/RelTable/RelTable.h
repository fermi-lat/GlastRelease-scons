
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
  ObjectList< Relation<T1,T2> >* getAllRelations() const;
  
private:
  
  /// Pointer to a collection of relations
  ObjectList < Relation<T1,T2> >* m_relations;
  

  void bindRelationFirst(Relation<T1,T2> *rel);
  void bindRelationSecond(Relation<T1,T2> *rel);

  void removeFirst(Relation<T1,T2> *rel);
  void removeSecond(Relation<T1,T2> *rel);


  };
  
  
  
  
  template <class T1,class T2>
    inline RelTable<T1,T2>::RelTable(ObjectList < Relation<T1,T2> >* rels) {
    
    m_relations = rels;
    
  }
  
  template <class T1,class T2>
    void RelTable<T1,T2>::addRelation(Relation<T1,T2>* rel) {
    // Purpose and Method:  This routine add a new relation to the collection.
    // Inputs:  rel is a pointer to the relation to be added.
    
    bindRelationFirst(rel);
    bindRelationSecond(rel);
    m_relations->push_back(rel);
  }
  
  
  template <class T1,class T2>
    std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelByFirst(const T1* pobj) const {
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
    std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelBySecond(const T2* pobj) const {
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
    void RelTable<T1,T2>::erase(Relation<T1,T2> *rel) {
    // Purpose: This method remove the given relation from the table

    removeFirst(rel);
    removeSecond(rel);

    m_relations->remove(rel);
    delete rel;
  }

  
  template <class T1,class T2>
    void RelTable<T1,T2>::changeFirst(Relation<T1,T2> *rel, T1 *pobj) {
    // Purpose: This method change the first data pointer of a relation with the
    // one given by the user

    removeFirst(rel);    
    rel->setFirst(pobj);
    bindRelationFirst(rel);
  }


  template <class T1,class T2>
    void RelTable<T1,T2>::changeSecond(Relation<T1,T2> *rel, T2 *pobj) {
    // Purpose: This method change the second data pointer of a relation with the
    // one given by the user

    removeSecond(rel);    
    rel->setSecond(pobj);
    bindRelationSecond(rel);
  }


  template <class T1,class T2>
    inline unsigned long RelTable<T1,T2>::size() const {
    // Purpose: This method returns the total number of relations contained in the
    // collection
    return m_relations->size();
  }
  
  template <class T1,class T2>
    inline ObjectList< Relation<T1,T2> >* RelTable<T1,T2>::getAllRelations() const {
    
    return m_relations;
    
  }

  

  template <class T1,class T2>
    inline void RelTable<T1,T2>::bindRelationFirst(Relation<T1,T2> *rel) {
    
    Relation<T1,T2>* temp;

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
        rel->m_first.setPrev(r);
      }
      else
      {
        temp = r->m_first.getSame();
        rel->m_first.setSame(temp);
        if (temp)
          temp->m_first.setPrev(rel);
        r->m_first.setSame(rel);
        rel->m_first.setPrev(r);
      }
    }
  }


  template <class T1,class T2>
    inline void RelTable<T1,T2>::bindRelationSecond(Relation<T1,T2> *rel) {
    Relation<T1,T2>* temp;

    if (m_relations->size())
    {
      SmartRef< Relation<T1,T2> > r = m_relations->front();

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
        rel->m_second.setPrev(r);
      }
      else
      {
        temp = r->m_second.getSame();
        rel->m_second.setSame(temp);
        if (temp) 
          temp->m_second.setPrev(rel);
        r->m_second.setSame(rel);
        rel->m_second.setPrev(r);
      }
    }    
  }


  template <class T1,class T2>
    inline void RelTable<T1,T2>::removeFirst(Relation<T1,T2> *rel) {

    Relation<T1,T2> *temp;

    temp = rel->m_first.getSame();
    if (temp)
      temp->m_first.setPrev(rel->m_first.getPrev());
    temp = rel->m_first.getPrev();
    if (temp)
    {
      if (temp->m_first.getFirst())
        temp->m_first.setFirst(rel->m_first.getSame());
      else
        temp->m_first.setSame(rel->m_first.getSame());
    }
  }


  template <class T1,class T2>
    inline void RelTable<T1,T2>::removeSecond(Relation<T1,T2> *rel) {

    Relation<T1,T2> *temp;

    temp = rel->m_second.getSame();
    if (temp)
      temp->m_second.setPrev(rel->m_second.getPrev());
    temp = rel->m_second.getPrev();
    if (temp)
    {
      if (temp->m_second.getFirst())
        temp->m_second.setFirst(rel->m_second.getSame());
      else
        temp->m_second.setSame(rel->m_second.getSame());
    }
  }



}

#endif // RELTABLE_H 
