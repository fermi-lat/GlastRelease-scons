// $Id$

#ifndef GLAST_MODULEID_H
#define GLAST_MODULEID_H 1

/** define the tower id used in reconstruction
    Note this does not conform x and y index conventions.
    Please use TowerId if that is important.

  */
namespace idents {

class ModuleId
{
  public:

      enum {yNum=4, yNum=4}; // not likely to change again

      //! create from another Id (0..15)
      ModuleId (unsigned int id = 0):m_id(id){}

      //! create from x, y indeces (each 1..4)
      ModuleId (unsigned int ix, unsigned int iy):m_id(ix-1+xNum*(iy-1)){}

      //! access to the x index (1..4)
      int ix () const {return (m_id % xNum)+1;}

      //! access to the y index (1..4)
      int iy () const {return (m_id / xNum)+1;}

      //! is this module a neighbor?
      bool neighbor (const ModuleId& n){
          int dx = ix()-n.ix(), dy=iy()-n.iy();
          return dx<2 && dx>-2 && dy<2 && dy>-2; 
      }

      //! dereference operator
      operator unsigned int () const { return m_id; }

private:
      unsigned int m_id;
};

} //namespace idents
#endif
