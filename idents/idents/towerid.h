
#ifndef GLAST_TOWERID_H
#define GLAST_TOWERID_H 1

/** define the tower id following the Ritz specs

  */
namespace idents {

class TowerId
{
  public:

      /*! Number of Y towers is used for the tower arithmetic
          Luckily, this works for the BTEM/BFEM, because the tower
          id is always zero anyway. 
          This won't necessarily work for the 2x2 configuration,
          unless we agree to number the towers 0, 1, 4, 5.
      */
      enum {yNum=4}; // not likely to change again

      //! create from another Id (0..15)
      TowerId (unsigned int id = 0):m_id(id){}

      //! create from x, y indices (each 0..3)
      TowerId (unsigned int ix, unsigned int iy):m_id(ix+yNum*iy){}

      //! access to the x index (0..3)
      int ix () const {return (m_id % yNum);}

      //! access to the y index (0..3)
      int iy () const {return (m_id / yNum);}

      //! is this module a neighbor?
      bool neighbor (const TowerId& n){
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
