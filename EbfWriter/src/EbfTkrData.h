#ifndef EBF_TKR_DATA_H
#define EBF_TKR_DATA_H

//
// At one time the arguments were passed by pointer,
// but now they are passed by reference. This means
// that one can no longer forward reference TkrDigiCol
//
#include "Event/Digi/TkrDigi.h"


/**
 *
 *  @class EbfTkrTowerData
 *  @brief Contains the data and methods for accessing the TKR data
 *         at the tower level.
 *
**/
class EbfTkrTowerData
{
  public:
    enum
    {
       NumTowers           = 16,
       NumCables           =  8,
       NumLayerEnds        = 72,
       NumLayerEndsPerCable=  9,
       MaxStripsPerLayerEnd= 64,
       MaxStripsPerCable   = MaxStripsPerLayerEnd*NumLayerEndsPerCable
    };


    enum eXY
    {
        X   = 0,  // Hit layer ends for X layers
        Y   = 1,  // Hit layer ends for Y layers
       XY   = 2  
    };

    enum eLOHI
    {   
       LO   = 0,
       HI   = 1,
       LOHI = 2
    };

    
    inline unsigned int mapXLO() const { return m_maps[X][LO]; }
    inline unsigned int mapXHI() const { return m_maps[X][HI]; }
    inline unsigned int mapYLO() const { return m_maps[Y][LO]; }
    inline unsigned int mapYHI() const { return m_maps[Y][HI]; }    
                                         
    
    unsigned int            m_maps[XY][LOHI]; /*!< Hit maps (18 bits per)  */
    unsigned int       m_nhits[NumLayerEnds]; /*!< Number of hits by layer */
    unsigned int        m_tots[NumLayerEnds]; /*!< TOTs by layer           */
    unsigned short int  m_data[NumLayerEnds][2*MaxStripsPerLayerEnd];
                                              /*!< The actual hit strips   */
};


/**
 *
 *  @class  EbfTkrData
 *  @brief  Contains the data and methods for assessing the TKR data at
 *          the instrument level.
 *
**/ 
class EbfTkrData
{
  public:
    void                 initialize ();
    void                       fill (const Event::TkrDigiCol &tkr);
    unsigned int            *format (unsigned int *dst)              const;
    unsigned int            *format (unsigned int *dst, int towerId) const;
    void print                      ();

    /**
     *
     *  @fn     unsigned int threeInARow () const
     *  @brief  Returns a bit list of the towers with a 3-in-a-row TKR
     *          trigger primitive set.
     *
     *  @return A bit list of the towers with a 3-in-a-row TKR trigger
     *          primitive set. The convention is Tower 0 = LSB.
     *
    **/
    inline unsigned int threeInARow () const { return m_trigger;  }


    /**
     *
     *  @fn     const EbfTkrTowerData *tower(int itower) const
     *  @brief  Returns a pointer to the tower level description of the
     *          TKR data.
     *
     *  @param itower  The tower to return the pointer for (0-15).
     *  @return        A pointer to the tower level description of the
     *                 TKR data.
     *
    **/
    inline const EbfTkrTowerData *tower(int itower) const
    {
        return m_towers + itower;
    }
    
            
    
  private:
    unsigned int                              m_trigger;
    EbfTkrTowerData m_towers[EbfTkrTowerData::NumTowers];
};




#endif




