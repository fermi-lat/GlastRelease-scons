#ifndef EBF_GEM_DATA_H
#define EBF_GEM_DATA_H


//#include "Event/TopLevel/Event.h"
#include "GaudiKernel/IInterface.h"
#include "LdfEvent/Gem.h"

class EbfAcdData;
class EbfTkrData;
class EbfCalData;
class EbfGemCounters;


/**
 *
 *  @class EbfGemData
 *  @brief Contains the data necessary to form the contribution of the
 *         GEM to the event.
 *
**/
class EbfGemData
{
  public:

    /**
     *
     *  @enum  TileFaces
     *  @brief Names the faces of the ACD by specifying the two axis
     *         which the plane lies in. For example, XY = Top plane.
     * 
    **/
    enum TileFaces
    {
      XZ = 0,  /*!< Tiles measuring X and Z dimensions                  */
      YZ = 1,  /*!< Tiles measuring Y and Z dimensions                  */
      XY = 2,  /*!< Tiles measuring X and Y dimensions, the TOP tiles   */
      RU = 3,  /*!< The 8 ribbons + the unused channels                 */
      NumMasks /*!< The number of 32 bit words needed to hold the masks */
    };

    /**
     *
     *  @enum  TReqTypes
     *  @brief Names the different trigger request types.
     *
     *   These determine the meaning of the bits in the trigger request
     *   list.
     *
    **/
    enum TReqTypes
    {
                           /* Physics triggers                          */
      TKR_NOT_VETOED =  0, /* Tracker 3 in a row, not vetoed by ACD     */
      TKR            =  1, /* Tracker 3 in a row                        */
      CAL_LO         =  2, /* CAL LO in any tower                       */
      CAL_HI         =  3, /* CAL HI in any tower                       */
      CNO            =  4, /* Any ACD CNO discriminator above threshold */      
      PERIODIC       =  5, /* Periodic (diagnostic trigger)             */
      SOLICITED      =  6, /* Trigger solicited by the CPU              */
      EXTERNAL       =  7  /* External Trigger                          */
    };
    
    
      
    inline void         initialize  ()                 { return; }
    
    void                fill        (LdfEvent::Gem*        gemTds,
                                     unsigned int       mc_number,
                                     double               mc_time,
                                     EbfGemCounters *lat_counters,
                                     const EbfAcdData        *acd,
                                     const EbfTkrData        *tkr,
                                     const EbfCalData       *cal);
    
    unsigned int       *format      (unsigned int *dst) const;
    
    void                print       () const;

    /**
     *
     *  @fn     unsigned int calHiLo () const
     *  @brief  Returns the bit mask of towers with the CAL HI or CAL LO
     *          trigger primitives active.
     *  @return Two concatenated bit lists, one representing the list of
     *          of towers with its CAL HI trigger primitive active and the
     *          other representing the list of towers with its CAL LO
     *          trigger primitive active.
     *
     *   This routine returns a bit mask of towers with the active CAL HI
     *   or CAL LO trigger primitives. The upper 16 bits contain the list
     *   of towers with CAL HI set, while the lower 16 bits contain the
     *   list of towers with CAL LO set. In both lists, TOWER LO corresponds
     *   to the LSB.
     *
    **/
    inline unsigned int calHiLo () const { return m_calHiLo; }


    /**
     *
     *  @fn     unsigned int thrTkr () const
     *  @brief  Returns the bit mask of towers that are shadowed by the
     *          ACD struck tiles (the so-called throttle mask) and a bit
     *          mask of the towers with a 3-in-a-row TKR trigger.
     *  @return Two bit lists, one representing the throttled towers and
     *          the other representing the list of towers with a 3-in-a-row
     *          TKR trigger.
     *
     *   This routine returns two concatenated 16 bit lists. The upper 16
     *   bits contain the list of towers shadowed by the struck ACD tiles
     *   in this event. The lower 16 bits contain the list of towers with
     *   a 3-in-a-row TKR trigger. In both lists TOWER 0 corresponds to the
     *   LSB.
     *
    **/
    inline unsigned int thrTkr      () const { return m_thrTkr;             }


    /**
     *
     *  @fn     unsigned int isTriggered () const
     *  @brief  Returns a value indicating whether this event resulted
     *          in the GEM declaring a trigger (TACK).
     *
     *  @retval 0  if this event did not result in a trigger
     *  @retval !0 if this event did result in a trigger
     *
     *   This routine returns a non-zero value if this event resulted in
     *   the GEM declaring a trigger. Note that this is not a boolean value,
     *   it is simply 0 if no trigger was declared and non-zero if a trigger
     *   was declared.
    **/
    inline unsigned int isTriggered () const { return (m_cnoReqvec & 0xff0000)>>16; }
    
    
  private:
    unsigned int          m_thrTkr; /*!< vetoed TKR and TKR                  */
    unsigned int         m_calHiLo; /*!< Cal trigger bits                    */
    unsigned int       m_cnoReqvec; /*!< TCNO bit mask and
                                         Trigger Request Vector              */
    unsigned int m_tiles[NumMasks]; /*!< Top, X+/X-, Y+/Y-                   */
    unsigned int        m_livetime; /*!< Livetime counter in sys clock ticks */
    unsigned int            m_busy; /*!< Deadtime counter in # of events     */
    unsigned int       m_prescaled; /*!< # of events pitched by prescale     */
    unsigned int            m_sent; /*!< # of events triggered               */
    unsigned int         m_evttime; /*!< Time of event in sys clock ticks    */
    unsigned int     m_lastEvtTime; /*!< Time of Previous Event */
    unsigned int       m_deltaTime; /*!< Delta Time between events */
    unsigned int         m_ppstime; /*!< System clock at last PPS plus (in
                                         the upper N bits PPS reference #    */
    unsigned int          m_number; /*!< The MC event number                 */
    double                  m_time; /*!< The MC event time                   */
};




#endif




