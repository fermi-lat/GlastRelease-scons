#ifndef EBF_CAL_DATA_H
#define EBF_CAL_DATA_H


//
// At one time the arguments were passed by pointer,
// but now they are passed by reference. This means
// that one can no longer forward reference CalDigiCol
//
#include "Event/Digi/CalDigi.h"

#include "CalXtalResponse/ICalTrigTool.h"

class IGlastDetSvc;


/**
 *   @class  EbfCalLayerData
 *   @brief  Describes the CAL data for an event at the layer level
**/
class EbfCalLayerData
{
  public:
    enum      { NumLogs = 12 };
    
    unsigned short int       m_lo;  /*!< Bit mask of lo discriminator hits */
    unsigned short int       m_hi;  /*!< Bit mask of hi discriminator hits */
    unsigned short int    m_nhits;  /*!< Number of hits                    */
    unsigned short int      m_msk;  /*!< Bit mask of crystals              */
    unsigned int m_xtals[NumLogs][4];/*!< The ADCs allow for 4 range readout*/
};




/**
 *   @class  EbfCalTowerData
 *   @brief  Describes the CAL data for an event at the tower level
**/
class EbfCalTowerData
{
   public:

     /**
      *
      *  @fn     unsigned short loTrigger () const
      *  @brief  Returns a bit mask of the layers of this tower with the
      *          CAL LO trigger primitive active
      *  @return A bit mask of the layers of this tower with the CAL LO
      *          LO trigger primitive active
      *
      *   Returns the mask of the layers of this tower with the CAL LO
      *   trigger primitive active. Layer 0 corresponds to the LSB.
      *
     **/
     inline unsigned short loTrigger () const { return m_lo; };


     /**
      *
      *  @fn     unsigned short hiTrigger () const
      *  @brief  Returns a bit mask of the layers of this tower with the
      *          CAL HI trigger primitive active
      *  @return A bit mask of the layers of this tower with the CAL HI
      *          HI trigger primitive active
      *
      *   Returns the mask of the layers of this tower with the CAL HI
      *   trigger primitive active. Layer 0 corresponds to the LSB.
      *
     **/
     inline unsigned short hiTrigger () const { return m_hi; };
 
     enum               { NumLayers = 8 };
     unsigned short int             m_lo; /*!< Mask of layers with CAL LO */
     unsigned short int             m_hi; /*!< Mask of layers with CAL HI */
     unsigned short int            m_msk; /*!< Mask of layers with data   */
     EbfCalLayerData m_layers[NumLayers]; /*!< The data for each layer    */
};




/**
 *   @class  EbfCalData
 *   @brief  Describes the CAL data for an event at the instrument level
**/
class EbfCalData
{
  public:
    enum { NumTowers = 16, NumLayers = EbfCalTowerData::NumLayers };
    
    
    void      initialize ();

    /// fill internal data arrays from CalDigi information, prepare
    /// for format() to EBF output.
    StatusCode            fill (const Event::CalDigiCol &calDigiCol,
                                ICalTrigTool &calTrigTool);

    /// \brief Fill internal data arrays with special bit pattern fro
    /// FET testing.
    /// 
    /// All events should to trigger, so default is to set
    /// all Cal trigger bits true
    void            fillEncode (const int encodeFlag, 
                                const int event,
                                const bool leTriggerBit=true,
                                const bool heTriggerBit=true
                                );
    /// create EBF file format output from internal flat arrays for
    /// whole cal
    unsigned int *format (unsigned int  *dst)              const;       
    /// create EBF file format output from internal flat arrays for
    /// single tower
    unsigned int *format (unsigned int  *dst, int towerId) const;

    void          print  ();

    /**
     *
     *  @fn     unsigned short loTrigger () const
     *  @brief  Returns a bit mask of the towers with CAL LO trigger
     *          primitive active
     *  @return A bit mask of the towers with the CAL LO trigger
     *          primitive active
     *
     *   Returns the mask of the towers with the CAL LO trigger primitive
     *   active. Tower 0 corresponds to the LSB.
     *
    **/
    inline unsigned short loTrigger () const { return m_lo; };


    /**
     *
     *  @fn     unsigned short hiTrigger () const
     *  @brief  Returns a bit mask of the towers with CAL HI trigger
     *          primitive active
     *  @return A bit mask of the towers with the CAL HI trigger
     *          primitive active
     *
     *   Returns the mask of the towers with the CAL HI trigger primitive
     *   active. Tower 0 corresponds to the LSB.
     *
     */
    inline unsigned short hiTrigger () const { return m_hi; };
    

    /**
     *
     *  @fn     const EbfCalTowerData *tower (int itower) const
     *  @brief  Returns a pointer to the CAL tower data for this event.
     *  @return A pointer to the CAL tower data for this event.
     *
    **/
    inline const EbfCalTowerData *tower (int itower) const
    {
        return m_towers + itower;
    }
            
    inline bool FourRangeReadOut () const
    {
        return m_range4;
    }    
    inline bool getCalStrobe () const
    {
        return m_calStrobe;
    }        
    
  private:
    unsigned short int             m_lo; /*!< Mask of towers with CAL LO  */
    unsigned short int             m_hi; /*!< Mask of towers with CAL HI  */
    unsigned short int            m_msk; /*!< Mask of towers with data    */
    EbfCalTowerData m_towers[NumTowers]; /*!< The data for all the towers */ 
    bool                       m_range4; /*!< TRUE=4 range readout enabled*/
    bool                    m_calStrobe; /*!< calStrobe set */

};



#endif




