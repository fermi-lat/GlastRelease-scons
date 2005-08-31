#ifndef EBF_CAL_DATA_H
#define EBF_CAL_DATA_H


//
// At one time the arguments were passed by pointer,
// but now they are passed by reference. This means
// that one can no longer forward reference CalDigiCol
//
#include "GaudiKernel/MsgStream.h"
#include "Event/Digi/CalDigi.h"

class IGlastDetSvc;



/**
 *  @class  EbfCalContants
 *  @brief  Contains the constants and methods to convert an CAL
 *          ADC value to an energy in MEV. Also contains the
 *          discriminator thresholds for CAL HI and CAL LO trigger
 *          primitives.
**/
class EbfCalConstants
{
  public:
     
    StatusCode initialize (IGlastDetSvc *detSvc, MsgStream &log);

    /**
     *
     *  @fn      double convertToMev (int range, unsigned int adc) const
     *  @brief   Converts an ADC value to the an energy in MeEV
     *
     *  @param   range The gain range to use in performing the conversion
     *  @param     adc The ADC to convert
     *  @return        The equivalent energy in MEV
     *
    **/
    inline double convertToMev (int range, unsigned int adc) const
    {
        return m_adcToMev[range] * (adc - m_pedestal);
    };

    
    double m_pedestal;       /*!< ADC pedestal value                 */
    double m_maxAdc;         /*!< Cross-over point for range switch  */
    double m_maxEnergy[4];   /*!< Maximum enery in each gain range   */
    double m_loThreshold;    /*!< CAL LO trigger primitive threshold */
    double m_hiThreshold;    /*!< CAL HI trigger primitive threshold */
    double m_adcToMev[4];    /*!< Effective gain for each range      */
};





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
    void            fill (const Event::CalDigiCol &calDigiCol,
                          const EbfCalConstants   &constants);
    void            fillEncode (int encodeFlag, const EbfCalConstants   &constants, int event);
    void            parseInput (unsigned int *contrib, unsigned int tower, unsigned int lcbWords, const EbfCalConstants *calCon);
    
    unsigned int *format (unsigned int  *dst)              const;       
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
    inline double getTotalEnergy() const
    {
        return m_TotalEnergy;
    }
    
  private:
    unsigned short int             m_lo; /*!< Mask of towers with CAL LO  */
    unsigned short int             m_hi; /*!< Mask of towers with CAL HI  */
    unsigned short int            m_msk; /*!< Mask of towers with data    */
    EbfCalTowerData m_towers[NumTowers]; /*!< The data for all the towers */ 
    bool                       m_range4; /*!< TRUE=4 range readout enabled*/
    bool                    m_calStrobe; /*!< calStrobe set */
    double                m_TotalEnergy; /*!< Total Energy deposited in cal */                            
};



#endif




