#ifndef ENUMS_XFCHFCSTATUS_H
#define ENUMS_XFCHFCSTATUS_H

namespace enums {


// The following definitions are found in the FSW OnboardFilter header files
// in XFC/HFC_status.h. They are repeated here to avoid unnecessary dependencies
// on the FSW packages. 

/* ---------------------------------------------------------------------- *//*!

  \enum    _HFC_STATUS_V
  \brief    Heavy Ion filter status, bit definitions
                                                                          *//*!
  \typedef  HFC_STATUS_V
  \brief    Typedef for enum _HFC_STATUS_V
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _HFC_STATUS_V
{
    HFC_STATUS_V_STAGE_GEM       =  0,
    /*!< GEM data analyzed          */
    
    HFC_STATUS_V_STAGE_DIR       =  1,
    /*!< DIR formed                 */
    
    HFC_STATUS_V_STAGE_CAL       =  2,
    /*!< CAL unpacked               */
    
    HFC_STATUS_V_STAGE_CAL_ECHK  =  3,
    /*!< CAL energy checked         */

    HFC_STATUS_V_STAGE_CAL_LCHK  =  4,
    /*!< CAL layer energy checked   */

    HFC_STATUS_V_MULTI_PKT       =  5,
    /*!< Multi-packet event         */

    HFC_STATUS_V_ERR_CAL         =  6,
    /*!< CAL data in error          */
    
    HFC_STATUS_V_ERR_CTB         =  7,
    /*!< Contributor in error       */
    
    HFC_STATUS_V_ERR_DIR         =  8,
    /*!< Directory in error         */
    
    HFC_STATUS_V_RSVD_09         =  9,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_10         = 10,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_11         = 11,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_12         = 12,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_13         = 13,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_14         = 14,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_15         = 15,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_16         = 16,
    /*!< Reserved for future use    */

    HFC_STATUS_V_RSVD_17         = 17,
    /*!< Reserved for future use    */

    HFC_STATUS_V_RSVD_18         = 18,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_19         = 19,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_20         = 20,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_21         = 21,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_22         = 22,
    /*!< Reserved for future use    */

    HFC_STATUS_V_RSVD_23         = 23,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_24         = 24,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_RSVD_25         = 25,
    /*!< Reserved for future use    */
    
    HFC_STATUS_V_LYR_ENERGY      = 26,
    /*!< Event CAL layer energies do not match closely enough */
    
    HFC_STATUS_V_LYR_COUNTS      = 27,
    /*!< Insufficient number of logs over threshold in first three layers */
    
    HFC_STATUS_V_GEM_NOTKR       = 28,
    /*!< No GEM 3-in-a-row trigger */
    
    HFC_STATUS_V_GEM_NOCALLO     = 29,
    /*!< CAL LO bit in the GEM is not set  */
    
    HFC_STATUS_V_GEM_NOCNO       = 30,
    /*!< CNO bit in the GEM is not set */
    
    HFC_STATUS_V_VETOED          = 31
    /*!< Any veto                   */
}
HFC_STATUS_V;
/* ---------------------------------------------------------------------- */




/* ---------------------------------------------------------------------- *//*!

  \enum    _HFC_STATUS_M
  \brief    Filter status, mask
                                                                          *//*!
  \typedef  HFC_STATUS_M
  \brief    Typedef for enum _HFC_STATUS_M
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _HFC_STATUS_M
{
    /* The GEM bits map directly onto the GEM condition summary word      */
    HFC_STATUS_M_STAGE_GEM       = 1 << HFC_STATUS_V_STAGE_GEM,
    /*!< GEM data analyzed                                                */

    HFC_STATUS_M_STAGE_DIR       = 1 << HFC_STATUS_V_STAGE_DIR,
    /*!< DIR formed                                                       */
    
    HFC_STATUS_M_STAGE_CAL       = 1 << HFC_STATUS_V_STAGE_CAL,
    /*!< CAL unpacked                                                     */

    HFC_STATUS_M_STAGE_CAL_ECHK  = 1 << HFC_STATUS_V_STAGE_CAL_ECHK,
    /*!< CAL energy checked                                               */

    HFC_STATUS_M_STAGE_CAL_LCHK  = 1 << HFC_STATUS_V_STAGE_CAL_LCHK,
    /*!< CAL layer checked                                                */
    
    HFC_STATUS_M_MULTI_PKT       = 1 << HFC_STATUS_V_MULTI_PKT,
    /*!< Multi-packet event                                               */

    HFC_STATUS_M_ERR_CAL         = 1 << HFC_STATUS_V_ERR_CAL,
    /*!< Reserved for future use                                          */

    HFC_STATUS_M_ERR_CTB         = 1 << HFC_STATUS_V_ERR_CTB,
    /*!< Contributor in error                                             */

    HFC_STATUS_M_ERR_DIR         = 1 << HFC_STATUS_V_ERR_DIR,
    /*!< Reserved for future use                                          */
        
    HFC_STATUS_M_RSVD_09         = 1 << HFC_STATUS_V_RSVD_09,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_10         = 1 << HFC_STATUS_V_RSVD_10,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_11         = 1 << HFC_STATUS_V_RSVD_11,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_12         = 1 << HFC_STATUS_V_RSVD_12,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_13         = 1 << HFC_STATUS_V_RSVD_13,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_14         = 1 << HFC_STATUS_V_RSVD_14,
    /*!< Reserved for future use                                          */

    HFC_STATUS_M_RSVD_15         = 1 << HFC_STATUS_V_RSVD_15,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_16         = 1 << HFC_STATUS_V_RSVD_16,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_17         = 1 << HFC_STATUS_V_RSVD_17,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_18         = 1 << HFC_STATUS_V_RSVD_18,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_19         = 1 << HFC_STATUS_V_RSVD_19,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_20         = 1 << HFC_STATUS_V_RSVD_20,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_21         = 1 << HFC_STATUS_V_RSVD_21,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_22         = 1 << HFC_STATUS_V_RSVD_22,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_23         = 1 << HFC_STATUS_V_RSVD_23,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_24         = 1 << HFC_STATUS_V_RSVD_24,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_RSVD_25         = 1 << HFC_STATUS_V_RSVD_25,
    /*!< Reserved for future use                                          */
    
    HFC_STATUS_M_LYR_ENERGY      = 1 << HFC_STATUS_V_LYR_ENERGY,
    /*!< Event CAL layer energies do not match closely enough             */
    
    HFC_STATUS_M_LYR_COUNTS      = 1 << HFC_STATUS_V_LYR_COUNTS,
    /*!< Insufficient number of logs over threshold in first three layers */
    
    HFC_STATUS_M_GEM_NOTKR       = 1 << HFC_STATUS_V_GEM_NOTKR,
    /*!< No GEM 3-in-a-row                                                */

    HFC_STATUS_M_GEM_NOCALLO     = 1 << HFC_STATUS_V_GEM_NOCALLO,
    /*!< CAL LO bit in the GEM is not set                                 */

    HFC_STATUS_M_GEM_NOCNO       = 1 << HFC_STATUS_V_GEM_NOCNO,
    /*!< CNO bit in the GEM is not set                                    */
    
    HFC_STATUS_M_VETOED          = 1 << HFC_STATUS_V_VETOED,
    /*!< Any veto                                                         */

    HFC_STATUS_M_VETO_DEF        =  HFC_STATUS_M_VETOED
                                 |  HFC_STATUS_M_GEM_NOCNO
                                 |  HFC_STATUS_M_GEM_NOCALLO
                                 |  HFC_STATUS_M_GEM_NOTKR
                                 |  HFC_STATUS_M_LYR_COUNTS
                                 |  HFC_STATUS_M_LYR_ENERGY
    /*!< Convenience symbol for all possible vetoes                       */

}
HFC_STATUS_M;
/* ---------------------------------------------------------------------- */

} // end enums namespace

#endif
