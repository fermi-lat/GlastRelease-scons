#ifndef ENUMS_XFCMIPSTATUS_H
#define ENUMS_XFCMIPSTATUS_H

namespace enums {


// The following definitions are found in the FSW OnboardFilter header files
// in XFC/HFC_status.h. They are repeated here to avoid unnecessary dependencies
// on the FSW packages. 

/* ---------------------------------------------------------------------- *//*!

  \enum    _MFC_STATUS_V
  \brief    Heavy Ion filter status, bit definitions
                                                                          *//*!
  \typedef  MFC_STATUS_V
  \brief    Typedef for enum _MFC_STATUS_V
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _MFC_STATUS_V
{
    MFC_STATUS_V_STAGE_GEM       =  0,
    /*!< GEM data analyzed          */
    
    MFC_STATUS_V_STAGE_ACD       =  1,
    /*!< ACD checked                */

    MFC_STATUS_V_STAGE_DIR       =  2,
    /*!< DIR formed                 */

    MFC_STATUS_V_STAGE_CAL       =  3,
    /*!< CAL unpacked               */

    MFC_STATUS_V_STAGE_XCAL      =  4,
    /*!< CAL checked                */

    MFC_STATUS_V_MULTI_PKT       =  5,
    /*!< Multi-packet event         */

    MFC_STATUS_V_ERR_CAL         =  6,
    /*!< CAL data in error          */
    
    MFC_STATUS_V_ERR_CTB         =  7,
    /*!< Contributor in error       */
    
    MFC_STATUS_V_ERR_DIR         =  8,
    /*!< Directory in error         */
    
    MFC_STATUS_V_RSVD_09         =  9,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_10         = 10,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_11         = 11,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_12         = 12,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_13         = 13,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_14         = 14,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_15         = 15,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_16         = 16,
    /*!< Reserved for future use    */

    MFC_STATUS_V_RSVD_17         = 17,
    /*!< Reserved for future use    */

    MFC_STATUS_V_RSVD_18         = 18,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_19         = 19,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_20         = 20,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_21         = 21,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_22         = 22,
    /*!< Reserved for future use    */

    MFC_STATUS_V_RSVD_23         = 23,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_24         = 24,
    /*!< Reserved for future use    */
    
    MFC_STATUS_V_RSVD_25         = 25,
    /*!< Reserved for future use    */

    MFC_STATUS_V_NO_TKR_ADJ      = 26,
    /*!< No 2 adjacent towers with a 3-in-a-row trigger */

    MFC_STATUS_V_LYR_COUNTS      = 27,
    /*!< Insufficient number of logs in the 1 MIP range */

    MFC_STATUS_V_NO_ACD_TKR      = 28,
    /*!< No ACD/TKR coincidence                         */
    
    MFC_STATUS_V_GEM_NOTKR       = 29,
    /*!< No GEM 3-in-a-row trigger */
    
    MFC_STATUS_V_GEM_CNO         = 30,
    /*!< CNO bit in the GEM is set */
    
    MFC_STATUS_V_VETOED          = 31
    /*!< Any veto                   */
}
MFC_STATUS_V;
/* ---------------------------------------------------------------------- */




/* ---------------------------------------------------------------------- *//*!

  \enum    _MFC_STATUS_M
  \brief    Filter status, mask
                                                                          *//*!
  \typedef  MFC_STATUS_M
  \brief    Typedef for enum _MFC_STATUS_M
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _MFC_STATUS_M
{
    /* The GEM bits map directly onto the GEM condition summary word      */
    MFC_STATUS_M_STAGE_GEM       = 1 << MFC_STATUS_V_STAGE_GEM,
    /*!< GEM data analyzed                                                */

    MFC_STATUS_M_STAGE_ACD       = 1 << MFC_STATUS_V_STAGE_ACD,
    /*!< ACD checked                                                      */

    MFC_STATUS_M_STAGE_DIR       = 1 << MFC_STATUS_V_STAGE_DIR,
    /*!< DIR formed                                                       */

    MFC_STATUS_M_STAGE_CAL       = 1 << MFC_STATUS_V_STAGE_CAL,
    /*!< CAL unpacked                                                     */

    MFC_STATUS_M_STAGE_XCAL      = 1 << MFC_STATUS_V_STAGE_XCAL,
    /*!< CAL checked                                                      */

    MFC_STATUS_M_MULTI_PKT       = 1 << MFC_STATUS_V_MULTI_PKT,
    /*!< Multi-packet event                                               */

    MFC_STATUS_M_ERR_CAL         = 1 << MFC_STATUS_V_ERR_CAL,
    /*!< Reserved for future use                                          */

    MFC_STATUS_M_ERR_CTB         = 1 << MFC_STATUS_V_ERR_CTB,
    /*!< Contributor in error                                             */

    MFC_STATUS_M_ERR_DIR         = 1 << MFC_STATUS_V_ERR_DIR,
    /*!< Reserved for future use                                          */
        
    MFC_STATUS_M_RSVD_09         = 1 << MFC_STATUS_V_RSVD_09,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_10         = 1 << MFC_STATUS_V_RSVD_10,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_11         = 1 << MFC_STATUS_V_RSVD_11,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_12         = 1 << MFC_STATUS_V_RSVD_12,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_13         = 1 << MFC_STATUS_V_RSVD_13,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_14         = 1 << MFC_STATUS_V_RSVD_14,
    /*!< Reserved for future use                                          */

    MFC_STATUS_M_RSVD_15         = 1 << MFC_STATUS_V_RSVD_15,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_16         = 1 << MFC_STATUS_V_RSVD_16,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_17         = 1 << MFC_STATUS_V_RSVD_17,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_18         = 1 << MFC_STATUS_V_RSVD_18,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_19         = 1 << MFC_STATUS_V_RSVD_19,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_20         = 1 << MFC_STATUS_V_RSVD_20,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_21         = 1 << MFC_STATUS_V_RSVD_21,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_22         = 1 << MFC_STATUS_V_RSVD_22,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_23         = 1 << MFC_STATUS_V_RSVD_23,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_24         = 1 << MFC_STATUS_V_RSVD_24,
    /*!< Reserved for future use                                          */
    
    MFC_STATUS_M_RSVD_25         = 1 << MFC_STATUS_V_RSVD_25,
    /*!< Reserved for future use                                          */

    MFC_STATUS_M_NO_TKR_ADJ      = 1 << MFC_STATUS_V_NO_TKR_ADJ,
    /*!< No 2 adjacent towers with a 3-in-a-row trigger */
    
    MFC_STATUS_M_LYR_COUNTS      = 1 << MFC_STATUS_V_LYR_COUNTS,
    /*!< Insufficient number of logs over threshold in first three layers */

    MFC_STATUS_M_NO_ACD_TKR      = 1 << MFC_STATUS_V_NO_ACD_TKR,
    /*!< No ACD/TKR coincidence                                           */
    
    MFC_STATUS_M_GEM_NOTKR       = 1 << MFC_STATUS_V_GEM_NOTKR,
    /*!< No GEM 3-in-a-row                                                */

    MFC_STATUS_M_GEM_CNO         = 1 << MFC_STATUS_V_GEM_CNO,
    /*!< CNO bit in the GEM is not set                                    */
    
    MFC_STATUS_M_VETOED          = 1 << MFC_STATUS_V_VETOED,
    /*!< Any veto                                                         */

    MFC_STATUS_M_VETO_DEF        =  MFC_STATUS_M_VETOED
                                 |  MFC_STATUS_M_GEM_CNO
                                 |  MFC_STATUS_M_NO_ACD_TKR
                                 |  MFC_STATUS_M_GEM_NOTKR
                                 |  MFC_STATUS_M_LYR_COUNTS
                                 |  MFC_STATUS_M_NO_TKR_ADJ
    /*!< Convenience symbol for all possible vetoes                       */

}
MFC_STATUS_M;
/* ---------------------------------------------------------------------- */

} // end enums namespace

#endif
