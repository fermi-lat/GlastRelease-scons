#ifndef ENUMS_XFCDFCSTATUS_H
#define ENUMS_XFCDFCSTATUS_H

namespace enums {


// The following definitions are found in the FSW OnboardFilter header files
// in XFC/HFC_status.h. They are repeated here to avoid unnecessary dependencies
// on the FSW packages. 

/* ---------------------------------------------------------------------- *//*!

  \enum    _DFC_STATUS_V
  \brief    Heavy Ion filter status, bit definitions
                                                                          *//*!
  \typedef  DFC_STATUS_V
  \brief    Typedef for enum _DFC_STATUS_V
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _DFC_STATUS_V
{
    DFC_STATUS_V_STAGE_GEM       =  0,
    /*!< GEM data analyzed          */

    DFC_STATUS_V_RSVD_01         =  1,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_02         =  2,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_03         =  3,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_04         =  4,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_05         =  5,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_06         =  6,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_07         =  7,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_08         =  8,
    /*!< Directory in error         */

    DFC_STATUS_V_RSVD_09         =  9,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_10         = 10,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_11         = 11,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_12         = 12,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_13         = 13,
    /*!< Reserved for future use    */

    DFC_STATUS_V_RSVD_14         = 14,
    /*!< Reserved for future use    */

    DFC_STATUS_V_GEM_CLASSES     = 15,
    /*!< Start of the 16 GEM classes */

    DFC_STATUS_V_VETOED          = 31
    /*!< Any veto                   */
}
DFC_STATUS_V;
/* ---------------------------------------------------------------------- */




/* ---------------------------------------------------------------------- *//*!

  \enum    _DFC_STATUS_M
  \brief    Filter status, mask
                                                                          *//*!
  \typedef  DFC_STATUS_M
  \brief    Typedef for enum _DFC_STATUS_M
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _DFC_STATUS_M
{
    DFC_STATUS_M_STAGE_GEM       = 1 << DFC_STATUS_V_STAGE_GEM,
    /*!< GEM data analyzed                                                */

    DFC_STATUS_M_RSVD_01         = 1 << DFC_STATUS_V_RSVD_01,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_02         = 1 << DFC_STATUS_V_RSVD_02,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_03         = 1 << DFC_STATUS_V_RSVD_03,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_04         = 1 << DFC_STATUS_V_RSVD_04,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_05         = 1 << DFC_STATUS_V_RSVD_05,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_06         = 1 << DFC_STATUS_V_RSVD_06,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_07         = 1 << DFC_STATUS_V_RSVD_07,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_08         = 1 << DFC_STATUS_V_RSVD_08,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_09         = 1 << DFC_STATUS_V_RSVD_09,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_10         = 1 << DFC_STATUS_V_RSVD_10,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_11         = 1 << DFC_STATUS_V_RSVD_11,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_12         = 1 << DFC_STATUS_V_RSVD_12,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_13         = 1 << DFC_STATUS_V_RSVD_13,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_RSVD_14         = 1 << DFC_STATUS_V_RSVD_14,
    /*!< Reserved for future use                                          */

    DFC_STATUS_M_GEM_CLASSES     = 0xffff << DFC_STATUS_V_GEM_CLASSES,
    /*!< The 16 GEM classes                                               */

    DFC_STATUS_M_VETOED          = 1 << DFC_STATUS_V_VETOED,
    /*!< Any veto                                                         */

    DFC_STATUS_M_VETO_DEF        =  DFC_STATUS_M_GEM_CLASSES
    /*!< Convenience symbol for all possible vetoes                       */

}
DFC_STATUS_M;
/* ---------------------------------------------------------------------- */

} // end enums namespace

#endif
