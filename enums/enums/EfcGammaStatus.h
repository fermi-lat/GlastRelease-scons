#ifndef ENUMS_EFCGAMMASTATUS_H
#define ENUMS_EFCGAMMASTATUS_H

namespace enums {


// The following definitions are found in the FSW OnboardFilter header files
// in EFC/GFC_status.h. They are repeated here to avoid unnecessary dependencies
// on the FSW packages. 
/* ---------------------------------------------------------------------- *//*!

  \enum    _GFC_STAGE_S
  \brief    Gamma filter analysis stage, bit definitions, size
                                                                          *//*!
  \typedef  GFC_STAGE_S
  \brief    Typedef for enum _GFC_STAGE_S

  \par
   This maps out the size of the 2 major fields, the STAGE bit field
   and the ENERGY bit field. Apologies for the poor naming conventions
   used here
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _GFC_STAGE_S
{
    GFC_STAGE_S_ENERGY = 24,  /*!< The size of the energy field        */
    GFC_STAGE_S_STAGE  = 8    /*!< The size of the stage bit field     */ 
}
GFC_STAGE_S;
/* ---------------------------------------------------------------------- */




/* ---------------------------------------------------------------------- *//*!

  \enum    _GFC_STAGE_V
  \brief    Gamma filter analysis stage, bit definitions
                                                                          *//*!
  \typedef  GFC_STAGE_V
  \brief    Typedef for enum _GFC_STAGE_V

  \par
   This maps out the analysis phases completed by the Gamma Filter. These
   bits occupy the upper 8 bits of the stage and energy word of the
   result vector. The energy occupies the lower 24 bits.
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _GFC_STAGE_V
{
   GFC_STAGE_V_ENERGY          =  0, 
   /*!< ENERGY bit field               */

   GFC_STAGE_V_STAGE           = 24,
   /*!< STAGE bit field                */

   GFC_STAGE_V_ACD             = 24,
   /*!< ACD was analyzed               */

   GFC_STAGE_V_DIR             = 25,
   /*!< DIR was decoded                */

   GFC_STAGE_V_ATF             = 26,
   /*!< ACD/TKR veto was analyzed      */

   GFC_STAGE_V_CAL1            = 27,
   /*!< CAL was analyzed, phase 1      */

   GFC_STAGE_V_TKR             = 28,
   /*!< TKR finding was done           */

   GFC_STAGE_V_RSVD_5          = 29,
   /*!< Reserved for future use        */

   GFC_STAGE_V_RSVD_6          = 30,
   /*!< Reserved for future use        */

   GFC_STAGE_V_MPKT            = 31
   /*!< Event was a multipacket event  */
}
GFC_STAGE_V;
/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- *//*!

  \enum    _GFC_STAGE_M
  \brief    Filter stage, mask
                                                                          *//*!
  \typedef  GFC_STAGE_M
  \brief    Typedef for enum _GFC_STAGE_M
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _GFC_STAGE_M
{
   GFC_STAGE_M_ENERGY         = 0x00ffffff,
   /*!< The energy portion of the stage and energy result word            */

   GFC_STAGE_M_ACD            = (1 << GFC_STAGE_V_ACD),
   /*!< ACD was analyzed                                                  */

   GFC_STAGE_M_DIR            = (1 << GFC_STAGE_V_DIR),
   /*!< DIR was decoded                                                   */

   GFC_STAGE_M_ATF            = (1 << GFC_STAGE_V_ATF),
   /*!< ACD/TKR veto was analyzed                                         */

   GFC_STAGE_M_CAL1           = (1 << GFC_STAGE_V_CAL1),
   /*!< CAL was analyzed, phase 1                                         */

   GFC_STAGE_M_TKR            = (1 << GFC_STAGE_V_TKR),
   /*!< TKR finding was done                                              */

   GFC_STAGE_M_RSVD_5         = (1 << GFC_STAGE_V_RSVD_5),
   /*!< Reserved for future use                                           */

   GFC_STAGE_M_RSVD_6         = (1 << GFC_STAGE_V_RSVD_6),
   /*!< Reserved for future use                                           */

   GFC_STAGE_M_MPKT           = (1 << GFC_STAGE_V_MPKT),
   /*!< Event was a multi-packet event                                    */


   GFC_STAGE_M_STAGE          = GFC_STAGE_M_ACD
                                 | GFC_STAGE_M_DIR
                                 | GFC_STAGE_M_ATF
                                 | GFC_STAGE_M_CAL1
                                 | GFC_STAGE_M_TKR
                                 | GFC_STAGE_M_RSVD_5
                                 | GFC_STAGE_M_RSVD_6
                                 | GFC_STAGE_M_MPKT

}
GFC_STAGE_M;

// This for FilterStatus.h - which we really really need to get rid of...

typedef enum _EFC_GAMMA_STAGE_M
{
   EFC_GAMMA_STAGE_M_ENERGY         = 0x00ffffff
   /*!< The energy portion of the stage and energy result word            */


}
EFC_GAMMA_STAGE_M;

/* ---------------------------------------------------------------------- */




/* ---------------------------------------------------------------------- *//*!

  \enum    _GFC_STATUS_V
  \brief    Gamma filter status, bit definitions
                                                                          *//*!
  \typedef  GFC_STATUS_V
  \brief    Typedef for enum _GFC_STATUS_V
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _GFC_STATUS_V
{
   GFC_STATUS_V_GEM_THROTTLE    =  0,
   /*!< Throttle bit set                                                  */

   GFC_STATUS_V_GEM_TKR         =  1,
   /*!< Have a 3-in-a-row trigger                                         */

   GFC_STATUS_V_GEM_CALLO       =  2,
   /*!< Cal Lo Trigger                                                    */

   GFC_STATUS_V_GEM_CALHI       =  3,
   /*!< Cal Hi Trigger                                                    */

   GFC_STATUS_V_GEM_CNO         =  4,
   /*!< Cal Hi Trigger                                                    */

   GFC_STATUS_V_ACD_TOP         =  5,
   /*!< ACD top  tile struck                                              */

   GFC_STATUS_V_ACD_SIDE        =  6,
    /*!< ACD side tile struck                                             */

   GFC_STATUS_V_ACD_SIDE_FILTER =  7,
   /*!< ACD      filter tile struck                                       */

   GFC_STATUS_V_TKR_EQ_1        =  8,
   /*!< Exactly 1 track                                                   */

   GFC_STATUS_V_TKR_GE_2        =  9,
   /*!< Greater or equal 2 tracks                                         */

   GFC_STATUS_V_RSVD_10         = 10,
   /*!< Reserved for future use                                           */

   GFC_STATUS_V_RSVD_11         = 11,
   /*!< Reserved for future use                                           */

   GFC_STATUS_V_ERR_CTB         = 12,
   /*!< Contributor in error                                              */

   GFC_STATUS_V_TKR_LT_2_ELO    = 13,
   /*!< Low energy, no 2 track evidence                                   */

   GFC_STATUS_V_TKR_SKIRT       = 14,
   /*!< Event into the skirt region                                       */

   GFC_STATUS_V_TKR_EQ_0        = 15,
   /*!< No tracks                                                         */

   GFC_STATUS_V_TKR_ROW2        = 16,
   /*!< Track Row 2 match                                                 */

   GFC_STATUS_V_TKR_ROW01       = 17,
   /*!< Track Row 0 or 1 match                                            */

   GFC_STATUS_V_TKR_TOP         = 18,
   /*!< Track Top match                                                   */

   GFC_STATUS_V_ZBOTTOM         = 19,
   /*!< No tracks into CAL with energy                                    */

   GFC_STATUS_V_EL0_ETOT_HI     = 20,
   /*!< E layer 0/ETOT > .90 (90% typical)                                */

   GFC_STATUS_V_EL0_ETOT_LO     = 21,
   /*!< E layer 0/ETOT < .01 ( 1% typical)                                */

   GFC_STATUS_V_HI_ENERGY       = 22,
   /*!< Energy > pass cut, typically > 20 Gev                             */

   GFC_STATUS_V_LO_ENERGY       = 23,
   /*!< Energy <= veto cut, typically <= 2-5 Mev                          */

   GFC_STATUS_V_SIDE            = 24,
   /*!< Event has a side face veto                                        */

   GFC_STATUS_V_TOP             = 25,
   /*!< Event has a top  face veto                                        */

   GFC_STATUS_V_SPLASH_1        = 26,
   /*!< Event has a splash veto                                           */

   GFC_STATUS_V_E350_FILTER_TILE= 27,
   /*!< Event <350Mev  + filter tiles                                     */

   GFC_STATUS_V_E0_TILE         = 28,
   /*!< Event 0 energy + tile hit                                         */

   GFC_STATUS_V_SPLASH_0        = 29,
   /*!< Event has a splash veto                                           */

   GFC_STATUS_V_NOCALLO_FILTER_TILE = 30,
   /*!< No CAL LO trigger + filter tile                                   */

   GFC_STATUS_V_VETOED          = 31
   /*!< Any veto                                                          */
}
GFC_STATUS_V;
/* ---------------------------------------------------------------------- */



/* ---------------------------------------------------------------------- *//*!

  \enum    _GFC_STATUS_M
  \brief    Filter status, mask
                                                                          *//*!
  \typedef  GFC_STATUS_M
  \brief    Typedef for enum _GFC_STATUS_M
                                                                          */
/* ---------------------------------------------------------------------- */
typedef enum _GFC_STATUS_M
{
   GFC_STATUS_M_GEM_THROTTLE    = (1 << GFC_STATUS_V_GEM_THROTTLE),
   /*!< Throttle bit set                                                  */

   GFC_STATUS_M_GEM_TKR         = (1 << GFC_STATUS_V_GEM_TKR),
   /*!< Have a 3-in-a-row trigger                                         */

   GFC_STATUS_M_GEM_CALLO       = (1 << GFC_STATUS_V_GEM_CALLO),
   /*!< Cal Lo Trigger                                                    */

   GFC_STATUS_M_GEM_CALHI       = (1 << GFC_STATUS_V_GEM_CALHI),
   /*!< Cal Hi Trigger                                                    */

   GFC_STATUS_M_GEM_CNO         = (1 << GFC_STATUS_V_GEM_CNO),
   /*!< Cal Hi Trigger                                                    */

   GFC_STATUS_M_GEM             = GFC_STATUS_M_GEM_THROTTLE
                                   | GFC_STATUS_M_GEM_TKR
                                   | GFC_STATUS_M_GEM_CALLO
                                   | GFC_STATUS_M_GEM_CALHI
                                   | GFC_STATUS_M_GEM_CNO,
   /* Mask of the GEM bits                                                */

   GFC_STATUS_M_ACD_TOP         = (1 << GFC_STATUS_V_ACD_TOP),
   /*!< ACD top  tile struck                                              */

   GFC_STATUS_M_ACD_SIDE        = (1 << GFC_STATUS_V_ACD_SIDE),
    /*!< ACD side tile struck                                             */

   GFC_STATUS_M_ACD_SIDE_FILTER = (1 <<GFC_STATUS_V_ACD_SIDE_FILTER),
   /*!< ACD      filter tile struck                                       */

   GFC_STATUS_M_TKR_EQ_1        = (1 << GFC_STATUS_V_TKR_EQ_1),
   /*!< Exactly 1 track                                                   */

   GFC_STATUS_M_TKR_GE_2        = (1 << GFC_STATUS_V_TKR_GE_2),
   /*!< Greater or equal 2 tracks                                         */

   GFC_STATUS_M_RSVD_10         = (1 << GFC_STATUS_V_RSVD_10),
   /*!< Reserved for future use                                           */

   GFC_STATUS_M_RSVD_11         = (1 << GFC_STATUS_V_RSVD_11),
   /*!< Reserved for future use                                           */

   /* ------------------------------------------------------------------- */
   /*                   VETOES / PASSES                                   */
   /* ------------------------------------------------------------------- */
   GFC_STATUS_M_ERR_CTB         = (1 << GFC_STATUS_V_ERR_CTB),
   /*!< Contributor in error                                              */

   GFC_STATUS_M_TKR_LT_2_ELO    = (1 << GFC_STATUS_V_TKR_LT_2_ELO),
   /*!< Low energy, no 2 track evidence                                   */

   GFC_STATUS_M_TKR_SKIRT       = (1 << GFC_STATUS_V_TKR_SKIRT),
   /*!< Have track(s) into the skirt region                               */

   GFC_STATUS_M_TKR_EQ_0        = (1 << GFC_STATUS_V_TKR_EQ_0),
   /*!< No tracks                                                         */

   GFC_STATUS_M_TKR_ROW2        = (1 << GFC_STATUS_V_TKR_ROW2),
   /*!< Track Row 2 match                                                 */

   GFC_STATUS_M_TKR_ROW01       = (1 << GFC_STATUS_V_TKR_ROW01),
   /*!< Track Row 0 or 1 match                                            */

   GFC_STATUS_M_TKR_TOP         = (1 << GFC_STATUS_V_TKR_TOP),
   /*!< Track Top match                                                   */

   GFC_STATUS_M_ZBOTTOM         = (1 << GFC_STATUS_V_ZBOTTOM),
   /*!< No tracks into CAL with energy                                    */

   GFC_STATUS_M_EL0_ETOT_HI     = (1 <<  GFC_STATUS_V_EL0_ETOT_HI),
   /*!< E layer 0/ETOT > .90 (90% typical)                                */

   GFC_STATUS_M_EL0_ETOT_LO     = (1 <<  GFC_STATUS_V_EL0_ETOT_LO),
   /*!< E layer 0/ETOT < .01 ( 1% typical)                                */

   GFC_STATUS_M_HI_ENERGY       = (1 <<  GFC_STATUS_V_HI_ENERGY),
   /*!< E layer 0/ETOT > .90 (90% typical)                                */

   GFC_STATUS_M_LO_ENERGY       = (1 <<  GFC_STATUS_V_LO_ENERGY),
   /*!< E layer 0/ETOT < .01 ( 1% typical)                                */

   GFC_STATUS_M_SIDE            = (1 <<  GFC_STATUS_V_SIDE),
   /*!< Event has a side face veto                                        */

   GFC_STATUS_M_TOP             = (1 <<  GFC_STATUS_V_TOP),
   /*!< Event has a top  face veto                                        */

   GFC_STATUS_M_SPLASH_1        = (1 <<  GFC_STATUS_V_SPLASH_1),
   /*!< Event has a splash veto                                           */

   GFC_STATUS_M_E350_FILTER_TILE=
                                (1 <<  GFC_STATUS_V_E350_FILTER_TILE),
   /*!< Event <350Mev  + filter tiles                                     */

   GFC_STATUS_M_E0_TILE         = (1 <<  GFC_STATUS_V_E0_TILE),
   /*!< Event 0 energy + tile hit                                         */

   GFC_STATUS_M_SPLASH_0        = (1 <<  GFC_STATUS_V_SPLASH_0),
   /*!< Event has a splash veto                                           */

   GFC_STATUS_M_NOCALLO_FILTER_TILE =
                              (1 << GFC_STATUS_V_NOCALLO_FILTER_TILE),
   /*!< No CAL LO trigger + filter tile                                   */

   GFC_STATUS_M_VETOED          = (1 << GFC_STATUS_V_VETOED),
   /*!< Any veto                                                          */


   /* ------------------------------------------------------------------- */
   /* The status bit by class                                             */
   GFC_STATUS_M_INFO            = GFC_STATUS_M_GEM_THROTTLE
                                   | GFC_STATUS_M_GEM_TKR
                                   | GFC_STATUS_M_GEM_CALLO
                                   | GFC_STATUS_M_GEM_CALHI
                                   | GFC_STATUS_M_GEM_CNO
                                   | GFC_STATUS_M_ACD_TOP
                                   | GFC_STATUS_M_ACD_SIDE
                                   | GFC_STATUS_M_ACD_SIDE_FILTER
                                   | GFC_STATUS_M_TKR_EQ_1
                                   | GFC_STATUS_M_TKR_GE_2,
   /*!< The set of all informational bits                                 */


   GFC_STATUS_M_PASSES          = GFC_STATUS_M_HI_ENERGY,
   /*!< The set of all pass bits                                          */


   GFC_STATUS_M_VETOES          = GFC_STATUS_M_ERR_CTB
                                   | GFC_STATUS_M_TKR_LT_2_ELO
                                   | GFC_STATUS_M_TKR_SKIRT
                                   | GFC_STATUS_M_TKR_EQ_0
                                   | GFC_STATUS_M_TKR_ROW2
                                   | GFC_STATUS_M_TKR_ROW01
                                   | GFC_STATUS_M_TKR_TOP
                                   | GFC_STATUS_M_ZBOTTOM
                                   | GFC_STATUS_M_NOCALLO_FILTER_TILE
                                   | GFC_STATUS_M_EL0_ETOT_HI
                                   | GFC_STATUS_M_EL0_ETOT_LO
                                   | GFC_STATUS_M_LO_ENERGY  
                                   | GFC_STATUS_M_SIDE
                                   | GFC_STATUS_M_TOP
                                   | GFC_STATUS_M_SPLASH_1
                                   | GFC_STATUS_M_E350_FILTER_TILE
                                   | GFC_STATUS_M_E0_TILE
                                   | GFC_STATUS_M_SPLASH_0
   /*!< The set of all veto bits                                          */
   /* ------------------------------------------------------------------- */

}
GFC_STATUS_M;
/* ---------------------------------------------------------------------- */

} // end enums namespace

#endif
