#ifndef ENUMS_TRIGGERBITS_H
#define ENUMS_TRIGGERBITS_H

namespace enums {
    typedef enum  {
        //! definition of  trigger bits

        b_ACDL =     1,  ///>  set if cover or side veto, low threshold
        b_Track=     2,  ///>  3 consecutive x-y layers hit
        b_LO_CAL=    4,  ///>  single log above low threshold
        b_HI_CAL=    8,  ///> single log above high threshold
        b_ACDH =    16,  ///>  cover or side veto, high threshold ("CNO")
        b_THROTTLE= 32,  ///> Ritz throttle

        number_of_trigger_bits = 6, ///> for size of table

        GEM_offset = 8  ///> offset to the GEM bits

    } TriggerBits;

}


#endif
