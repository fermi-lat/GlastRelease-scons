#ifndef ENUMS_GEMCONDITIONSUMMARY_H
#define ENUMS_GEMCONDITIONSUMMARY_H

namespace enums {
    
    /// Defines the 8 bits in the Condition Summary word
    /// Please see Section 1.7.2 The Condition Summary in the 
    /// "The GLT Electronics Module" available at:
    /// http://www-glast.slac.stanford.edu/IntegrationTest/ONLINE/docs/GEM.pdf
    typedef enum {
            ROI = 1,
            TKR = 2,
            CALLE = 4,
            CALHE = 8,
            CNO = 16,
            PERIODIC = 32,
            SOLICITED = 64,
            EXTERNAL = 128 
    } GemConditionSummary;


}


#endif
