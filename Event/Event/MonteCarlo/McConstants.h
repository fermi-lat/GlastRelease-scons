// $Header :$
#ifndef GlastEvent_McConstants_H
#define GlastEvent_McConstants_H 1


namespace GlastEvent {
    namespace McConstants {
        /// Status flags
        // McParticle status
        const unsigned long  PRIMARY            = 1 << 0;
        const unsigned long  CALOSHOWER         = 1 << 1;
        // Hit flags
        const unsigned long  ORIGIN_PRIMARY     = 1 << 0;
        const unsigned long  ORIGIN_CALOSHOWER  = 1 << 1;
        const unsigned long  NEED_DIGI          = 1 << 8;
    }
}

#endif // GlastEvent_McConstants_H
