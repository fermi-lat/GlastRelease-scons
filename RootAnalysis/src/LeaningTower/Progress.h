#ifndef __LEANINGTOWER_PROGRESS__
#define __LEANINGTOWER_PROGRESS__

class Progress {
 public:
    Progress() : m_fraction(0) {}
    void Go(int, int);
 private:
    int m_fraction;
};

inline void Progress::Go(int i, int n) { // some progress bar
    int fraction = (int)(100.0 * i / n);
    if ( fraction > m_fraction ) {
        m_fraction = fraction;
        if ( fraction >= 10 )
            m_fraction = m_fraction/10*10 + 9;
        std::cout << fraction << "% complete: event "<< i << std::endl;
    }
}

#endif
