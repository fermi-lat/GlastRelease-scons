#ifndef __LEANINGTOWER_PROGRESS__
#define __LEANINGTOWER_PROGRESS__

class Progress {
 public:
    Progress() : m_fraction(0) {}
    void Go(int, int, int=0);
 private:
    int m_fraction;
};

inline void Progress::Go(int i, int n, int f) { // some progress bar
    int fraction = static_cast<int>(100.0 * ( i - f ) / n );
    if ( fraction > m_fraction ) {
        m_fraction = fraction;
        if ( fraction >= 10 )
            m_fraction = m_fraction/10*10 + 9;
        std::cout << fraction << "% complete: entry "<< i << std::endl;
    }
}

#endif
