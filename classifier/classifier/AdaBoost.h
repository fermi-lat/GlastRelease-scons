/** @file  AdaBoost.h
    @brief declaration of class AdaBoost

    $Header$

*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifndef AdaBoost_h
#define AdaBoost_h
#include "classifier/Classifier.h"

/** @class AdaBoost
@brief Implement the AdaBoost boosting algorithm


*/
class AdaBoost {
public:

    /** @brief constructor 
    */
    AdaBoost(Classifier::Table& data, double beta=0.5);
    double operator()(const Classifier& tree);
    Classifier::Table& data() {return m_data;};
    ~AdaBoost();

    static inline void set_purity (double p)
    {
        s_purity = p;
    }

private:
    Classifier::Table& m_data;
    double m_beta;
	static double s_purity;
};
#endif

