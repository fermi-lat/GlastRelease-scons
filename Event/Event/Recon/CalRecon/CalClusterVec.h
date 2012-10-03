#ifndef CalClusterVec_H
#define CalClusterVec_H


#include <vector>
#include <string>
#include "Event/Recon/CalRecon/CalCluster.h"


/** 
* @class CalClusterVec
*
* @brief TDS class representing a std::vector of pointers to Event::CalCluster
* obects.
*
* @author Luca Baldini.
*/


namespace Event { //Namespace Event

  class CalClusterVec : virtual public std::vector<CalCluster*>
    {
    public:
      CalClusterVec();
      virtual ~CalClusterVec() {};

      /// Overloaded push_back() method (we have to keep track of the
      /// clusters with the highest energy/gamma probability).
      void push_back(CalCluster*);

      /// Access methods.
      inline CalCluster* getHighestEnergyCluster()
        const { return m_highestEnergyCluster; }
      inline CalCluster* getHighestGamProbCluster()
        const { return m_highestGamProbCluster; }
      
    private:
      /// Pointer to the cluster in the vector with the highest energy.
      CalCluster* m_highestEnergyCluster;
      /// Pointer to the cluster in the vector with the highest gam prob.
      CalCluster* m_highestGamProbCluster;
  };


}; //Namespace Event

#endif        

