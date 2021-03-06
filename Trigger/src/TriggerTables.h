/** @file TriggerTables.h
  *  @brief Declaration of the class TriggerTables
  *
  *  $Header$
*/

#ifndef Trigger_TriggerTables_h
#define Trigger_TriggerTables_h

#include "Engine.h"
#include <vector>
#include <iostream>

namespace Trigger {

/** @class Engine
    @brief encapsulate the data and operation of the table-driven trigger scheme

    @author T. Burnett <tburnett@u.washington.edu>

    This is a list of engines
*/

class TriggerTables : public std::vector<Engine> {
public:

    /// ctor -- expect to create the engines
    /// @param configuration specify a trigger table 
    /// @param prescale vector of prescale factors: if empty, use default in table
    TriggerTables(std::string configuration, const std::vector<int>& prescale);

    /// for a gltword, return associated engine, or -1 if disabled
    int operator()(int gltword)const;


    /// make a table of the current trigger table
    void print(std::ostream& out = std::cout)const;

private:
    int engineNumber(int gltword)const;
    std::vector<const Engine*> m_table; ///< table of engine  for a bit pattern
};
}

#endif
