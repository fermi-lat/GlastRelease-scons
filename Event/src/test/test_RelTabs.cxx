// $Header$


#include "Event/RelTable/Relation.h"
#include "Event/RelTable/RelTable.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectList.h"
#include <vector>
#include <iostream>
#include <string>

using Event::RelTable;
using Event::Relation; 





// First I define two dummy classes to be used to test a relational table

static const CLID CLID_FakeOne = 5101;
static const CLID CLID_FakeTwo = 5102;



class FakeOne: public ContainedObject {
public:

    virtual const CLID& clID() const   { return FakeOne::classID(); }
    static const CLID& classID()       { return CLID_FakeOne; }

  FakeOne(std::string s, int number): name(s), ssn(number) {}
  
  // Name
  std::string name;
  // Social Securuty Number
  int ssn;
};

class FakeTwo: public ContainedObject {
public:

    virtual const CLID& clID() const   { return FakeTwo::classID(); }
    static const CLID& classID()       { return CLID_FakeTwo; }

  FakeTwo(std::string add, int fl): address(add), floor(fl) {}
  
  // Address of an Hotel
  std::string address;
  // Floor in the Hotel
  int floor;
};


int main() 
{
  FakeOne *person1, *person2, *person3;
  FakeTwo *location1, *location2;
  

  // First we create some contained object 
  person1 = new FakeOne("Jack Tripper", 654);
  person2 = new FakeOne("Richard Kanningam", 456);
  person3 = new FakeOne("Dana Scully", 231);

  location1 = new FakeTwo("531 Stanford Avenue", 2);
  location2 = new FakeTwo("520 Cowper Street", 3);


  // Then we relate them
  Relation<FakeOne, FakeTwo> *rel1 = new Relation<FakeOne, FakeTwo>(person1, location1, "06/08/2002");
  Relation<FakeOne, FakeTwo> *rel2 = new Relation<FakeOne, FakeTwo>(person2, location1, "06/08/2002");
  Relation<FakeOne, FakeTwo> *rel3 = new Relation<FakeOne, FakeTwo>(person3, location2, "10/08/2002");
  Relation<FakeOne, FakeTwo> *rel4 = new Relation<FakeOne, FakeTwo>(person1, location2, "20/09/2002");
  Relation<FakeOne, FakeTwo> *rel5 = new Relation<FakeOne, FakeTwo>(person2, location1, "09/10/2002");

  // Using the TDS, we should only require an ObjectList of Relations, ie:
  // typdef ObjectList< Relation<FakeOne, FakeTwo> > ListRelations;
  // ListRelations* rels = SmartDataPtr<ListRelations > (SomeSvc(), "/SomePath");


  // We add them to a table of relations
  RelTable<FakeOne, FakeTwo> tab;
  tab.init();

  // The table is empty and so the following query should return an emtpy vector
  std::vector< Relation<FakeOne, FakeTwo>* > empty = tab.getRelByFirst(person1);
  std::cout << "Querying and empty table the size of the returned vector is " 
            << empty.size() << std::endl;  
  
  
  tab.addRelation(rel1);
  tab.addRelation(rel2);
  tab.addRelation(rel3);
  tab.addRelation(rel4);
  tab.addRelation(rel5);

  // Using the TDS, the table is directly initialized by the ObjectList of relations: 
  // RelTable<FakeOne, FakeTwo> tab(rels);

  // Only to verify if the fillStream method works
  rel1->fillStream(std::cout);

  // Now lets do some queries

  std::vector< Relation<FakeOne, FakeTwo>* >::iterator i;

  // What are all hotels visited by person1 (Jack Tripper)?
  std::vector< Relation<FakeOne, FakeTwo>* > locs = tab.getRelByFirst(person1);
  
  std::cout << std::endl << person1->name << std::endl;
  for (i = locs.begin(); i != locs.end(); i++)
    {
      std::cout << "Address: "    << (*i)->getSecond()->address
                << "     Floor: " << (*i)->getSecond()->floor   
                << "     Date : " << (*i)->getInfos()[0]     << std::endl;
    }

  // And all persons that visited location2 ?
  std::vector< Relation<FakeOne, FakeTwo>* > pers = tab.getRelBySecond(location2);
   
  std::cout << std::endl << location2->address << std::endl;
  for (i = pers.begin(); i != pers.end(); i++)
    {
      std::cout << "Name: "     << (*i)->getFirst()->name 
                << "     ssn: " << (*i)->getFirst()->ssn   << std::endl;
    }


  // Now we change a relation already inserted
 

  locs = tab.getRelByFirst(person3);
  for (i = locs.begin(); i != locs.end(); i++)
    {
      tab.changeSecond(*i,location1);
    }
  
  std::cout << std::endl << "persion3 location changed" << std::endl;
  std::cout << std::endl << location2->address << std::endl;

  pers = tab.getRelBySecond(location1);
  for (i = pers.begin(); i != pers.end(); i++)
    {
      std::cout << "Name: "     << (*i)->getFirst()->name 
                << "     ssn: " << (*i)->getFirst()->ssn   << std::endl;
    }


  
  locs = tab.getRelByFirst(person1);
 
  for (i = locs.begin(); i != locs.end(); i++)
    {
      tab.changeFirst(*i,0);
    }
  
  std::cout << std::endl << "Jack Tripper set to null"  << std::endl;  

  // Let's see the result oh the previous query after this change. We have to
  // take care about null pointers now.

  pers = tab.getRelBySecond(location2);

  std::cout << "Number of relations with location2 = " << pers.size() << std::endl;   

  std::cout << std::endl << location2->address << std::endl;
  for (i = pers.begin(); i != pers.end(); i++)
    {
      if ((*i)->getFirst())
        std::cout << "Name: "     << (*i)->getFirst()->name 
                  << "     ssn: " << (*i)->getFirst()->ssn   << std::endl;
    }
  

  // Now let's remove the two relations with the null pointer

  pers = tab.getRelByFirst(0);
  for (i = pers.begin(); i != pers.end(); i++)
  {
    tab.erase(*i);
  }

  std::cout << std::endl << "Removed relations with null pointer" << std::endl;   
  std::cout << "Number of relations = " << tab.size() << std::endl;   
  std::cout << std::endl << location2->address << std::endl;

  pers = tab.getRelBySecond(location1);
  for (i = pers.begin(); i != pers.end(); i++)
    {
      std::cout << "Name: "     << (*i)->getFirst()->name 
                << "     ssn: " << (*i)->getFirst()->ssn   << std::endl;
    }



 std::cout << std::endl << person2->name << std::endl;


  unsigned int index;
  // Now, we verify how many times person2 has been in a location
  locs = tab.getRelByFirst(person2);
  for (i = locs.begin(); i != locs.end(); i++)
    {
      std::cout << "Address: " << (*i)->getSecond()->address << std::endl;
      std::cout << "Floor: " << (*i)->getSecond()->floor << std::endl;  
      std::cout << "Dates : " << std::endl;
      for ( index = 0; index < (*i)->getInfos().size(); index++)
        {
          std::cout << (*i)->getInfos()[index] << std::endl;
        }
    }
 

  // Using the TDS, one can finally register the collection of relations:
  // SomeSvc()->registerObject("/Event/MC/RelFakeOneFakeTwo",
  //   tab.getAllRelations());
  
  return 0;
}
