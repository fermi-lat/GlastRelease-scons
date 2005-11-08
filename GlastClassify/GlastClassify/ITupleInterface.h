/** @file ITupleInterface.h
@brief  Declare abstract class ITupleInterface

$Header$
*/


#ifndef GlastClassify_ITupleInterface_h
#define GlastClassify_ITupleInterface_h

#include <string>

namespace GlastClassify{


    /** @class Item

    Nested abstract class definition
    */
    class Item {
    public:
        virtual ~Item(){};
        virtual operator double()const = 0;
    protected:
        Item(){};
    };

    /** @class ITupleInterface
    @brief nested class definition interface that must be implemented by client
    */
    class ITupleInterface 
    {
    public:
        virtual ~ITupleInterface(){}

        // set an object that evaluates at run-time to the current value of the item
        virtual const Item* getItem(const std::string& name)const=0;

        /// create a new item (float only for now) in the tuple, which will take the given value
        virtual void addItem(const std::string& name, float & value){};
    protected:
        ITupleInterface(){};
    };

}
#endif

