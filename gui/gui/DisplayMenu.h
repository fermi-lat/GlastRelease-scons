// $Header$
#ifndef D0VE_DISPLAYMENU_HPP
#define D0VE_DISPLAYMENU_HPP
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "gui/DisplayControl.h"
#include "gui/DisplayRep.h"

namespace gui {

// class that should be below, but KAI compiler chokes
class AnonymousRep : public DisplayRep { public: void update(){}; void clear(){};};

template <class Element, class Key, class Rep>class DisplayMenu 
// Purpose: combine static display elements according to some criterion,
//          with individual show/hide capability on a named button.
// template arguments:
//     Element: a class with 3D info to display
//     Key:     a class that classifies instances of Element: must have constructor
//              Key(const Element&), and be usable as a key: e.g., an operator int(),
//              and a member function std::string name(). (Which could be the key.)
//     Rep:     a DisplayRep that has constructor Rep(const Element&);
//  all Elements with the same Key will be appended to a rep on the button with label name.
{
public:
    DisplayMenu(DisplayControl::SubMenu& menu):_menu(menu),_count(0){}
    // ctor: specify menu to add buttons to.

    void add(const Element& element)
    // add a new instance of the class Element to the display. Will create a new button
    // if Key returns a new value, otherwist append its display to previous with that
    // value.
    {
//KAI chokes: cannot be here
//class AnonymousRep : public DisplayRep { public: void update(){}; void clear(){};};
        ++_count;
        Key key(element);
        if( _map.find(key)==_map.end()) {
            _map[key]=new AnonymousRep;
            _menu.add(_map[key], key.name());
        }
	Rep r(element);
        _map[key]->append(r);
    }
    int count()const {return _count;}
    // number of times added to

    int size()const { return _map.size();}
    // number of different keys found
private:
    std::map<Key, DisplayRep*> _map;
    DisplayControl::SubMenu& _menu;
    int _count;
};

} // namespace gui
#endif // D0VE_DISPLAYMENU_HPP


