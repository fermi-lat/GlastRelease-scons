// $Header$
//  Author: Toby Burnett


#include "gui/Scene.h"
#include <algorithm>

namespace gui {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Scene::Scene()
    :m_enabled(true)
{}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Scene::~Scene()
{
    RepList::iterator irep = m_replist.begin();
    while( irep != m_replist.end() ) {
	delete (*irep++);
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::add(DisplayRep* pRep)
//-----------------------------------------
{
    m_replist.push_back(pRep);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::remove(DisplayRep* pRep)
//--------------------------------------------
{
     RepList::iterator it;
     it= std::find(m_replist.begin(), m_replist.end(), pRep);
     if( it != m_replist.end() )	m_replist.erase(it);

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::clear()
//---------------------
{
    if(! m_enabled)  return;

    RepList::iterator irep = m_replist.begin();
    while( irep != m_replist.end() ) {
	(*irep++)->clear();
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::update()
//-----------------------
{
    if( !m_enabled) return;
    RepList::iterator irep = m_replist.begin();
    bool updated = false;
    while( irep != m_replist.end() ) {
	if( !(*irep)->hidden() ) {
	    (*irep)->update();
	    updated = true;
	}
	++irep;
    }
    //(need something like this?) if( updated) redisplay();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::enable()
{
    m_enabled = true;
    update();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::disable()
{
    m_enabled = true;
    clear();
    m_enabled = false;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Scene::draw(ViewPort* v)
//--------------------------------
{
    RepList::iterator it=m_replist.begin();
    while( it != m_replist.end() ) (*it++)->draw(v);
}

} // namespace gui
