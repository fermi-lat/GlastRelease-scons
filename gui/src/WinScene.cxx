//  $Header$
//  Author: T. Burnett
//
#ifdef WIN32
// Implement the WinScene widget
#include "WinScene.h"
#include "WinDraw.h"
#include "gui/Projector.h"
using namespace gui;

WinScene::WinScene(HWND win, float size, int initial_view )
: SceneControl(m_proj = new gui::Projector(m_wd= new WinDraw(win)),size, initial_view )
{
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WinScene::~WinScene()
{
    delete m_wd;
    delete m_proj;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinScene::redisplay()
// force redraw of this window--schedules a WM_PAINT message
{
    m_wd->invalidate();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void  WinScene::redisplay(HDC hdc)
// respond to WM_PAINT message to main window
{   
    m_wd->setHDC(hdc);
    SceneControl::redisplay();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinScene::print_commands(std::ostream& out)
{
   out << 
"\n"
"Mouse commands\n"
"==============================\n"
"left botton click: \t center at mouse position\n"
"left botton drag: \t select viewport\n"
"right botton click: \t display menu\n"
"\n"
"Special Keys\n"
"============\n"
"arrows:    move reference point in selected view\n"
"PgUp/PgDn: increase/decrease magnification\n"
"Home:      reset selected view\n" 
;
    SceneControl::print_commands(out);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void WinScene::mouseDown(int ix, int iy)
{
    m_x = m_wd->xnorm(ix);
    m_y = m_wd->ynorm(iy);
}
void WinScene::mouseUp(int ix, int iy)
{
    SceneControl::setViewPort(m_x,m_y, m_wd->xnorm(ix), m_wd->ynorm(iy));
    redisplay();
}

#endif
