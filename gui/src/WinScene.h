//     $Id$
//  Author: T. Burnett
//
// It declares the class WinScene, a specific  display
// using Windows

#define STRICT
#include <windows.h>

#include "gui/SceneControl.h"
#include "gui/Projector.h"
using namespace gui;



class WinDraw;

class WinScene  : public gui::SceneControl
{
public:

    WinScene(HWND win, float size=100, int initial_view=1);
    // constuctor

    ~WinScene();
    // dtor

    virtual void redisplay();
    // set to redisplay

    void redisplay(HDC hdc);
    // special for windows

    virtual void print_commands(std::ostream& );
    // override to describe mouse commands, special keys

    void mouseDown( int, int);
    void mouseUp( int, int);
    // process mouse click
    // process mouse click

private:

    WinDraw * m_wd;
    float m_x,m_y; // coords of mouse down if dragging
    // associated Windraw object

    gui::Projector* m_proj; // the projector object for this scene


};


