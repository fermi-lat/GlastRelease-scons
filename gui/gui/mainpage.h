// $Header$
// Mainpage for doxygen

/*! \mainpage package gui
The primary user interface is through a singleton GuiMgr object. It manages:
- gui::GUI the primary interface to the window system. GUI is abstract, must be actually implemented by 
      either WinGUI or MotifGUI
- gui::Menu - Encapsulates the menu bar for an application, a collection of pull-down menus.
- gui::SubMenu - Represents a pull-down menu, GuiMgr maintains one for the File menu.
- gui::DisplayControl - Manages 3-D graphics
- gui::PrintControl- Allows user control of printout to an associated console window

<br>
See the <a href=../index.html>Document</a> describing its use.
  <hr>
  \section notes release notes
  \include release.notes
  \section requirements requirements file
  \include requirements
  \section ChangeLog ChangeLog
  \include ChangeLog

*/

