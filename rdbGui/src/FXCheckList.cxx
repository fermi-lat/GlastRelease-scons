/******************************************************************************
*                                                                             *
*                            L i s t   O b j e c t                            *
*                                                                             *
*******************************************************************************
* Copyright (C) 1997,2001 by Jeroen van der Zijp.   All Rights Reserved.      *
*******************************************************************************
* This library is free software; you can redistribute it and/or               *
* modify it under the terms of the GNU Lesser General Public                  *
* License as published by the Free Software Foundation; either                *
* version 2.1 of the License, or (at your option) any later version           *
*                                                                             *
* This library is distributed in the hope that it will be useful              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           *
* Lesser General Public License for more details.                             *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public            *
* License along with this library; if not, write to the Free Software         *
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.  *
*******************************************************************************
* $Id$                *
******************************************************************************/
#include "xincs.h"
#include "fxver.h"
#include "fxdefs.h"
#include "fxkeys.h"
#include "FXHash.h"
#include "FXThread.h"
#include "FXStream.h"
#include "FXString.h"
#include "FXSize.h"
#include "FXPoint.h"
#include "FXRectangle.h"
#include "FXRegistry.h"
#include "FXApp.h"
#include "FXDCWindow.h"
#include "FXFont.h"
#include "FXIcon.h"
#include "FXScrollBar.h"

#include "FXCheckList.h"

namespace FX{

/*

  Notes:
  - Draw stuff differently when disabled.
  - Move active item by means of focus navigation, or first letters typed.
  - PageUp/PageDn should also change currentitem.
  - Should support borders in ScrollWindow & derivatives.
  - Should issue callbacks from keyboard also.
  - Need distinguish various callbacks better:
     - Selection changes (all selected/unselected items or just changes???)
     - Changes of currentitem
     - Clicks on currentitem
  - Return key simulates double click.
  - Autoselect mode selects whatever is under the cursor, and gives callback
    on mouse release.
  - Sortfunc's will be hard to serialize, and hard to write w/o secretly 
    #including the FXCheckListItem header!
  - Should replace hard-wired SIDE_SPACING etc by flexible padding; we'll
    do this when FXScrollArea gets derived from FXFrame, so we'll have internal
    padding, and inter-item padding.
  - When in single/browse mode, current item always the one that's selected
  - Get rid of marking stuff.
  - Generate callbacks when selected, deselected, etc.
  - Add flag (default=false) to issue callback for setCurrentItem() etc.
  - It may be convenient to have ways to move items around.
  - Need insertSorted() API to add item in the right place based on current
    sort function.
  - Need to add support for arbitrary icon sizes, as FXTreeList already has.
  - Need to add support for justification (similar to FXTableItem).  Perhaps
    multi-line labels too while we're at it.
*/


#define ICON_SPACING             4    // Spacing between icon and label
#define SIDE_SPACING             6    // Left or right spacing between items
#define LINE_SPACING             4    // Line spacing between items

#define SELECT_MASK (CHECKLIST_SINGLESELECT|CHECKLIST_BROWSESELECT)
#define LIST_MASK   (SELECT_MASK|CHECKLIST_AUTOSELECT)



/*****************************************************************************/


// Object implementation
FXIMPLEMENT(FXCheckListItem,FXObject,NULL,0)


// Draw item
  void FXCheckListItem::draw(const FXCheckList* list,FXDC& dc,FXint x,FXint y,
                             FXint w,FXint h){
  FXint ih=0,th=0;
  if(icon) ih=icon->getHeight();
  if(!label.empty()) th=list->getFont()->getFontHeight();
  if(isSelected())
    dc.setForeground(list->getSelBackColor());
  else
    dc.setForeground(list->getBackColor());
  dc.fillRectangle(x,y,w,h);
  if(hasFocus()){
    drawFocus(list,dc,x,y,w,h);
  }
  x+=SIDE_SPACING/2;

  FXint ix,iy;
  FXbool check=(state&CHECKED)!=0;
  ix=x; iy=y+(h-13)/2;

  dc.setForeground(list->getApp()->getShadowColor());
  dc.fillRectangle(ix,iy,12,1);
  dc.fillRectangle(ix,iy,1,12);

  dc.setForeground(list->getApp()->getBorderColor());
  dc.fillRectangle(ix+1,iy+1,10,1);
  dc.fillRectangle(ix+1,iy+1,1,10);

  dc.setForeground(list->getApp()->getHiliteColor());
  dc.fillRectangle(ix,iy+12,13,1);
  dc.fillRectangle(ix+12,iy,1,13);

  dc.setForeground(list->getApp()->getBaseColor());
  dc.fillRectangle(ix+1,iy+11,11,1);
  dc.fillRectangle(ix+11,iy+1,1,11);

  if(!isEnabled())
    dc.setForeground(list->getApp()->getBaseColor());
  else
    dc.setForeground(list->getBackColor());
  dc.fillRectangle(ix+2,iy+2,9,9);

  if(check!=FALSE){
    FXSegment seg[6];
#ifndef WIN32
    seg[0].x1=3+ix; seg[0].y1=5+iy; seg[0].x2=5+ix; seg[0].y2=7+iy;
    seg[1].x1=3+ix; seg[1].y1=6+iy; seg[1].x2=5+ix; seg[1].y2=8+iy;
    seg[2].x1=3+ix; seg[2].y1=7+iy; seg[2].x2=5+ix; seg[2].y2=9+iy;
    seg[3].x1=5+ix; seg[3].y1=7+iy; seg[3].x2=9+ix; seg[3].y2=3+iy;
    seg[4].x1=5+ix; seg[4].y1=8+iy; seg[4].x2=9+ix; seg[4].y2=4+iy;
    seg[5].x1=5+ix; seg[5].y1=9+iy; seg[5].x2=9+ix; seg[5].y2=5+iy;
#else
    seg[0].x1=3+ix; seg[0].y1=5+iy; seg[0].x2=5+ix; seg[0].y2=7+iy;
    seg[1].x1=3+ix; seg[1].y1=6+iy; seg[1].x2=5+ix; seg[1].y2=8+iy;
    seg[2].x1=3+ix; seg[2].y1=7+iy; seg[2].x2=5+ix; seg[2].y2=9+iy;
    seg[3].x1=5+ix; seg[3].y1=7+iy; seg[3].x2=10+ix; seg[3].y2=2+iy;
    seg[4].x1=5+ix; seg[4].y1=8+iy; seg[4].x2=10+ix; seg[4].y2=3+iy;
    seg[5].x1=5+ix; seg[5].y1=9+iy; seg[5].x2=10+ix; seg[5].y2=4+iy;
#endif
    if(isEnabled()){
      if(check==MAYBE)
        dc.setForeground(list->getApp()->getShadowColor());
      else
        dc.setForeground(list->getTextColor());
    }
    else{
      dc.setForeground(list->getApp()->getShadowColor());
    }
    dc.drawLineSegments(seg,6);
  }

  x+=13+SIDE_SPACING/2;
  ///////////////////////
  if(icon){
    dc.drawIcon(icon,x,y+(h-ih)/2);
    x+=ICON_SPACING+icon->getWidth();
  }

  if(!label.empty()){
    dc.setFont(list->getFont());
    if(!isEnabled())
      dc.setForeground(makeShadowColor(list->getBackColor()));
    else if(state&SELECTED)
      dc.setForeground(list->getSelTextColor());
    else
      dc.setForeground(list->getTextColor());
    dc.drawText(x,y+(h-th)/2+list->getFont()->getFontAscent(),label.text(),label.length());
  }
}


  // Draw dotted rectangle for focus
  void FXCheckListItem::drawFocus(const FXCheckList* list,FXDC& dc,FXint x,FXint y,FXint w,FXint h) const {
    dc.setFillStyle(FILL_OPAQUESTIPPLED);
    dc.setStipple(STIPPLE_GRAY,x,y);
    dc.setForeground(list->getTextColor());
    dc.setBackground(list->getBackColor());
    dc.drawRectangle(x+1,y+1,w-3,h-3);
    dc.setFillStyle(FILL_SOLID);
    dc.setStipple(STIPPLE_NONE);
  }



  FXint FXCheckListItem::hitItem(const FXCheckList* list,FXint x,FXint y) const {
    register FXint iw=0,ih=0,tw=0,th=0,ix,iy,cx,tx,ty,h;
    register FXFont *font=list->getFont();
    if(icon){
      iw=icon->getWidth();
      ih=icon->getHeight();
    }
    if(!label.empty()){
      tw=4+font->getTextWidth(label.text(),label.length());
      th=4+font->getFontHeight();
    }
    h=LINE_SPACING+FXMAX(th,ih);
    ix=SIDE_SPACING/2;
    if (state & VISIBLE_CHECKBOX)
      cx=2*ix+13;
    else
      cx = ix;
    tx=0;
    if(iw) tx+=iw+ICON_SPACING;
    iy=(h-ih)/2;
    ty=(h-th)/2;

    // In check?
    if(ix<=x && ty<=y && x<cx && y<ty+th) return 3;  
  
    // In icon?
    if(cx<=x && iy<=y && x<cx+tx && y<iy+ih) return 1;

    // In text?
    if(cx<=x && ty<=y && x<cx+tx+tw && y<ty+th) return 2;

    // Outside
    return 0;
  }  
  
  // Set or kill focus
  void FXCheckListItem::setFocus(FXbool focus){
    if(focus) state|=FOCUS; else state&=~FOCUS;
  }

  // Select or deselect item
  void FXCheckListItem::setSelected(FXbool selected){
    if(selected) state|=SELECTED; else state&=~SELECTED;
  }

  // Check or uncheck item
  void FXCheckListItem::setChecked(FXbool checked){
    if(checked) state|=CHECKED; else state&=~CHECKED;
  }


  // Enable or disable the item
  void FXCheckListItem::setEnabled(FXbool enabled){
    if(enabled) state&=~DISABLED; else state|=DISABLED;
  }

  // Enable or disable the checkbox widget
  void FXCheckListItem::setCheckboxVisible(FXbool visible){
    if (visible) state |= VISIBLE_CHECKBOX; else state &= ~VISIBLE_CHECKBOX;
  }


  // Icon is draggable
  void FXCheckListItem::setDraggable(FXbool draggable){
    if(draggable) state|=DRAGGABLE; else state&=~DRAGGABLE;
  }


  // Icons owner by item
  void FXCheckListItem::setIconOwned(FXuint owned){
    state=(state&~ICONOWNED)|(owned&ICONOWNED);
  }


  // Create icon
  void FXCheckListItem::create(){
    if(icon) icon->create();
  }


  // Destroy icon
  void FXCheckListItem::destroy(){
    if((state&ICONOWNED) && icon) icon->destroy();
  }


  // Detach from icon resource
  void FXCheckListItem::detach(){
    if(icon) icon->detach();
  }


  // Get width of item
  FXint FXCheckListItem::getWidth(const FXCheckList* list) const {
    FXint w=0;
    if(icon) w=icon->getWidth();
    if(!label.empty()){ w+=list->getFont()->getTextWidth(label.text(),label.length()); if(icon) w+=ICON_SPACING; }
    return SIDE_SPACING+w;
  }


  // Get height of item
  FXint FXCheckListItem::getHeight(const FXCheckList* list) const {
    FXint th,ih;
    th=ih=0;
    if(!label.empty()) th=list->getFont()->getFontHeight();
    if(icon) ih=icon->getHeight();
    return LINE_SPACING+FXMAX(th,ih);
  }


  // Save data
  void FXCheckListItem::save(FXStream& store) const {
    FXObject::save(store);
    store << label;
    store << icon;
    store << state;
  }


  // Load data
  void FXCheckListItem::load(FXStream& store){
    FXObject::load(store);
    store >> label;
    store >> icon;
    store >> state;
  }


  // Delete icon if owned
  FXCheckListItem::~FXCheckListItem(){
    if(state&ICONOWNED) delete icon;
  }


  /*******************************************************************************/


  // Map
  FXDEFMAP(FXCheckList) FXCheckListMap[]={
    FXMAPFUNC(SEL_PAINT,0,FXCheckList::onPaint),
    FXMAPFUNC(SEL_ENTER,0,FXCheckList::onEnter),
    FXMAPFUNC(SEL_LEAVE,0,FXCheckList::onLeave),
    FXMAPFUNC(SEL_MOTION,0,FXCheckList::onMotion),
    FXMAPFUNC(SEL_TIMEOUT,FXWindow::ID_AUTOSCROLL,FXCheckList::onAutoScroll),
    FXMAPFUNC(SEL_TIMEOUT,FXCheckList::ID_TIPTIMER,FXCheckList::onTipTimer),
    FXMAPFUNC(SEL_TIMEOUT,FXCheckList::ID_LOOKUPTIMER,FXCheckList::onLookupTimer),
    FXMAPFUNC(SEL_COMMAND,FXCheckList::ID_CHECKSEL,FXCheckList::onCheckSel),
    FXMAPFUNC(SEL_COMMAND,FXCheckList::ID_UNCHECKSEL,FXCheckList::onUncheckSel),
    FXMAPFUNC(SEL_COMMAND,FXCheckList::ID_CHECKALL,FXCheckList::onCheckAll),
    FXMAPFUNC(SEL_COMMAND,FXCheckList::ID_UNCHECKALL,FXCheckList::onUncheckAll),
    FXMAPFUNC(SEL_UNGRABBED,0,FXCheckList::onUngrabbed),
    FXMAPFUNC(SEL_LEFTBUTTONPRESS,0,FXCheckList::onLeftBtnPress),
    FXMAPFUNC(SEL_LEFTBUTTONRELEASE,0,FXCheckList::onLeftBtnRelease),
    FXMAPFUNC(SEL_RIGHTBUTTONPRESS,0,FXCheckList::onRightBtnPress),
    FXMAPFUNC(SEL_RIGHTBUTTONRELEASE,0,FXCheckList::onRightBtnRelease),
    FXMAPFUNC(SEL_KEYPRESS,0,FXCheckList::onKeyPress),
    FXMAPFUNC(SEL_KEYRELEASE,0,FXCheckList::onKeyRelease),
    FXMAPFUNC(SEL_FOCUSIN,0,FXCheckList::onFocusIn),
    FXMAPFUNC(SEL_FOCUSOUT,0,FXCheckList::onFocusOut),
    FXMAPFUNC(SEL_CLICKED,0,FXCheckList::onClicked),
    FXMAPFUNC(SEL_DOUBLECLICKED,0,FXCheckList::onDoubleClicked),
    FXMAPFUNC(SEL_TRIPLECLICKED,0,FXCheckList::onTripleClicked),
    FXMAPFUNC(SEL_COMMAND,0,FXCheckList::onCommand),
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    FXMAPFUNC(SEL_UPDATE,FXWindow::ID_QUERY_TIP,FXCheckList::onQueryTip),
    FXMAPFUNC(SEL_UPDATE,FXWindow::ID_QUERY_HELP,FXCheckList::onQueryHelp),
#elif FOX_MAJOR >=1 && FOX_MINOR >=4
    FXMAPFUNC(SEL_QUERY_TIP,0,FXCheckList::onQueryTip),
    FXMAPFUNC(SEL_QUERY_HELP,0,FXCheckList::onQueryHelp),
#endif

    FXMAPFUNC(SEL_COMMAND,FXWindow::ID_SETVALUE,FXCheckList::onCmdSetValue),
    FXMAPFUNC(SEL_COMMAND,FXWindow::ID_SETINTVALUE,FXCheckList::onCmdSetIntValue),
    FXMAPFUNC(SEL_COMMAND,FXWindow::ID_GETINTVALUE,FXCheckList::onCmdGetIntValue),
  };


  // Object implementation
  FXIMPLEMENT(FXCheckList,FXScrollArea,FXCheckListMap,ARRAYNUMBER(FXCheckListMap))


    // List
    FXCheckList::FXCheckList(){
    flags|=FLAG_ENABLED;
    items=NULL;
    nitems=0;
    anchor=-1;
    current=-1;
    extent=-1;
    cursor=-1;
    font=(FXFont*)-1;
    textColor=0;
    selbackColor=0;
    seltextColor=0;
    itemWidth=1;
    itemHeight=1;
    visible=0;
    sortfunc=NULL;
    grabx=0;
    graby=0;
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    timer=NULL;
    lookuptimer=NULL;
#endif
    state=FALSE;
  }


  // List
  FXCheckList::FXCheckList(FXComposite *p,FXint nvis,FXObject* tgt,FXSelector sel,FXuint opts,FXint x,FXint y,FXint w,FXint h):
    FXScrollArea(p,opts,x,y,w,h){
    flags|=FLAG_ENABLED;
    target=tgt;
    message=sel;
    items=NULL;
    nitems=0;
    anchor=-1;
    current=-1;
    extent=-1;
    cursor=-1;
    font=getApp()->getNormalFont();
    textColor=getApp()->getForeColor();
    selbackColor=getApp()->getSelbackColor();
    seltextColor=getApp()->getSelforeColor();
    itemWidth=1;
    itemHeight=1;
    visible=FXMAX(nvis,0);
    sortfunc=NULL;
    grabx=0;
    graby=0;
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    timer=NULL;
    lookuptimer=NULL;
#endif
    state=FALSE;
  }


  // Create window
  void FXCheckList::create(){
    register FXint i;
    FXScrollArea::create();
    for(i=0; i<nitems; i++){items[i]->create();}
    font->create();
  }


  // Detach window
  void FXCheckList::detach(){
    register FXint i;
    FXScrollArea::detach();
    for(i=0; i<nitems; i++){items[i]->detach();}
    font->detach();
  }


  // Can have focus
#if FOX_MAJOR >=1 && FOX_MINOR < 6
  FXbool
#elif FOX_MAJOR >=1 && FOX_MINOR >= 6
  bool
#endif
  FXCheckList::canFocus() const { return TRUE; }

  // Get default width
  FXint FXCheckList::getDefaultWidth(){
    return FXScrollArea::getDefaultWidth();
  }


  // Get default height
  FXint FXCheckList::getDefaultHeight(){
    if(visible) return visible*(LINE_SPACING+font->getFontHeight());
    return FXScrollArea::getDefaultHeight();
  }


  // Propagate size change
  void FXCheckList::recalc(){
    FXScrollArea::recalc();
    flags|=FLAG_RECALC;
    cursor=-1;
  }


  // List is multiple of nitems
  void FXCheckList::setNumVisible(FXint nvis){
    if(nvis<0) nvis=0;
    if(visible!=nvis){
      visible=nvis;
      recalc();
    }
  }


  // Recompute interior
  void FXCheckList::recompute(){
    register FXint w,h,i;
    itemWidth=1;
    itemHeight=1;
    for(i=0; i<nitems; i++){
      w=items[i]->getWidth(this);
      h=items[i]->getHeight(this);
      if(w>itemWidth) itemWidth=w;
      if(h>itemHeight) itemHeight=h;
    }
    flags&=~FLAG_RECALC;
  }


  // Determine content width of list
  FXint FXCheckList::getContentWidth(){
    if(flags&FLAG_RECALC) recompute();
    return itemWidth;
  }


  // Determine content height of list
  FXint FXCheckList::getContentHeight(){
    if(flags&FLAG_RECALC) recompute();
    return nitems*itemHeight;
  }


  // Recalculate layout determines item locations and sizes
  void FXCheckList::layout(){

    // Force repaint if content changed
    //if(flags&FLAG_RECALC) update();

    // Calculate contents
    FXScrollArea::layout();

    // Determine line size for scroll bars
    vertical->setLine(itemHeight);
    horizontal->setLine(1);

    update();

    // No more dirty
    flags&=~FLAG_DIRTY;
  }


  // Change item text
  void FXCheckList::setItemText(FXint index,const FXString& text){
    if(index<0 || nitems<=index){ fxerror("%s::setItemText: index out of range.\n",getClassName()); }
    items[index]->setText(text);
    recalc();
  }


  // Get item text
  FXString FXCheckList::getItemText(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::getItemText: index out of range.\n",getClassName()); }
    return items[index]->getText();
  }


  // Set item icon
  void FXCheckList::setItemIcon(FXint index,FXIcon* icon){
    if(index<0 || nitems<=index){ fxerror("%s::setItemIcon: index out of range.\n",getClassName()); }
    items[index]->setIcon(icon);
    recalc();
  }


  // Get item icon
  FXIcon* FXCheckList::getItemIcon(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::getItemIcon: index out of range.\n",getClassName()); }
    return items[index]->getIcon();
  }


  // Set item data
  void FXCheckList::setItemData(FXint index,void* ptr){
    if(index<0 || nitems<=index){ fxerror("%s::setItemData: index out of range.\n",getClassName()); }
    items[index]->setData(ptr);
  }


  // Get item data
  void* FXCheckList::getItemData(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::getItemData: index out of range.\n",getClassName()); }
    return items[index]->getData();
  }


  // True if item is selected
  FXbool FXCheckList::isItemSelected(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::isItemSelected: index out of range.\n",getClassName()); }
    return items[index]->isSelected();
  }


  // True if item is checked
  FXbool FXCheckList::isItemChecked(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::isItemChecked: index out of range.\n",getClassName()); }
    return items[index]->isChecked();
  }


  // True if item is current
  FXbool FXCheckList::isItemCurrent(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::isItemCurrent: index out of range.\n",getClassName()); }
    return index==current;
  }


  // True if item is enabled
  FXbool FXCheckList::isItemEnabled(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::isItemEnabled: index out of range.\n",getClassName()); }
    return items[index]->isEnabled();
  }


  // True if item (partially) visible
  FXbool FXCheckList::isItemVisible(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::isItemVisible: index out of range.\n",getClassName()); }
    return (0 < (pos_y+index*itemHeight+itemHeight)) && ((pos_y+index*itemHeight) < viewport_h);
  }


  // Make item fully visible
  void FXCheckList::makeItemVisible(FXint index){
    if(xid==0) return;
    // FIXME maybe force layout first???
    if(0<=index && index<nitems){
      if((pos_y+index*itemHeight) < 0){
        setPosition(pos_x,-index*itemHeight);
      }
      else if(viewport_h<=(pos_y+index*itemHeight+itemHeight)){
        setPosition(pos_x,viewport_h-index*itemHeight-itemHeight);
      }
    }
  }


  // Get item at position x,y
  FXint FXCheckList::getItemAt(FXint,FXint y) const {
    FXint index=(y-pos_y)/itemHeight;
    if(index<0 || index>=nitems) index=-1;
    return index;
  }


  // Did we hit the item, and which part of it did we hit (0=outside, 1=icon, 2=text)
  // 3=check
  FXint FXCheckList::hitItem(FXint index,FXint x,FXint y) const {
    FXint ix,iy,hit=0;
    if(0<=index && index<nitems){
      x-=pos_x;
      y-=pos_y;
      ix=0;
      iy=itemHeight*index;
      hit=items[index]->hitItem(this,x-ix,y-iy);
    }
    return hit;
  }


  // Repaint
  void FXCheckList::updateItem(FXint index){
    if(0<=index && index<nitems){
      update(0,pos_y+index*itemHeight,viewport_w,itemHeight);
    }
  }


  // Enable one item
  FXbool FXCheckList::enableItem(FXint index){
    if(index<0 || nitems<=index){ fxerror("%s::enableItem: index out of range.\n",getClassName()); }
    if(!items[index]->isEnabled()){
      items[index]->setEnabled(TRUE);
      updateItem(index);
      return TRUE;
    }
    return FALSE;
  }


  // Disable one item
  FXbool FXCheckList::disableItem(FXint index){
    if(index<0 || nitems<=index){ fxerror("%s::disableItem: index out of range.\n",getClassName()); }
    if(items[index]->isEnabled()){
      items[index]->setEnabled(FALSE);
      updateItem(index);
      return TRUE;
    }
    return FALSE;
  }

  void FXCheckList::enableAllItems(){
    register FXint i;
    for(i=0; i<nitems; i++)
      enableItem(i);
  }
  
  
  void FXCheckList::disableAllItems(){
    register FXint i;
    for(i=0; i<nitems; i++)
      disableItem(i);
  }


  /// Make checkbox visible
  FXbool FXCheckList::setItemCheckboxVisible(FXint index){
    if(index<0 || nitems<=index){ fxerror("%s::setItemCheckboxVisible: index out of range.\n",getClassName()); }
    if(!items[index]->isCheckboxVisible()){
      items[index]->setCheckboxVisible(TRUE);
      updateItem(index);
      return TRUE;
    }
    return FALSE;
  }

  /// Make checkbox hidden
  FXbool FXCheckList::setItemCheckboxHidden(FXint index){
    if(index<0 || nitems<=index){ fxerror("%s::setItemCheckboxHidden: index out of range.\n",getClassName()); }
    if(items[index]->isCheckboxVisible()){
      items[index]->setCheckboxVisible(FALSE);
      updateItem(index);
      return TRUE;
    }
    return FALSE;
  }

  /// Make all checkbox visible
  void FXCheckList::setAllCheckboxVisible(){
    register FXint i;
    for(i=0; i<nitems; i++)
      setItemCheckboxVisible(i);
  }

  /// Make all checkbox hidden
  void FXCheckList::setAllCheckboxHidden(){
    register FXint i;
    for(i=0; i<nitems; i++)
      setItemCheckboxHidden(i);
  }

  // Select one item
  FXbool FXCheckList::selectItem(FXint index,FXbool notify){
    if(index<0 || nitems<=index){ fxerror("%s::selectItem: index out of range.\n",getClassName()); }
    if(!items[index]->isSelected()){
      switch(options&SELECT_MASK){
      case CHECKLIST_SINGLESELECT:
      case CHECKLIST_BROWSESELECT:
        killSelection(notify);
      case CHECKLIST_EXTENDEDSELECT:
      case CHECKLIST_MULTIPLESELECT:
        items[index]->setSelected(TRUE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*) index);}
        break;
      }
      return TRUE;
    }
    return FALSE;
  }


  // Deselect one item
  FXbool FXCheckList::deselectItem(FXint index,FXbool notify){
    if(index<0 || nitems<=index){ fxerror("%s::deselectItem: index out of range.\n",getClassName()); }
    if(items[index]->isSelected()){
      switch(options&SELECT_MASK){
      case CHECKLIST_EXTENDEDSELECT:
      case CHECKLIST_MULTIPLESELECT:
      case CHECKLIST_SINGLESELECT:
        items[index]->setSelected(FALSE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*)index);}
        break;
      }
      return TRUE;
    }
    return FALSE;
  }


  // check one item
  FXbool FXCheckList::checkItem(FXint index,FXbool notify){
    if(index<0 || nitems<=index){ fxerror("%s::checkItem: index out of range.\n",getClassName()); }
    if(!items[index]->isChecked()){
      items[index]->setChecked(TRUE);
      updateItem(index);
      if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*) index);}
      return TRUE;
    }
    return FALSE;
  }


  // uncheck one item
  FXbool FXCheckList::uncheckItem(FXint index,FXbool notify){
    if(index<0 || nitems<=index){ fxerror("%s::uncheckItem: index out of range.\n",getClassName()); }
    if(items[index]->isChecked()){
      items[index]->setChecked(FALSE);
      updateItem(index);
      if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*) index);}
      return TRUE;
    }
    return FALSE;
  }


  // Toggle check item
  FXbool FXCheckList::toggleCheckItem(FXint index,FXbool notify){
    if(index<0 || nitems<=index){ fxerror("%s::toggleCheckItem: index out of range.\n",getClassName()); }
    if(items[index]->isChecked()){
      items[index]->setChecked(FALSE);
      updateItem(index);
      if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*) index);}
      return TRUE;
    }
    else{
      items[index]->setChecked(TRUE);
      updateItem(index);
      if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*) index);}
      return TRUE;
    }
    return FALSE;
  }


  // Toggle one item
  FXbool FXCheckList::toggleItem(FXint index,FXbool notify){
    if(index<0 || nitems<=index){ fxerror("%s::toggleItem: index out of range.\n",getClassName()); }
    switch(options&SELECT_MASK){
    case CHECKLIST_BROWSESELECT:
      if(!items[index]->isSelected()){
        killSelection(notify);
        items[index]->setSelected(TRUE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*) index);}
      }
      break;
    case CHECKLIST_SINGLESELECT:
      if(!items[index]->isSelected()){
        killSelection(notify);
        items[index]->setSelected(TRUE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*) index);}
      }
      else{
        items[index]->setSelected(FALSE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*) index);}
      }
      break;
    case CHECKLIST_EXTENDEDSELECT:
    case CHECKLIST_MULTIPLESELECT:
      if(!items[index]->isSelected()){
        items[index]->setSelected(TRUE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*) index);}
      }
      else{
        items[index]->setSelected(FALSE);
        updateItem(index);
        if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*) index);}
      }
      break;
    }
    return TRUE;
  }


  // Update value from a message
  long FXCheckList::onCmdSetValue(FXObject*,FXSelector,void* ptr){
    setCurrentItem((FXint)(long)ptr);
    return 1;
  }


  // Obtain value from list
  long FXCheckList::onCmdGetIntValue(FXObject*,FXSelector,void* ptr){
    *((FXint*)ptr)=getCurrentItem();
    return 1;
  }


  // Update value from a message
  long FXCheckList::onCmdSetIntValue(FXObject*,FXSelector,void* ptr){
    setCurrentItem(*((FXint*)ptr));
    return 1;
  }


  // Enter window
  long FXCheckList::onEnter(FXObject* sender,FXSelector sel,void* ptr){
    FXScrollArea::onEnter(sender,sel,ptr);
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    if(!timer) {
      timer=getApp()->addTimeout(this,ID_TIPTIMER,getApp()->getMenuPause());
    }
#endif
    cursor=-1;
    return 1;
  }


  // Leave window
  long FXCheckList::onLeave(FXObject* sender,FXSelector sel,void* ptr){
    FXScrollArea::onLeave(sender,sel,ptr);
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    if(timer){getApp()->removeTimeout(timer);timer=NULL;}
#elif FOX_MAJOR >=1 && FOX_MINOR >= 4
    getApp()->removeTimeout(this, ID_TIPTIMER);
#endif
    cursor=-1;
    return 1;
  }


  // Gained focus
  long FXCheckList::onFocusIn(FXObject* sender,FXSelector sel,void* ptr){
    FXScrollArea::onFocusIn(sender,sel,ptr);
    if(0<=current){
      FXASSERT(current<nitems);
      items[current]->setFocus(TRUE);
      updateItem(current);
    }
    return 1;
  }


  // Lost focus
  long FXCheckList::onFocusOut(FXObject* sender,FXSelector sel,void* ptr){
    FXScrollArea::onFocusOut(sender,sel,ptr);
    if(0<=current){
      FXASSERT(current<nitems);
      items[current]->setFocus(FALSE);
      updateItem(current);
    }
    return 1;
  }


  // Draw item list
  long FXCheckList::onPaint(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    FXDCWindow dc(this,event);
    FXint indexlo=(event->rect.y-pos_y)/itemHeight;
    FXint indexhi=(event->rect.y+event->rect.h-pos_y)/itemHeight;
    FXint i,y;
    if(indexlo<0) indexlo=0;
    if(indexhi>=nitems) indexhi=nitems-1;
    y=pos_y+indexlo*itemHeight;
    for(i=indexlo; i<=indexhi; i++){
      items[i]->draw(this,dc,pos_x,y,FXMAX(itemWidth,viewport_w),itemHeight);
      y+=itemHeight;
    }
    if(y<event->rect.y+event->rect.h){
      dc.setForeground(backColor);
      dc.fillRectangle(event->rect.x,y,event->rect.w,event->rect.y+event->rect.h-y);
    }
    return 1;
  }


  // We timed out, i.e. the user didn't move for a while
  long FXCheckList::onTipTimer(FXObject*,FXSelector,void*){
    FXTRACE((250,"%s::onTipTimer %p\n",getClassName(),this));
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    timer=NULL;
#endif
    flags|=FLAG_TIP;
    return 1;
  }


  // Zero out lookup string
  long FXCheckList::onLookupTimer(FXObject*,FXSelector,void*){
    lookup=FXString::null;
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    lookuptimer=NULL;
#endif
    return 1;
  }

  // Check selected items
  long FXCheckList::onCheckSel(FXObject*,FXSelector,void*){
    for (FXint i = 0; i < nitems; i++){
      if (isItemSelected(i))      checkItem(i, true);  
    }

    return 1;
  }
  
  // Uncheck selected items  
  long FXCheckList::onUncheckSel(FXObject*,FXSelector,void*){
    for (FXint i = 0; i < nitems; i++){
      if (isItemSelected(i))      uncheckItem(i, true);  
    }

    return 1;
  }

  // Check all items
  long FXCheckList::onCheckAll(FXObject*,FXSelector,void*){
    for (FXint i = 0; i < nitems; i++){
      checkItem(i, true);
    }

    return 1;
  }

  // Uncheck all items
  long FXCheckList::onUncheckAll(FXObject*,FXSelector,void*){
    for (FXint i = 0; i < nitems; i++){
      if (isItemEnabled(i))  uncheckItem(i, true);
    }

    return 1;
  }

  // We were asked about tip text
  long FXCheckList::onQueryTip(FXObject* sender,FXSelector,void*){
    if((flags&FLAG_TIP) && !(options&CHECKLIST_AUTOSELECT)){   // No tip when autoselect!
      if(0<=cursor){
        FXString string=items[cursor]->getTipText();
        sender->handle(this,MKUINT(ID_SETSTRINGVALUE,SEL_COMMAND),(void*)&string);
        return 1;
      }
    }
    return 0;
  }


  // We were asked about status text
  long FXCheckList::onQueryHelp(FXObject* sender,FXSelector,void*){
    if(!help.empty() && (flags&FLAG_HELP)){
      FXTRACE((250,"%s::onQueryHelp %p\n",getClassName(),this));
      sender->handle(this,MKUINT(ID_SETSTRINGVALUE,SEL_COMMAND),(void*)&help);
      return 1;
    }
    return 0;
  }


  // Key Press
  long FXCheckList::onKeyPress(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    FXint index=current;
    flags&=~FLAG_TIP;
    if(!isEnabled()) return 0;
    if(target && target->handle(this,MKUINT(message,SEL_KEYPRESS),ptr)) return 1;
    if(index<0) index=0;
    switch(event->code){
    case KEY_Control_L:
    case KEY_Control_R:
    case KEY_Shift_L:
    case KEY_Shift_R:
    case KEY_Alt_L:
    case KEY_Alt_R:
      if(flags&FLAG_DODRAG){handle(this,MKUINT(0,SEL_DRAGGED),ptr);}
      return 1;
    case KEY_Page_Up:
    case KEY_KP_Page_Up:
      lookup=FXString::null;
      return 1;
    case KEY_Page_Down:
    case KEY_KP_Page_Down:
      lookup=FXString::null;
      return 1;
    case KEY_Up:
    case KEY_KP_Up:
      index-=1;
      goto hop;
    case KEY_Down:
    case KEY_KP_Down:
      index+=1;
      goto hop;
    case KEY_Home:
    case KEY_KP_Home:
      index=0;
      goto hop;
    case KEY_End:
    case KEY_KP_End:
      index=nitems-1;
    hop:  lookup=FXString::null;
      if(0<=index && index<nitems){
        setCurrentItem(index,TRUE);
        makeItemVisible(index);
        if(items[index]->isEnabled()){
          if((options&SELECT_MASK)==CHECKLIST_EXTENDEDSELECT){
            if(event->state&SHIFTMASK){
              if(0<=anchor){
                selectItem(anchor,TRUE);
                extendSelection(index,TRUE);
              }
              else{
                selectItem(index,TRUE);
                setAnchorItem(index);
              }
            }
            else if(!(event->state&CONTROLMASK)){
              killSelection(TRUE);
              selectItem(index,TRUE);
              setAnchorItem(index);
            }
          }
        }
      }
      handle(this,MKUINT(0,SEL_CLICKED),(void*)current);
      if(0<=current && items[current]->isEnabled()){
        handle(this,MKUINT(0,SEL_COMMAND),(void*)current);
      }
      return 1;
    case KEY_space:
    case KEY_KP_Space:
      lookup=FXString::null;
      if(0<=index && items[index]->isEnabled()){
        switch(options&SELECT_MASK){
        case CHECKLIST_EXTENDEDSELECT:
          if(event->state&SHIFTMASK){
            if(0<=anchor){
              selectItem(anchor,TRUE);
              extendSelection(index,TRUE);
            }
            else{
              selectItem(index,TRUE);
            }
          }
          else if(event->state&CONTROLMASK){
            toggleItem(index,TRUE);
          }
          else{
            killSelection(TRUE);
            selectItem(index,TRUE);
            toggleCheckItem(index,TRUE);
          }
          break;
        case CHECKLIST_MULTIPLESELECT:
        case CHECKLIST_SINGLESELECT:
          toggleItem(index,TRUE);
          toggleCheckItem(index,TRUE);
          break;
        }
        setAnchorItem(index);
      }
      handle(this,MKUINT(0,SEL_CLICKED),(void*)current);
      if(0<=current && items[current]->isEnabled()){
        handle(this,MKUINT(0,SEL_COMMAND),(void*)current);
      }
      return 1;
    case KEY_Return:
    case KEY_KP_Enter:
      lookup=FXString::null;
      handle(this,MKUINT(0,SEL_DOUBLECLICKED),(void*)current);
      if(0<=current && items[current]->isEnabled()){
        handle(this,MKUINT(0,SEL_COMMAND),(void*)current);
      }
      return 1;
    default:
      if((event->state&(CONTROLMASK|ALTMASK)) || ((FXuchar)event->text[0]<32)) return 0;
      lookup.append(event->text);
#if FOX_MAJOR >=1 && FOX_MINOR < 4
      if(lookuptimer) getApp()->removeTimeout(lookuptimer);
      lookuptimer=getApp()->addTimeout(this,ID_LOOKUPTIMER,getApp()->getTypingSpeed());
#elif FOX_MAJOR >=1 && FOX_MINOR >= 4
      getApp()->addTimeout(this,ID_LOOKUPTIMER,getApp()->getTypingSpeed());
#endif
      index=findItem(lookup,current,SEARCH_FORWARD|SEARCH_WRAP|SEARCH_PREFIX);
      if(0<=index){
	setCurrentItem(index,TRUE);
	makeItemVisible(index);
	if((options&SELECT_MASK)==CHECKLIST_EXTENDEDSELECT){
	  if(items[index]->isEnabled()){
	    killSelection(TRUE);
	    selectItem(index,TRUE);
          }
        }
	setAnchorItem(index);
      }
      handle(this,MKUINT(0,SEL_CLICKED),(void*)current);
      if(0<=current && items[current]->isEnabled()){
	handle(this,MKUINT(0,SEL_COMMAND),(void*)current);
      }
      return 1;
    }
    return 0;
  }


  // Key Release
  long FXCheckList::onKeyRelease(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    if(!isEnabled()) return 0;
    if(target && target->handle(this,MKUINT(message,SEL_KEYRELEASE),ptr)) return 1;
    switch(event->code){
    case KEY_Shift_L:
    case KEY_Shift_R:
    case KEY_Control_L:
    case KEY_Control_R:
    case KEY_Alt_L:
    case KEY_Alt_R:
      if(flags&FLAG_DODRAG){handle(this,MKUINT(0,SEL_DRAGGED),ptr);}
      return 1;
    }
    return 0;
  }


  // Automatic scroll
  long FXCheckList::onAutoScroll(FXObject* sender,FXSelector sel,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    FXint index;

    // Scroll the window
    FXScrollArea::onAutoScroll(sender,sel,ptr);

    // Drag and drop mode
    if(flags&FLAG_DODRAG){
      handle(this,MKUINT(0,SEL_DRAGGED),ptr);
      return 1;
    }

    // In autoselect mode, stop scrolling when mouse outside window
    if((flags&FLAG_PRESSED) || (options&CHECKLIST_AUTOSELECT)){

      // Validated position
      FXint xx=event->win_x; 

      if(xx<0) xx=0; 
      else if(xx>=viewport_w) xx=viewport_w-1;

      FXint yy=event->win_y; if(yy<0) yy=0; else if(yy>=viewport_h) yy=viewport_h-1;

      // Find item
      index=getItemAt(xx,yy);

      // Got item and different from last time
      if(0<=index && index!=current){

        // Make it the current item
        setCurrentItem(index,TRUE);

        // Extend the selection
        if((options&SELECT_MASK)==CHECKLIST_EXTENDEDSELECT){
          state=FALSE;
          extendSelection(index,TRUE);
        }
      }
      return 1;
    }
    return 0;
  }


  // Mouse moved
  long FXCheckList::onMotion(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    FXint oldcursor=cursor;
    FXuint flg=flags;

    // Kill the tip
    flags&=~FLAG_TIP;

    // Kill the tip timer
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    if(timer) timer=getApp()->removeTimeout(timer);
#elif FOX_MAJOR >=1 && FOX_MINOR >= 4
    getApp()->removeTimeout(this, ID_TIPTIMER);
#endif

    // Right mouse scrolling
    if(flags&FLAG_SCROLLING){
      setPosition(event->win_x-grabx,event->win_y-graby);
      return 1;
    }

    // Drag and drop mode
    if(flags&FLAG_DODRAG){
      if(startAutoScroll(event,TRUE)) return 1;
      handle(this,MKUINT(0,SEL_DRAGGED),ptr);
      return 1;
    }

    // Tentative drag and drop
    if((flags&FLAG_TRYDRAG) && event->moved){
      flags&=~FLAG_TRYDRAG;
      if(handle(this,MKUINT(0,SEL_BEGINDRAG),ptr)){
        flags|=FLAG_DODRAG;
      }
      return 1;
    }

    // Normal operation
    if((flags&FLAG_PRESSED) || (options&CHECKLIST_AUTOSELECT)){

      // Start auto scrolling?
      if(startAutoScroll(event,FALSE)) return 1;

      // Find item
      FXint index=getItemAt(event->win_x,event->win_y);

      // Got an item different from before
      if(0<=index && index!=current){

        // Make it the current item
        setCurrentItem(index,TRUE);

        // Extend the selection
        if((options&SELECT_MASK)==CHECKLIST_EXTENDEDSELECT){
          state=FALSE;
          extendSelection(index,TRUE);
        }
        return 1;
      }
    }

    // Reset tip timer if nothing's going on
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    timer=getApp()->addTimeout(this,ID_TIPTIMER,getApp()->getMenuPause());
#elif FOX_MAJOR >=1 && FOX_MINOR >= 4
    getApp()->addTimeout(this,ID_TIPTIMER,getApp()->getMenuPause());
#endif

    // Get item we're over
    cursor=getItemAt(event->win_x,event->win_y);

    // Force GUI update only when needed
    return (cursor!=oldcursor)||(flg&FLAG_TIP);
  }


  // Pressed a button
  long FXCheckList::onLeftBtnPress(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    FXint index,code;
    flags&=~FLAG_TIP;
    handle(this,MKUINT(0,SEL_FOCUS_SELF),ptr);
    if(isEnabled()){
      grab();
      flags&=~FLAG_UPDATE;

      // First change callback
      if(target && target->handle(this,MKUINT(message,SEL_LEFTBUTTONPRESS),ptr)) return 1;

      // Autoselect mode
      if(options&CHECKLIST_AUTOSELECT) return 1;

      // Locate item
      index=getItemAt(event->win_x,event->win_y);
      code=hitItem(index,event->win_x,event->win_y);

      // No item
      if(index<0) return 1;

      // Find out where hit
      code=hitItem(index,event->win_x,event->win_y);

      // Change current item
      setCurrentItem(index,TRUE);

      if(code==3){
        toggleCheckItem(index,TRUE);
      }
      else{
        // Change item selection
        state=items[index]->isSelected();
        switch(options&SELECT_MASK){
        case CHECKLIST_EXTENDEDSELECT:
          if(event->state&SHIFTMASK){
            if(0<=anchor){
              if(items[anchor]->isEnabled()) selectItem(anchor,TRUE);
              extendSelection(index,TRUE);
            }
            else{
              if(items[index]->isEnabled()) selectItem(index,TRUE);
              setAnchorItem(index);
            }
          }
          else if(event->state&CONTROLMASK){
            if(items[index]->isEnabled() && !state) selectItem(index,TRUE);
            setAnchorItem(index);
          }
          else{
            if(items[index]->isEnabled() && !state){ 
              killSelection(TRUE); selectItem(index,TRUE);
            }
            setAnchorItem(index);
          }
          break;
        case CHECKLIST_MULTIPLESELECT:
        case CHECKLIST_SINGLESELECT:
          if(items[index]->isEnabled() && !state){
            selectItem(index,TRUE);
          }
          break;
        }
      }

      // Start drag if actually pressed text or icon only
      if(code && items[index]->isSelected() && items[index]->isDraggable()){
        flags|=FLAG_TRYDRAG;
      }

      flags|=FLAG_PRESSED;
      return 1;
    }
    return 0;
  }


  // Released button
  long FXCheckList::onLeftBtnRelease(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    FXuint flg=flags;
    if(isEnabled()){
      ungrab();
      stopAutoScroll();
      flags|=FLAG_UPDATE;
      flags&=~(FLAG_PRESSED|FLAG_TRYDRAG|FLAG_DODRAG);

      // First chance callback
      if(target && target->handle(this,MKUINT(message,SEL_LEFTBUTTONRELEASE),ptr)) return 1;

      // No activity
      if(!(flg&FLAG_PRESSED) && !(options&CHECKLIST_AUTOSELECT)) return 1;

      // Was dragging
      if(flg&FLAG_DODRAG){
        handle(this,MKUINT(0,SEL_ENDDRAG),ptr);
        return 1;
      }

      // Selection change
      switch(options&SELECT_MASK){
      case CHECKLIST_EXTENDEDSELECT:
        if(0<=current && items[current]->isEnabled()){
          if(event->state&CONTROLMASK){
            if(state) deselectItem(current,TRUE);
          }
          else if(!(event->state&SHIFTMASK)){
            if(state){ killSelection(TRUE); selectItem(current,TRUE); }
          }
        }
        break;
      case CHECKLIST_MULTIPLESELECT:
      case CHECKLIST_SINGLESELECT:
        if(0<=current && items[current]->isEnabled()){
          if(state) deselectItem(current,TRUE);
        }
        break;
      }

      // Scroll to make item visibke
      makeItemVisible(current);

      // Update anchor
      setAnchorItem(current);

      // Generate clicked callbacks
      if(event->click_count==1){
        handle(this,MKUINT(0,SEL_CLICKED),(void*)current);
      }
      else if(event->click_count==2){
        handle(this,MKUINT(0,SEL_DOUBLECLICKED),(void*)current);
      }
      else if(event->click_count==3){
        handle(this,MKUINT(0,SEL_TRIPLECLICKED),(void*)current);
      }

      // Command callback only when clicked on item
      if(0<=current && items[current]->isEnabled()){
        handle(this,MKUINT(0,SEL_COMMAND),(void*)current);
      }
      return 1;
    }
    return 0;
  }


  // Pressed right button
  long FXCheckList::onRightBtnPress(FXObject*,FXSelector,void* ptr){
    FXEvent* event=(FXEvent*)ptr;
    flags&=~FLAG_TIP;
    handle(this,MKUINT(0,SEL_FOCUS_SELF),ptr);
    if(isEnabled()){
      grab();
      flags&=~FLAG_UPDATE;
      if(target && target->handle(this,MKUINT(message,SEL_RIGHTBUTTONPRESS),ptr)) return 1;
      flags|=FLAG_SCROLLING;
      grabx=event->win_x-pos_x;
      graby=event->win_y-pos_y;
      return 1;
    }
    return 0;
  }


  // Released right button
  long FXCheckList::onRightBtnRelease(FXObject*,FXSelector,void* ptr){
    if(isEnabled()){
      ungrab();
      flags&=~FLAG_SCROLLING;
      flags|=FLAG_UPDATE;
      if(target && target->handle(this,MKUINT(message,SEL_RIGHTBUTTONRELEASE),ptr)) return 1;
      return 1;
    }
    return 0;
  }


  // The widget lost the grab for some reason
  long FXCheckList::onUngrabbed(FXObject* sender,FXSelector sel,void* ptr){
    FXScrollArea::onUngrabbed(sender,sel,ptr);
    flags&=~(FLAG_DODRAG|FLAG_TRYDRAG|FLAG_PRESSED|FLAG_CHANGED|FLAG_SCROLLING);
    flags|=FLAG_UPDATE;
    stopAutoScroll();
    return 1;
  }


  // Command message
  long FXCheckList::onCommand(FXObject*,FXSelector,void* ptr){
    return target && target->handle(this,MKUINT(message,SEL_COMMAND),ptr);
  }


  // Clicked in list
  long FXCheckList::onClicked(FXObject*,FXSelector,void* ptr){
    return target && target->handle(this,MKUINT(message,SEL_CLICKED),ptr);
  }


  // Double clicked in list; ptr may or may not point to an item
  long FXCheckList::onDoubleClicked(FXObject*,FXSelector,void* ptr){
    return target && target->handle(this,MKUINT(message,SEL_DOUBLECLICKED),ptr);
  }


  // Triple clicked in list; ptr may or may not point to an item
  long FXCheckList::onTripleClicked(FXObject*,FXSelector,void* ptr){
    return target && target->handle(this,MKUINT(message,SEL_TRIPLECLICKED),ptr);
  }


  // Extend selection
  FXbool FXCheckList::extendSelection(FXint index,FXbool notify){
    register FXbool changes=FALSE;
    FXint i1,i2,i3,i;
    if(0<=index && 0<=anchor && 0<=extent){

      // Find segments
      i1=index;
      if(anchor<i1){i2=i1;i1=anchor;}
      else{i2=anchor;}
      if(extent<i1){i3=i2;i2=i1;i1=extent;}
      else if(extent<i2){i3=i2;i2=extent;}
      else{i3=extent;}

      // First segment
      for(i=i1; i<i2; i++){

        // item===extent---anchor
        // item===anchor---extent
        if(i1==index){
          if(!items[i]->isSelected()){
            items[i]->setSelected(TRUE);
            updateItem(i);
            changes=TRUE;
            if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*)i);}
          }
        }

        // extent===anchor---item
        // extent===item-----anchor
        else if(i1==extent){
          if(items[i]->isSelected()){
            items[i]->setSelected(FALSE);
            updateItem(i);
            changes=TRUE;
            if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*)i);}
          }
        }
      }

      // Second segment
      for(i=i2+1; i<=i3; i++){

        // extent---anchor===item
        // anchor---extent===item
        if(i3==index){
          if(!items[i]->isSelected()){
            items[i]->setSelected(TRUE);
            updateItem(i);
            changes=TRUE;
            if(notify && target){target->handle(this,MKUINT(message,SEL_SELECTED),(void*)i);}
          }
        }

        // item-----anchor===extent
        // anchor---item=====extent
        else if(i3==extent){
          if(items[i]->isSelected()){
            items[i]->setSelected(FALSE);
            updateItem(i);
            changes=TRUE;
            if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*)i);}
          }
        }
      }
      extent=index;
    }
    return changes;
  }


  // Kill selection
  FXbool FXCheckList::killSelection(FXbool notify){
    register FXbool changes=FALSE;
    register FXint i;
    for(i=0; i<nitems; i++){
      if(items[i]->isSelected()){
        items[i]->setSelected(FALSE);
        updateItem(i);
        changes=TRUE;
        if(notify && target){target->handle(this,MKUINT(message,SEL_DESELECTED),(void*)i);}
      }
    }
    return changes;
  }


  // Sort items in ascending order
  FXint FXCheckList::ascending(const FXCheckListItem* a,const FXCheckListItem* b){
    return compare(a->label,b->label);
  }


  // Sort items in descending order
  FXint FXCheckList::descending(const FXCheckListItem* a,const FXCheckListItem* b){
    return compare(b->label,a->label);
  }


  // Sort the items based on the sort function
  void FXCheckList::sortItems(){
    register FXCheckListItem *v;
    register FXint i,j,h;
    register FXbool exch=FALSE;
    if(sortfunc){
      for(h=1; h<=nitems/9; h=3*h+1);
      if(hasFocus() && 0<=current) items[current]->setFocus(FALSE);
      for(; h>0; h/=3){
        for(i=h+1;i<=nitems;i++){
          v=items[i-1];
          j=i;
          while(j>h && sortfunc(items[j-h-1],v)>0){
            items[j-1]=items[j-h-1];
            exch=TRUE;
            j-=h;
          }
          items[j-1]=v;
        }
      }
      if(hasFocus() && 0<=current) items[current]->setFocus(TRUE);
      if(exch) update();
    }
  }


  // Set current item
  void FXCheckList::setCurrentItem(FXint index,FXbool notify){
    if(index<-1 || nitems<=index){ fxerror("%s::setCurrentItem: index out of range.\n",getClassName()); }
    if(index!=current){

      // Deactivate old item
      if(0<=current){

        // No visible change if it doen't have the focus
        if(hasFocus()){
          items[current]->setFocus(FALSE);
          updateItem(current);
        }
      }

      current=index;

      // Activate new item
      if(0<=current){

        // No visible change if it doen't have the focus
        if(hasFocus()){
          items[current]->setFocus(TRUE);
          updateItem(current);
        }
      }

      // Notify item change
      if(notify && target){target->handle(this,MKUINT(message,SEL_CHANGED),(void*)current);}
    }

    // In browse select mode, select this item
    if((options&SELECT_MASK)==CHECKLIST_BROWSESELECT && 0<=current && items[current]->isEnabled()){
      selectItem(current,notify);
    }
  }


  // Set anchor item
  void FXCheckList::setAnchorItem(FXint index){
    if(index<-1 || nitems<=index){ fxerror("%s::setAnchorItem: index out of range.\n",getClassName()); }
    anchor=index;
    extent=index;
  }


  // Create custom item
  FXCheckListItem *FXCheckList::createItem(const FXString& text,FXIcon* icon,void* ptr){
    return new FXCheckListItem(text,icon,ptr);
  }


  // Retrieve item
  FXCheckListItem *FXCheckList::retrieveItem(FXint index) const {
    if(index<0 || nitems<=index){ fxerror("%s::retrieveItem: index out of range.\n",getClassName()); }
    return items[index];
  }


  // Replace item with another
  FXint FXCheckList::replaceItem(FXint index,FXCheckListItem* item,FXbool notify){

    // Must have item
    if(!item){ fxerror("%s::replaceItem: item is NULL.\n",getClassName()); }

    // Must be in range
    if(index<0 || nitems<=index) { 
      fxerror("%s::replaceItem: index out of range.\n",getClassName()); 
    }

    // Notify item will be replaced
    if (notify && target) {
      target->handle(this,MKUINT(message,SEL_REPLACED),(void*)index);
    }

    // Copy the state over
    item->state=items[index]->state;

    // Delete old
    delete items[index];

    // Add new
    items[index]=item;

    // Redo layout
    recalc();
    return index;
  }


  // Replace item with another
  FXint FXCheckList::replaceItem(FXint index,const FXString& text,FXIcon *icon,void* ptr,FXbool notify){
    return replaceItem(index,createItem(text,icon,ptr),notify);
  }


  // Insert item
  FXint FXCheckList::insertItem(FXint index,FXCheckListItem* item,FXbool notify){
    register FXint old=current;

    // Must have item
    if(!item) { fxerror("%s::insertItem: item is NULL.\n",getClassName()); }

    // Must be in range
    if(index<0 || nitems<index) { 
      fxerror("%s::insertItem: index out of range.\n",getClassName()); 
    }

    // Add item to list
    FXRESIZE(&items,FXCheckListItem*,nitems+1);
    memmove(&items[index+1],&items[index],
            sizeof(FXCheckListItem*)*(nitems-index));
    items[index]=item;
    nitems++;

    // Adjust indices
    if(anchor>=index) anchor++;
    if(extent>=index) extent++;
    if(current>=index) current++;
    if(current<0 && nitems==1) current=0;

    // Notify item has been inserted
    if(notify && target) {
      target->handle(this,MKUINT(message,SEL_INSERTED),(void*)index);
    }

    // Current item may have changed
    if(old!=current){
      if(notify && target) {
        target->handle(this,MKUINT(message,SEL_CHANGED),(void*) current);
      }
    }

    // Was new item
    if(0<=current && current==index){
      if(hasFocus()){
        items[current]->setFocus(TRUE);
      }
      if ((options&SELECT_MASK)==CHECKLIST_BROWSESELECT && 
          items[current]->isEnabled()) {
        selectItem(current,notify);
      }
    }

    // Redo layout
    recalc();
    return index;
  }


  // Insert item
  FXint FXCheckList::insertItem(FXint index,const FXString& text,FXIcon *icon,
                                void* ptr,FXbool notify){
    return insertItem(index,createItem(text,icon,ptr),notify);
  }


  // Append item
  FXint FXCheckList::appendItem(FXCheckListItem* item,FXbool notify){
    return insertItem(nitems,item,notify);
  }


  // Append item
  FXint FXCheckList::appendItem(const FXString& text,FXIcon *icon,void* ptr,
                                FXbool notify){
    return insertItem(nitems,createItem(text,icon,ptr),notify);
  }


  // Prepend item
  FXint FXCheckList::prependItem(FXCheckListItem* item,FXbool notify){
    return insertItem(0,item,notify);
  }

  // Prepend item
  FXint FXCheckList::prependItem(const FXString& text,FXIcon *icon,void* ptr,
                                 FXbool notify){
    return insertItem(0,createItem(text,icon,ptr),notify);
  }


  // Remove node from list
  void FXCheckList::removeItem(FXint index,FXbool notify){
    register FXint old=current;

    // Must be in range
    if(index<0 || nitems<=index) { 
      fxerror("%s::removeItem: index out of range.\n",getClassName()); 
    }

    // Notify item will be deleted
    if(notify && target) {
      target->handle(this,MKUINT(message,SEL_DELETED),(void*)index);
    }

    // Remove item from list
    nitems--;
    delete items[index];
    memmove(&items[index],&items[index+1],
            sizeof(FXCheckListItem*)*(nitems-index));

    // Adjust indices
    if(anchor>index || anchor>=nitems)  anchor--;
    if(extent>index || extent>=nitems)  extent--;
    if(current>index || current>=nitems) current--;

    // Current item has changed
    if(old!=current){
      if(notify && target) {
        target->handle(this,MKUINT(message,SEL_CHANGED),(void*)current);
      }
    }

    // Deleted current item
    if(0<=current && index==old){
      if(hasFocus()){
        items[current]->setFocus(TRUE);
      }
      if((options&SELECT_MASK)==CHECKLIST_BROWSESELECT && 
         items[current]->isEnabled())      {
        selectItem(current,notify);
      }
    }

    // Redo layout
    recalc();
  }


  // Remove all items
  void FXCheckList::clearItems(FXbool notify){
    register FXint old=current;

    // Delete items
    for(FXint index=0; index<nitems; index++){
      if(notify && target) {
        target->handle(this,MKUINT(message,SEL_DELETED),(void*)index);
      }
      delete items[index];
    }

    // Free array
    FXFREE(&items);
    nitems=0;

    // Adjust indices
    current=-1;
    anchor=-1;
    extent=-1;

    // Current item has changed
    if(old!=current){
      if (notify && target) {
        target->handle(this,MKUINT(message,SEL_CHANGED),(void*)-1);
      }
    }

    // Redo layout
    recalc();
  }


  typedef FXint (*FXCompareFunc)(const FXString&,const FXString &,FXint);


  // Get item by name
  FXint FXCheckList::findItem(const FXString& text,FXint start,
                              FXuint flags) const {
    register FXCompareFunc comparefunc;
    register FXint index,len;
    if(0<nitems){
      comparefunc=(flags&SEARCH_IGNORECASE) ? (FXCompareFunc)comparecase : 
        (FXCompareFunc)compare;
      len=(flags&SEARCH_PREFIX)?text.length():2147483647;
      if(!(flags&SEARCH_BACKWARD)){
        if(start<0) start=0;
        for(index=start; index<nitems; index++){
          if((*comparefunc)(items[index]->label,text,len)==0) return index;
        }
        if(!(flags&SEARCH_WRAP)) return -1;
        for(index=0; index<start; index++){
          if((*comparefunc)(items[index]->label,text,len)==0) return index;
        }
      }
      else{
        if(start<0) start=nitems-1;
        for(index=start; 0<=index; index--){
          if((*comparefunc)(items[index]->label,text,len)==0) return index;
        }
        if(!(flags&SEARCH_WRAP)) return -1;
        for(index=nitems-1; start<index; index--){
          if((*comparefunc)(items[index]->label,text,len)==0) return index;
        }
      }
    }
    return -1;
  }


  // Change the font
  void FXCheckList::setFont(FXFont* fnt){
    if(!fnt){ fxerror("%s::setFont: NULL font specified.\n",getClassName()); }
    if(font!=fnt){
      font=fnt;
      recalc();
      update();
    }
  }


  // Set text color
  void FXCheckList::setTextColor(FXColor clr){
    textColor=clr;
    update();
  }


  // Set select background color
  void FXCheckList::setSelBackColor(FXColor clr){
    selbackColor=clr;
    update();
  }


  // Set selected text color
  void FXCheckList::setSelTextColor(FXColor clr){
    seltextColor=clr;
    update();
  }


  // Change list style
  void FXCheckList::setListStyle(FXuint style){
    options=(options&~LIST_MASK) | (style&LIST_MASK);
  }


  // Get list style
  FXuint FXCheckList::getListStyle() const {
    return (options&LIST_MASK);
  }


  // Change help text
  void FXCheckList::setHelpText(const FXString& text){
    help=text;
  }


  // Save data
  void FXCheckList::save(FXStream& store) const {
    register FXint i;
    FXScrollArea::save(store);
    store << nitems;
    for(i=0; i<nitems; i++){store<<items[i];}
    store << anchor;
    store << current;
    store << extent;
    store << textColor;
    store << selbackColor;
    store << seltextColor;
    store << itemWidth;
    store << itemHeight;
    store << visible;
    store << font;
    store << help;
  }


  // Load data
  void FXCheckList::load(FXStream& store){
    register FXint i;
    FXScrollArea::load(store);
    store >> nitems;
    FXRESIZE(&items,FXCheckListItem*,nitems);
    for(i=0; i<nitems; i++){store>>items[i];}
    store >> anchor;
    store >> current;
    store >> extent;
    store >> textColor;
    store >> selbackColor;
    store >> seltextColor;
    store >> itemWidth;
    store >> itemHeight;
    store >> visible;
    store >> font;
    store >> help;
  }


  // Clean up
  FXCheckList::~FXCheckList(){
#if FOX_MAJOR >=1 && FOX_MINOR < 4
    if(timer) {getApp()->removeTimeout(timer);}
    if(lookuptimer) {getApp()->removeTimeout(lookuptimer);}
    timer=(FXTimer*)-1;
    lookuptimer=(FXTimer*)-1;
#elif FOX_MAJOR >=1 && FOX_MINOR >= 4
    getApp()->removeTimeout(this, ID_TIPTIMER);
    getApp()->removeTimeout(this, ID_LOOKUPTIMER);
#endif
    clearItems(FALSE);
    items=(FXCheckListItem**)-1;
    font=(FXFont*)-1;
  }


}
