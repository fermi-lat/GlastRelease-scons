
#include "ConnectionDialog.h"
#include "fxkeys.h"

// Message Map ConnectionDialog class
FXDEFMAP(ConnectionDialog) ConnectionDialogMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,   ConnectionDialog::ID_LIST,          ConnectionDialog::onSelectProfile),
  FXMAPFUNC(SEL_KEYPRESS,  ConnectionDialog::ID_LIST,          ConnectionDialog::onKeyPress),
  FXMAPFUNC(SEL_COMMAND,   ConnectionDialog::ID_NEW,           ConnectionDialog::onNewProfile), 
  FXMAPFUNC(SEL_COMMAND,   ConnectionDialog::ID_SAVE,          ConnectionDialog::onSaveProfile),
  FXMAPFUNC(SEL_COMMAND,   ConnectionDialog::ID_DELETE,        ConnectionDialog::onDeleteProfile),
  FXMAPFUNC(SEL_KEYPRESS,  ConnectionDialog::ID_DBNAME,        ConnectionDialog::onKeyPress),
  FXMAPFUNC(SEL_KEYPRESS,  ConnectionDialog::ID_HOST,          ConnectionDialog::onKeyPress),
  FXMAPFUNC(SEL_KEYPRESS,  ConnectionDialog::ID_USER,          ConnectionDialog::onKeyPress),
  FXMAPFUNC(SEL_KEYPRESS,  ConnectionDialog::ID_PASSWORD,      ConnectionDialog::onKeyPress)
  };

// Object implementation
FXIMPLEMENT(ConnectionDialog,FXDialogBox,ConnectionDialogMap,ARRAYNUMBER(ConnectionDialogMap))


ConnectionDialog::ConnectionDialog(FXWindow *owner):
  FXDialogBox(owner, "Open connection", DECOR_TITLE|DECOR_BORDER, 0, 0, 0, 0, 0,
      0, 0, 0, 4, 4)
{
  FXHorizontalFrame *uiContent = new FXHorizontalFrame(this, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  
  // Left part - profiles list
  FXHorizontalFrame *uiListFrame = new FXHorizontalFrame(uiContent, FRAME_THICK|
      LAYOUT_FIX_WIDTH|LAYOUT_FILL_Y,0,0,100);
  m_uiList = new FXList(uiListFrame, this, ID_LIST, FRAME_SUNKEN|FRAME_SUNKEN|
      LAYOUT_FILL_X|LAYOUT_FILL_Y|LIST_SINGLESELECT);
  // Right panel - connection settings
  FXVerticalFrame *uiRpanel = new FXVerticalFrame(uiContent, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  // Right panel - buttons
  FXHorizontalFrame *uiProfilebox = new FXHorizontalFrame(uiRpanel, LAYOUT_TOP|
      LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  FXButton *newBtn = new FXButton(uiProfilebox, "&New", NULL, this, ID_NEW, FRAME_RAISED|
      FRAME_THICK, 0, 0, 0, 0, 20, 20);
  FXButton *saveBtn = new FXButton(uiProfilebox, "&Save", NULL, this, ID_SAVE,
      FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
  FXButton *deleteBtn = new FXButton(uiProfilebox, "&Delete", NULL, this, ID_DELETE,
      FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
  // Right panel - settings
  FXMatrix *uiMatrix = new FXMatrix(uiRpanel, 2, LAYOUT_FILL_X|LAYOUT_FILL_Y|
      FRAME_RAISED|MATRIX_BY_COLUMNS);
  // Profile name
  new FXLabel(uiMatrix, "profile name");
  m_uiProfile = new FXTextField(uiMatrix, 30, NULL, 0, FRAME_SUNKEN|FRAME_THICK,
      0, 0, 0, 0);
  // Profile description
  new FXLabel(uiMatrix, "profile description");
  m_uiDescription = new FXTextField(uiMatrix, 30, NULL, 0, FRAME_SUNKEN|FRAME_THICK, 
      0, 0, 0, 0);
  // Database name
  new FXLabel(uiMatrix, "database name");
  m_uiDbname = new FXTextField(uiMatrix, 30, this, ID_DBNAME, FRAME_SUNKEN|FRAME_THICK, 
      0, 0, 0, 0);
  // Host and port
  new FXLabel(uiMatrix, "host:port");
  m_uiHost = new FXTextField(uiMatrix, 30, this, ID_HOST, FRAME_SUNKEN|FRAME_THICK|
      LAYOUT_FILL_ROW, 0, 0, 0, 0);
  // Username
  new FXLabel(uiMatrix, "username");
  m_uiUser = new FXTextField(uiMatrix, 30, NULL, ID_USER, FRAME_SUNKEN|FRAME_THICK, 
      0, 0, 0, 0);
  // Password
  new FXLabel(uiMatrix, "password");
  m_uiPass = new FXTextField(uiMatrix, 30, NULL, ID_PASSWORD, FRAME_SUNKEN|
      FRAME_THICK|TEXTFIELD_PASSWD, 0, 0, 0, 0);
  // Optional password saving
  new FXLabel(uiMatrix, "save password");
  m_uiSavepass = new FXCheckButton(uiMatrix, "", NULL, 0, 0, 0, 0, 0, 0);
  // Bottom part
  FXHorizontalFrame *uiClosebox = new FXHorizontalFrame(uiRpanel, LAYOUT_BOTTOM
      |LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  new FXButton(uiClosebox, "&OK", NULL, this, FXDialogBox::ID_ACCEPT,
      BUTTON_INITIAL|LAYOUT_RIGHT|FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
  new FXButton(uiClosebox, "&Cancel", NULL, this, FXDialogBox::ID_CANCEL,
      LAYOUT_RIGHT|FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
  // Read saved profiles
  readProfiles();
  // Now there is a list of profiles with no one selected  
   
}

// Dialog settings for each execution
FXuint ConnectionDialog::execute(FXuint opts)
{
    if (m_uiList->getCurrentItem() >= 0)
      m_uiList->setFocus();
    else
      m_uiDbname->setFocus();
      
    if (m_uiSavepass->getCheck())
      m_uiPass->setText("");
    FXDialogBox::execute(opts);
}

// Quick finish if the return key was pressed, pass other keys
long ConnectionDialog::onKeyPress(FXObject *sender,FXSelector sel, void* ptr)
{  
  FXEvent *event = (FXEvent*) ptr;
  if (event->code == KEY_Return || event->code == KEY_KP_Enter)
    return 0;
  //handle(this, MKUINT(ID_ACCEPT, SEL_COMMAND), NULL);
  return FXDialogBox::onKeyPress(sender, sel, ptr);
}


std::vector<FXString> ConnectionDialog::splitString(const FXString& str, const FXchar* delim)
{
  std::vector<FXString> elems; 
  FXint start, end;
  start = 0;
  end = str.length();
  while (start < end)
    {
      int pos = str.find(delim,1,start);
      if (pos == -1)
        {
          pos = end;
        }
      FXString elem = (pos > start) ? str.mid(start, pos-start).trim(): "";
      elems.push_back(elem.mid(0, elem.length()));
      start = pos+1;
    }
  return elems;
}

// Read all profiles and populate list
void ConnectionDialog::readProfiles()
{
  FXString record = getApp()->reg().readStringEntry("PROFILES", "VALUE", "");
  if (record.empty())
    return;
  
  std::vector<FXString> lines = splitString(record, "\n"); 

  for (int i = 0; i < lines.size(); i++)
    {
      FXString *line = new FXString(lines[i]);
      std::vector<FXString> settings = splitString(lines[i],"|");
      m_uiList->appendItem(settings[0], NULL, line);
    }
}


// Save profiles
void ConnectionDialog::saveProfiles()
{  
  FXString settings = "";
  for (int i = 0; i < m_uiList->getNumItems(); i++)
    {
      settings += ((FXString*)m_uiList->getItemData(i))->text();
      settings += "\n";
    }
  getApp()->reg().writeStringEntry("PROFILES", "VALUE", settings.text());
}


// Change settings based on profile
long ConnectionDialog::onSelectProfile(FXObject *,FXSelector, void*)
{
  m_uiList->selectItem(m_uiList->getCurrentItem()); // ensure current == selected
  FXString *settingsStr = (FXString*)m_uiList->getItemData(m_uiList->getCurrentItem());
  std::vector<FXString> settings = splitString(*settingsStr, "|");
  m_uiProfile->setText(settings[0]);
  m_uiDescription->setText(settings[1]);
  m_uiDbname->setText(settings[2]);
  m_uiHost->setText(settings[3]);
  m_uiUser->setText(settings[4]);
  m_uiPass->setText(settings[5]);
  m_uiSavepass->setCheck(settings[6] == "Y");
  return 1;
}


// Save settings to include new profile or just settings to connect
long ConnectionDialog::onNewProfile(FXObject *,FXSelector, void*)
{
  m_uiList->killSelection(false);
  m_uiProfile->setText("");
  m_uiDescription->setText("");
  m_uiDbname->setText("");
  m_uiHost->setText("");
  m_uiUser->setText("");
  m_uiPass->setText("");
  m_uiSavepass->setCheck(false);
  m_uiDbname->setFocus();
  return 1;
}


// Save profile settings
long ConnectionDialog::onSaveProfile(FXObject *,FXSelector, void*)
{
  if (m_uiProfile->getText().empty())
    {
      FXMessageBox::question(this, MBOX_OK, "Save profile", "Missing name for this profile.");
      m_uiProfile->setFocus();
      return 1;
    }
  FXString *settings = new FXString();
  *settings += m_uiProfile->getText() + "|";
  *settings += m_uiDescription->getText() + "|";
  *settings += m_uiDbname->getText() + "|";
  *settings += m_uiHost->getText() + "|";
  *settings += m_uiUser->getText() + "|";
  *settings += ((m_uiSavepass->getCheck()) ? m_uiPass->getText():"") + "|";
  *settings += ((m_uiSavepass->getCheck()) ? "Y":"N");
  // Replace item if selected, append otherwise
  if ((m_uiList->getCurrentItem() >=0) && (m_uiList->isItemSelected(m_uiList->getCurrentItem())))
    {
      FXListItem *item = m_uiList->getItem(m_uiList->getCurrentItem());
      item->setText(m_uiProfile->getText());
      item->setData(settings);
    }
  else
      m_uiList->appendItem(m_uiProfile->getText(), NULL, settings);
  
  saveProfiles();
  return 1;
}


// Delete selected profile
long ConnectionDialog::onDeleteProfile(FXObject *,FXSelector, void*)
{
  // Check if there is at least one item in the list and it is selected
  if ((m_uiList->getCurrentItem() >=0) && (m_uiList->isItemSelected(m_uiList->getCurrentItem())))
    {
      FXString message = "Delete profile" + m_uiList->getItemText(m_uiList->getCurrentItem()) + "?";
      
      if (MBOX_CLICKED_YES == FXMessageBox::question(getOwner(), MBOX_YES_NO, "Delete profile", message.text()))
        {
          m_uiList->removeItem(m_uiList->getCurrentItem());
          saveProfiles();
        }
    }
  else
    FXMessageBox::question(this, MBOX_OK, "Delete profile", "No profile selected.");

  return 1;
}

// Return array with connection attributes
std::vector<FXString> ConnectionDialog::getConnectionData()
{
  std::vector<FXString> data, host;
  data.push_back(m_uiDbname->getText());
  host = splitString(m_uiHost->getText(), ":");
  data.insert(data.end(), host.begin(), host.end());
  data.push_back(m_uiUser->getText());
  data.push_back(m_uiPass->getText());
  return data;
}
