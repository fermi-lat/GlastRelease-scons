
#ifndef CONNECTIONDIALOG_H
#define CONNECTIONDIALOG_H

#include "fx.h"
#include <vector>

class ConnectionDialog: public FXDialogBox
{
  FXDECLARE(ConnectionDialog)
 public:
 
   enum{
    ID_LIST=FXDialogBox::ID_LAST,
    ID_NEW,
    ID_SAVE,
    ID_DELETE,
    ID_DBNAME,
    ID_HOST,
    ID_USER,
    ID_PASSWORD
    };

 
  ConnectionDialog(FXWindow *owner);
  
  FXuint execute(FXuint);
  
  long onKeyPress(FXObject*,FXSelector,void*);
  
  void readProfiles();
  void saveProfiles(); 
  
  long onSelectProfile(FXObject*,FXSelector,void*);
  long onNewProfile(FXObject*,FXSelector,void*);
  long onSaveProfile(FXObject*,FXSelector,void*);
  long onDeleteProfile(FXObject*,FXSelector,void*);
  
  std::vector<FXString> getConnectionData();
  
 protected:
  ConnectionDialog(){}
  ConnectionDialog(const ConnectionDialog&);
  
 private:
 
  FXList *m_uiList;  // List of connections profiles
  FXTextField *m_uiProfile;
  FXTextField *m_uiDescription;
  FXTextField *m_uiDbname;
  FXTextField *m_uiHost;
  FXTextField *m_uiUser;
  FXTextField *m_uiPass;
  FXCheckButton *m_uiSavepass;
  
  std::vector<FXString> splitString(const FXString& str, const FXchar* delim);
};


#endif
