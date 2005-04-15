/** @file InsertDialog.cxx
* @brief InsertDialog is a non-modal dialog to insert new row in the calib
* database.
* @author Riccardo Giannitrapani
* @author Marco Frailis
* $Header$
*/

#include "InsertDialog.h"
#include "LogText.h"
#include "Icons.h"
#include "fxkeys.h"
#include "facilities/Util.h"
#include "rdbModel/Tables/Datatype.h"
#include "ColWidgetFactory.h"
#include <iostream>
#include <cstdio>


// Message Map ConnectionDialog class
FXDEFMAP(InsertDialog) InsertDialogMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,  InsertDialog::ID_GO,        InsertDialog::onGoPress),
  };

// Object implementation
FXIMPLEMENT(InsertDialog,FXDialogBox,InsertDialogMap,ARRAYNUMBER(InsertDialogMap))


InsertDialog::InsertDialog(FXApp *owner):
  FXDialogBox(owner, "Insert",DECOR_TITLE|DECOR_BORDER|DECOR_RESIZE)
{ 
  m_matrix = 0;
 
  m_uiRpanel = new FXVerticalFrame(this, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  // Bottom part
  FXHorizontalFrame *uiClosebox = new FXHorizontalFrame(m_uiRpanel,LAYOUT_BOTTOM|LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  new FXButton(uiClosebox, "&Cancel", NULL, this, FXDialogBox::ID_CANCEL,
      LAYOUT_RIGHT|FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
  new FXButton(uiClosebox, "&Send", NULL, this, ID_GO,
               BUTTON_INITIAL|LAYOUT_RIGHT|FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
               
  m_connection = 0;
  m_factory = new ColWidgetFactory();

  m_lastTblName = "";
  // -1 means no row has yet been inserted in this session
  m_lastRow = -1;

  // Let start in insert mode
  m_insertMode = 1;

  multi->hide();  
}


void InsertDialog::fillWithRowByKey(std::string primKeyVal)
{
  unsigned int i;
  
  std::vector<std::string> colNames;
  std::vector<std::string> colValues;

  std::string name;
   
  for(i=0;i<m_widgets.size();i++)
  {
    ColWidget* temp = m_widgets[i]; 
    
    name = temp->getColumn()->getName();
    facilities::Util::trimTrailing(&name);
    colNames.push_back(name); 
  }

  // Set up WHERE clause, always the same
  rdbModel::Assertion::Operator* serEquals = 
    new rdbModel::Assertion::Operator(rdbModel::OPTYPEequal, m_primKey,
                                      primKeyVal, false, true);
  rdbModel::Assertion* whereSer = new rdbModel::Assertion(rdbModel::Assertion::WHENwhere, serEquals);


  if(m_connection)
    m_result = m_connection->select(m_tableName, colNames, colNames, whereSer);


  m_result->getRow(colValues, 0, 1);

  for(i=0;i<m_widgets.size();i++)
  {
    ColWidget* temp = m_widgets[i]; 
    temp->setValue(colValues[i].c_str());
  }
  m_selRow = primKeyVal;
}

void InsertDialog::fillStickyWithLastRow()
{
  std::string primKeyVal;
  primKeyVal = FXStringVal(m_lastRow).text();

  unsigned int i;
  
  std::vector<std::string> colNames;
  std::vector<std::string> colValues;

  std::string name;
   
  for(i=0;i<m_widgets.size();i++)
  {
    ColWidget* temp = m_widgets[i]; 
    
    name = temp->getColumn()->getName();
    facilities::Util::trimTrailing(&name);
    colNames.push_back(name); 
  }

  // Set up WHERE clause, always the same
  rdbModel::Assertion::Operator* serEquals = 
    new rdbModel::Assertion::Operator(rdbModel::OPTYPEequal, m_primKey,
                                      primKeyVal, false, true);
  rdbModel::Assertion* whereSer = new rdbModel::Assertion(rdbModel::Assertion::WHENwhere, serEquals);


  if(m_connection)
    m_result = m_connection->select(m_tableName, colNames, colNames, whereSer);


  m_result->getRow(colValues, 0, 1);

  for(i=0;i<m_widgets.size();i++)
  {
    ColWidget* temp = m_widgets[i]; 
    if (temp->getColumn()->stickyInsert())
      temp->setValue(colValues[i].c_str());
  }
  m_selRow = primKeyVal;
  
}


// Get the last row from the database and fill the dialog form
void InsertDialog::fillWithLastRow()
{
  std::string primKeyVal;
  primKeyVal = FXStringVal(m_lastRow).text();
  fillWithRowByKey(primKeyVal);
  m_selRow = primKeyVal;
}

// Try to insert the new row
long InsertDialog::onGoPress(FXObject *,FXSelector, void*)
{     
  unsigned int i;
  
  std::vector<std::string> colNames;
  std::vector<std::string> values;
  std::vector<std::string> nullValues;
  std::string name;
  std::string value;
  
#ifdef WIN32
  const char* usrName = "USERNAME";
#else
  const char* usrName = "USER";
#endif
   
  for(i=0;i<m_widgets.size();i++)
  {
    ColWidget* temp = m_widgets[i]; 
    
    name = temp->getColumn()->getName();
    facilities::Util::trimTrailing(&name);
    value = temp->getValue();
    facilities::Util::trimTrailing(&value);
    
    if (value == "")
      nullValues.push_back(name);
    else
    {
      colNames.push_back(name); 
      values.push_back(value);
    }      
  }
   
  for (i = 0; i < m_fromService.size(); i++)
    if (m_fromService[i]->getContentsType() == rdbModel::Column::CONTENTSserviceName)
      {
        colNames.push_back(m_fromService[i]->getName());
        values.push_back("rdbGui");
      }
    else if (m_fromService[i]->getContentsType() == rdbModel::Column::CONTENTSusername)
      {
        colNames.push_back(m_fromService[i]->getName());
        values.push_back(::getenv(usrName));        
      }
      
  if (m_connection)
  {
    if(m_insertMode) //If in insert mode
      {
        m_connection->insertRow(m_tableName, colNames, values, &m_lastRow, &nullValues);
        m_lastTblName = m_tableName;
      }  
    else //Otherwise we are in the Update Last Row mode
    {
      std::string serialStr;
      rdbModel::Assertion::Operator* serEquals = 
       new rdbModel::Assertion::Operator(rdbModel::OPTYPEequal, m_primKey,
                                      m_selRow, false, true);
      rdbModel::Assertion* whereSer = new rdbModel::Assertion(rdbModel::Assertion::WHENwhere, serEquals);

      m_connection->update(m_tableName, colNames, values, whereSer, &nullValues);  
    }
  }
    
      
    
  m_uiLog->update();
  this->handle(this, MKUINT(ID_ACCEPT, SEL_COMMAND),NULL);
  return 1;
}


/// Visitors to build the insert widgets  
rdbModel::Visitor::VisitorState InsertDialog::visitRdb(rdbModel::Rdb *)
{
  if (m_matrix != 0)
  {
    m_matrix->destroy();
    delete m_matrix;
  }
  
  m_matrix = new FXMatrix(m_uiRpanel,2,MATRIX_BY_COLUMNS|LAYOUT_FILL_X,0,0,0,0,0,0,0,0,0,4);
  m_matrix->create();

  // Let's delete the previous widgets
  unsigned int i;
  for(i=0;i<m_widgets.size();i++)
  {
    delete m_widgets[i];
    m_widgets.clear();      
  }
  
  m_fromService.clear();
    
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitTable(rdbModel::Table *tab)
{
  if (tab->getName()!=m_tableName)
    return rdbModel::Visitor::VBRANCHDONE;
  else
    return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitColumn(rdbModel::Column *column)
{
  rdbModel::Datatype* dt = column->getDatatype();
  rdbModel::Column::FROM source = column->getSourceType();
  
  if (column->isPrimaryKey())
    m_primKey = column->getName();
  
  std::string min;
  std::string max;

  if ((source == rdbModel::Column::FROMdefault) || 
      (source == rdbModel::Column::FROMendUser))
  {
    //    if (dt->getInterval(min, max))
    //      std::cout << column->getName() << " " << min << " " << max << std::endl;
  
    FXLabel* label = new FXLabel(m_matrix,column->getName().c_str());
    label->create();
    label->setTipText(column->getComment().c_str());
    
    ColWidget* colWidget;

    switch(dt->getType()){
    case rdbModel::Datatype::TYPEdatetime :
    {
      	// label->setIcon(new FXICOIcon(getApp(), calImg));	
      	colWidget = m_factory->createDateWidget(m_matrix, column);
        break;
    };
    case rdbModel::Datatype::TYPEtimestamp :
    {
      	//label->setIcon(new FXICOIcon(getApp(), clockImg));
        colWidget = m_factory->createStringWidget(m_matrix, column);
        break;
    };
    case rdbModel::Datatype::TYPEenum :
    {
      	colWidget = m_factory->createEnumWidget(m_matrix, column);
        break;
    };
    case rdbModel::Datatype::TYPEreal :
    case rdbModel::Datatype::TYPEdouble :
    {
        colWidget = m_factory->createRealWidget(m_matrix, column);
        break;
    };
    case rdbModel::Datatype::TYPEint :
    case rdbModel::Datatype::TYPEmediumint :
    case rdbModel::Datatype::TYPEsmallint :
    {
        colWidget = m_factory->createIntWidget(m_matrix, column);
        break;
    };
    case rdbModel::Datatype::TYPEvarchar :
    case rdbModel::Datatype::TYPEchar :
    {
      if (dt->getRestrict() == rdbModel::Datatype::RESTRICTenum)
        colWidget = m_factory->createEnumWidget(m_matrix, column);
      else
        colWidget = m_factory->createStringWidget(m_matrix, column);
      break;
    };
    case rdbModel::Datatype::TYPEnotFound :
    {
      colWidget = NULL;
      break;
    };
    } 

    if (!column->nullAllowed())
      label->setTextColor(FXRGB(255,0,0));

    if ((source == rdbModel::Column::FROMdefault)) 
      colWidget->setValue(column->getDefault());
      
    m_widgets.push_back(colWidget);
  } 
  
  if (source == rdbModel::Column::FROMprogram)
    m_fromService.push_back(column);
    
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitIndex(rdbModel::Index *)
{
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitAssertion(rdbModel::Assertion *)
{
  return rdbModel::Visitor::VCONTINUE;
}


