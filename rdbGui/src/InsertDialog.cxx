#include "InsertDialog.h"
#include "Icons.h"
#include "fxkeys.h"
#include "rdbModel/Tables/Datatype.h"
#include "ColWidgetFactory.h"
#include <iostream>

// Message Map ConnectionDialog class
FXDEFMAP(InsertDialog) InsertDialogMap[]={

  //__Message_Type_____________ID________________________Message_Handler_____
  FXMAPFUNC(SEL_COMMAND,  InsertDialog::ID_GO,        InsertDialog::onGoPress),
  };

// Object implementation
FXIMPLEMENT(InsertDialog,FXDialogBox,InsertDialogMap,ARRAYNUMBER(InsertDialogMap))


InsertDialog::InsertDialog(FXWindow *owner):
  FXDialogBox(owner, "Insert", DECOR_TITLE|DECOR_BORDER, 0, 0, 0, 0, 0,
      0, 0, 0, 4, 4)
{
  FXHorizontalFrame *uiContent = new FXHorizontalFrame(this, LAYOUT_FILL_X|LAYOUT_FILL_Y);
 
  m_matrix = 0;
  
  m_uiRpanel = new FXVerticalFrame(uiContent, LAYOUT_FILL_X|LAYOUT_FILL_Y);
  // Bottom part
  FXHorizontalFrame *uiClosebox = new FXHorizontalFrame(m_uiRpanel, LAYOUT_BOTTOM
      |LAYOUT_FILL_X|PACK_UNIFORM_WIDTH);
  new FXButton(uiClosebox, "&Cancel", NULL, this, FXDialogBox::ID_CANCEL,
      LAYOUT_RIGHT|FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);
  new FXButton(uiClosebox, "&Send", NULL, this, ID_GO,
               BUTTON_INITIAL|LAYOUT_RIGHT|FRAME_RAISED|FRAME_THICK, 0, 0, 0, 0, 20, 20);

  m_factory = new ColWidgetFactory();
}

// Quick finish if the return key was pressed, pass other keys
long InsertDialog::onGoPress(FXObject *sender,FXSelector sel, void* ptr)
{     
  unsigned int i;
  
  for(i=0;i<m_widgets.size();i++)
  {
    ColWidget* temp = m_widgets[i]; 
    std::cout << temp->getColumn()->getName() << " -> ";    
    std::cout << temp->getValue() << std::endl;
  }
    
  this->handle(this, MKUINT(ID_ACCEPT, SEL_COMMAND),NULL);
  return 1;
}


/// Visitors to build the insert widgets  
rdbModel::Visitor::VisitorState InsertDialog::visitRdb(rdbModel::Rdb *rdb)
{
  if (m_matrix != 0)
  {
    m_matrix->destroy();
    delete m_matrix;
  }
  
  m_matrix = new FXMatrix(m_uiRpanel,2,MATRIX_BY_COLUMNS,0,0,0,0,0,0,0,0,0,4);
  m_matrix->create();

  // Let's delete the previous widgets
  unsigned int i;
  for(i=0;i<m_widgets.size();i++)
  {
    delete m_widgets[i];
    m_widgets.clear();      
  }
    
  
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitTable(rdbModel::Table *table)
{
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitColumn(rdbModel::Column *column)
{
  rdbModel::Datatype* dt = column->getDatatype();
  rdbModel::Column::FROM source = column->getSourceType();

  std::string min;
  std::string max;

  if ((source == rdbModel::Column::FROMdefault) || 
      (source == rdbModel::Column::FROMendUser))
  {
    //    if (dt->getInterval(min, max))
    //      std::cout << column->getName() << " " << min << " " << max << std::endl;
  
    FXLabel* label = new FXLabel(m_matrix,column->getName().c_str());
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
        colWidget = m_factory->createStringWidget(m_matrix, column);
    };
    }
    
    if (source != rdbModel::Column::FROMdefault)
        label->setTextColor(FXRGB(255,0,0));
    else
        {
            colWidget->setValue(column->getDefault());
        }

    m_widgets.push_back(colWidget);
  }  
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitIndex(rdbModel::Index *index)
{
  return rdbModel::Visitor::VCONTINUE;
}

rdbModel::Visitor::VisitorState InsertDialog::visitAssertion(rdbModel::Assertion *assertion)
{
  return rdbModel::Visitor::VCONTINUE;
}


