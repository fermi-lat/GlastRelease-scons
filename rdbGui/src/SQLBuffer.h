
#ifndef SQLBUFFER_H
#define SQLBUFFER_H

#include "fx.h"
#include "FXRex.h"
#include <map>
#include <vector>

class SQLBuffer : public FXText
{

 private:
  
  static const FXRex SQLCMD;

  /// Number of subgroups into the SQLCMD regular expression
  static const int m_regexGroups;

  std::map<FXString, unsigned char> KEYWCOLS;

  std::vector<FXString> KEYWORDS;

  FXRex m_shregexp;

  FXuint m_shsum;


 public:

  SQLBuffer(FXComposite *p,FXObject* tgt=NULL,FXSelector sel=0,FXuint opts=0,FXint x=0,FXint y=0,FXint w=0,FXint h=0);

  /// Get characters from the start of the current word up to the cursor
  FXString getWordBefore();

  /// Find boundaries of the SQL query in the buffer starting at the position
  FXint* getSQLBoundaries(FXint from = 0);

  /// Find boundaries of the SQL query under the cursor
  FXint* getCurrentSQL();

  /// Get text of the SQL query under the cursor
  FXString getCurrentSQLText(FXint *bound, FXbool atbegin = false);

  /// Select text of the SQL query under the cursor
  void selectCurrentSQL(FXbool atbegin = false);

  /// Read color settings from registry or set defaults
  std::vector<FXColor> readColors();

  /// Set the editor to styled mode and create styles from given colors
  void createStyles(std::vector<FXColor> &colors);

  /// Prepare the SQL buffer for syntax highlighting
  void initSyntaxHighlighting();

  /// Regexp through the visible part of the text buffer
  /// and set the appropriate highlighting style
  void highlightSyntax(FXbool complete = false);


};


#endif  // end of SQLBUFFER_H
