
#include "SQLBuffer.h"


// Regular expression for the complete SQL query (hopefully covers 99.9% of cases)
const FXRex SQLBuffer::SQLCMD("(--.*$|[^-;\'\"]+|-|\'(\\\\\'|[^\'])*(\'|\\z)|\"(\\\\\"|[^\"])*(\"|\\z))*(;|\\z)(--.*|[ ]+)*",
			      REX_CAPTURE);

// Number of subgroups into the SQLCMD regular expression
const int SQLBuffer::m_regexGroups = 7;

SQLBuffer::SQLBuffer(FXComposite *p,FXObject* tgt,FXSelector sel,FXuint opts,FXint x,FXint y,FXint w,FXint h):
  FXText(p, tgt, sel, opts, x, y, w, h) 
{ 

  // SQL keywords and their color styles
  KEYWCOLS["explain"] = 1, KEYWCOLS["select"]= 1, KEYWCOLS["distinct"] = 2;
  KEYWCOLS["from"] = 2, KEYWCOLS["where"] = 2, KEYWCOLS["order by"] = 2; 
  KEYWCOLS["ascending"] = 2, KEYWCOLS["descending"] = 2, KEYWCOLS["group by"] = 2;
  KEYWCOLS["having"] = 2, KEYWCOLS["as"] = 2, KEYWCOLS["limit"] = 2; 
  KEYWCOLS["offset"] = 2, KEYWCOLS["top"] = 2, KEYWCOLS["full"] = 2; 
  KEYWCOLS["left"] = 2, KEYWCOLS["right"] = 2, KEYWCOLS["natural"] = 2; 
  KEYWCOLS["outer"] = 2, KEYWCOLS["join"] = 2, KEYWCOLS["on"] = 2;
  KEYWCOLS["union"] = 1, KEYWCOLS["intersect"] = 1, KEYWCOLS["except"] = 1; 
  KEYWCOLS["minus"] = 1, KEYWCOLS["all"] = 1, KEYWCOLS["insert"] = 1; 
  KEYWCOLS["into"] = 2, KEYWCOLS["values"] = 2, KEYWCOLS["update"] = 1; 
  KEYWCOLS["set"] = 1, KEYWCOLS["delete"] = 1, KEYWCOLS["declare"] = 1; 
  KEYWCOLS["cursor"] = 2, KEYWCOLS["open"] = 1, KEYWCOLS["close"] = 1; 
  KEYWCOLS["fetch"] = 1, KEYWCOLS["move"] = 1,KEYWCOLS["create"] = 1;
  KEYWCOLS["replace"] = 1, KEYWCOLS["alter"] = 1, KEYWCOLS["drop"] = 1;
  KEYWCOLS["database"] = 2, KEYWCOLS["schema"] = 2, KEYWCOLS["table"] = 2;
  KEYWCOLS["view"] = 2, KEYWCOLS["unique"] = 2, KEYWCOLS["index"] = 2;
  KEYWCOLS["sequence"] = 2, KEYWCOLS["constraint"] = 2, KEYWCOLS["rule"] = 2;
  KEYWCOLS["aggregate"] = 2, KEYWCOLS["type"] = 2, KEYWCOLS["operator"] = 2;
  KEYWCOLS["user"] = 2, KEYWCOLS["group"] = 2, KEYWCOLS["password"] = 2;
  KEYWCOLS["trigger"] = 2, KEYWCOLS["before"] = 2, KEYWCOLS["after"] = 2;
  KEYWCOLS["primary"] = 2, KEYWCOLS["foreign"] = 2, KEYWCOLS["key"] = 2;
  KEYWCOLS["references"] = 2, KEYWCOLS["default"] = 2, KEYWCOLS["check"] = 2;
  KEYWCOLS["grant"] = 1, KEYWCOLS["revoke"] = 1, KEYWCOLS["start"] = 1;
  KEYWCOLS["transaction"] = 2, KEYWCOLS["isolation level"] = 2, KEYWCOLS["read committed"] = 2;
  KEYWCOLS["serializable"] = 2, KEYWCOLS["commit"] = 1, KEYWCOLS["rollback"] = 1;
  KEYWCOLS["abort"] = 1, KEYWCOLS["checkpoint"] = 1, KEYWCOLS["lock"] = 1;
  KEYWCOLS["vacuum"] = 1, KEYWCOLS["analyze"] = 2, KEYWCOLS["verbose"] = 2;
  KEYWCOLS["cluster"] = 1, KEYWCOLS["comment"] = 1, KEYWCOLS["show"] = 1;
  KEYWCOLS["truncate"] = 1, KEYWCOLS["is"] = 3, KEYWCOLS["in"] = 3;
  KEYWCOLS["between"] = 3, KEYWCOLS["null"] = 3, KEYWCOLS["exists"] = 3;
  KEYWCOLS["and"] = 3, KEYWCOLS["or"] = 3, KEYWCOLS["not"] = 3, KEYWCOLS["like"] = 3;
  KEYWCOLS["count"] = 3, KEYWCOLS["sum"] = 3, KEYWCOLS["max"] = 3;
  KEYWCOLS["min"] = 3, KEYWCOLS["avg"] = 3, KEYWCOLS["procedure"] = 2;
  KEYWCOLS["function"] = 2, KEYWCOLS["package"] = 2, KEYWCOLS["body"] = 2;
  KEYWCOLS["returns"] = 2, KEYWCOLS["language"] = 2, KEYWCOLS["with"] = 2;
  KEYWCOLS["using"] = 2, KEYWCOLS["alias"] = 1, KEYWCOLS["begin"] = 1; 
  KEYWCOLS["end"] = 1, KEYWCOLS["return"] = 1, KEYWCOLS["exit"] = 1;
  KEYWCOLS["raise"] = 1, KEYWCOLS["exception"] = 2, KEYWCOLS["case"] = 1;
  KEYWCOLS["when"] = 1, KEYWCOLS["if"] = 1, KEYWCOLS["then"] = 1;
  KEYWCOLS["elsif"] = 1, KEYWCOLS["else"] = 1, KEYWCOLS["for"] = 1;
  KEYWCOLS["loop"] = 1, KEYWCOLS["int"] = 3, KEYWCOLS["integer"] = 3;
  KEYWCOLS["serial"] = 3,KEYWCOLS["boolean"] = 3, KEYWCOLS["number"] = 3;
  KEYWCOLS["numeric"] = 3,KEYWCOLS["decimal"] = 3, KEYWCOLS["float"] = 3;
  KEYWCOLS["double"] = 3, KEYWCOLS["opaque"] = 3, KEYWCOLS["char"] = 3;
  KEYWCOLS["varchar"] = 3, KEYWCOLS["variable"] = 3, KEYWCOLS["character"] = 3;
  KEYWCOLS["time"] = 3, KEYWCOLS["date"] = 3, KEYWCOLS["timestamp"] = 3;
  KEYWCOLS["datetime"] = 3, KEYWCOLS["by"] = 2;

  std::map<FXString, unsigned char>::const_iterator i;

  for (i = KEYWCOLS.begin(); i != KEYWCOLS.end(); i++)
    KEYWORDS.push_back(i->first);
  
}
 
// Get characters from the start of the current word up to the cursor
FXString SQLBuffer::getWordBefore()
{
  FXint before = leftWord(cursorpos);
  FXchar *text=0;
  extractText(text, before, cursorpos - before);
  FXString textString = text;
  delete text;
  return textString.trim();
}


// Find boundaries of the SQL query in the buffer starting at the position
FXint* SQLBuffer::getSQLBoundaries(FXint from)
{
  FXint beg[m_regexGroups], end[m_regexGroups];  
  FXString text = getText().mid(from, length-from);
  SQLCMD.match(text, beg, end, REX_FORWARD, m_regexGroups);
  FXint* bound = new FXint[2];
  bound[0] = from + beg[0] + ((text[beg[0]] == 10) ? 1 : 0);
  bound[1] = from + end[0] + 1;
  return bound;
}


// Find boundaries of the SQL query under the cursor
FXint* SQLBuffer::getCurrentSQL()
{
  FXint *bound = getSQLBoundaries();
  while  (!((cursorpos >= bound[0]) && (cursorpos <= bound[1])))
    bound = getSQLBoundaries(bound[1]);
  return bound;
}


// Get text of the SQL query under the cursor
FXString SQLBuffer::getCurrentSQLText(FXint *bound, FXbool atbegin)
{
  bound = (atbegin) ? getSQLBoundaries(cursorpos) : getCurrentSQL();
  FXchar *text=0;
  extractText(text, bound[0], bound[1] - bound[0]);
  FXString out = text;
  delete text;
  return out;
}


// Select text of the SQL query under the cursor
void SQLBuffer::selectCurrentSQL(FXbool atbegin)
{
  FXint *bound = (atbegin) ? getSQLBoundaries(cursorpos) : getCurrentSQL();
  setSelection(bound[0], bound[1] - bound[0]);
  delete bound;
}


// Read color settings from registry or set defaults
std::vector<FXColor> SQLBuffer::readColors()
{
  std::vector<FXColor> colors;
  colors.push_back(getApp()->reg().readColorEntry("SETTINGS", "SHLiteralsColor", fxcolorfromname("DarkGreen")));
  colors.push_back(getApp()->reg().readColorEntry("SETTINGS", "SHCommentsColor", fxcolorfromname("DarkGoldenrod")));
  colors.push_back(getApp()->reg().readColorEntry("SETTINGS", "SHKeyword1Color", fxcolorfromname("DarkViolet")));
  colors.push_back(getApp()->reg().readColorEntry("SETTINGS", "SHKeyword2Color", fxcolorfromname("RoyalBlue")));
  colors.push_back(getApp()->reg().readColorEntry("SETTINGS", "SHKeyword3Color", fxcolorfromname("NavyBlue")));
  return colors;
}


// Set the editor to styled mode and create styles from given colors
void SQLBuffer::createStyles(std::vector<FXColor> &colors)
{
  // Create styles array
  FXHiliteStyle *styles = new FXHiliteStyle[colors.size()];
  
  unsigned int i;
  for (i = 0; i < colors.size(); i++)
    {
      styles[i].normalForeColor = colors[i];
      styles[i].normalBackColor = backColor;
      styles[i].selectForeColor = seltextColor;
      styles[i].selectBackColor = selbackColor;
      styles[i].hiliteForeColor = hilitetextColor;
      styles[i].hiliteBackColor = hilitebackColor;
      styles[i].activeBackColor = activebackColor;
      styles[i].style = 0;
    }
  // Enable the style buffer
  setStyled();
  // Set the styles ( a possible little memory leak)
  setHiliteStyles(styles);
}


// Prepare the SQL buffer for syntax highlighting
void SQLBuffer::initSyntaxHighlighting()
{
  // Construct the reqular expression for SQL literals,
  FXString regexpsrc = "^\' |\'$|\'(\\\\\'|[^\'])*(\'|\\z)|\"(\\\\\"|[^\"])*(\"|\\z)";
  // numbers,
  regexpsrc += "|\\b\\d+\\.?\\d*([eE][+-]?\\d+)?\\b";
  // comments,
  regexpsrc += "|--.*$";
  // and keywords 
  std::vector<FXString>::const_iterator i;
  for (i = KEYWORDS.begin(); i != KEYWORDS.end(); i++)
    regexpsrc += "|\\b" + *i + "\\b";
  // Compile the regular expression
  m_shregexp = FXRex(regexpsrc, REX_CAPTURE|REX_ICASE);
  // Create styles (colors for literals, comments and three sets of keywords)
  std::vector<FXColor> colors = readColors();
  createStyles(colors);
  // Initialize the checksum for detecting editor window changes
  m_shsum = 0;
}


// Regexp through the visible part of the text buffer
// and set the appropriate highlighting style
void SQLBuffer::highlightSyntax(FXbool complete)
{
  // Return if no action is necessary
  if (length == 0)
    return;

  FXString text = getText();
  FXString buf;
  FXint i;
  if (complete)
    {
      // entire file
      buf = text;
      i = 0;
    }
  else
    {
      // Return if no action is necessary (visible part is unchanged)
      if ((text.mid(getTopLine(), getBottomLine() - getTopLine()).hash() ^ getTopLine()) == m_shsum)
	return;
      m_shsum = text.mid(getTopLine(), getBottomLine() - getTopLine()).hash() ^ getTopLine();
      buf = text.mid(getTopLine(), getBottomLine()-getTopLine()), i =  getTopLine();  // visible part
    }
  
  // While a pattern is found
  FXint beg, end;
  unsigned short s;
  while (m_shregexp.match(buf, &beg, &end))
    {
      // decide and set the right style
      if ((buf[beg] == 34) || (buf[beg] == 39) || ((buf[beg] >= 48) && (buf[beg] <= 57)))  // ' " 0-9
        s = 1;
      else if ( buf[beg] == 45)  // -- SQL comment
        s = 2;
      else
	s = KEYWCOLS[buf.mid(beg, end - beg).lower()] + 2;  // SQL keyword has defined style
      if (beg > 0)
	changeStyle(i, beg, 0);           // normal style
      changeStyle(i + beg, end - beg, s);  // highlight
	// Move on in the buffer
      i += end;
      buf = buf.mid(end, buf.length()-1);
    }
  // The rest of the buffer is normal
  if (buf.length() > 0 )
    changeStyle(i, buf.length(), 0);  // normal
}
