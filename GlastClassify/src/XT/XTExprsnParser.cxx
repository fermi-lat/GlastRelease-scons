/**@file XTExprsnParser.cxx

@brief implementation of class XTExprsnParser

$Header$
*/

#include "XTExprsnParser.h"
#include "facilities/Util.h"
    
XTExprsnParser::XTExprsnParser(XTtupleMap& tuple) : m_tuple(tuple) 
{
    m_delimiters.clear();

    m_delimiters.push_back(DelimPair("|",7));
    m_delimiters.push_back(DelimPair("&",6));
    m_delimiters.push_back(DelimPair("==",5));
    m_delimiters.push_back(DelimPair("!=",5));
    m_delimiters.push_back(DelimPair(">=",4));
    m_delimiters.push_back(DelimPair("<=",4));
    m_delimiters.push_back(DelimPair(">",4));
    m_delimiters.push_back(DelimPair("<",4));
    m_delimiters.push_back(DelimPair("+",3));
    m_delimiters.push_back(DelimPair("-",3));
    m_delimiters.push_back(DelimPair("*",2));
    m_delimiters.push_back(DelimPair("/",2));
    m_delimiters.push_back(DelimPair("^",1));

    m_delimMap["("]  = " ";
    m_delimMap[")"]  = " ";
    m_delimMap["&"]  = "B";
    m_delimMap["|"]  = "B";
    m_delimMap[">"]  = "B";
    m_delimMap[">="] = "B";
    m_delimMap["<"]  = "B";
    m_delimMap["<="] = "B";
    m_delimMap["=="] = "B";
    m_delimMap["!="] = "B";
    m_delimMap["^"]  = "D";
    m_delimMap["*"]  = "D";
    m_delimMap["/"]  = "D";
    m_delimMap["+"]  = "D";
    m_delimMap["-"]  = "D";
}

IXTExprsnNode* XTExprsnParser::parseExpression(std::string& expression)
{
    expression = trimTrailing(expression);
    expression = trimCharacters(expression, ' ');

    // Can we get rid of the embedded special characters?
    for(int idx = 0; idx < 4; idx++)
    {
        static const char tabChar[]={0x9,0xA,0xC,0xD};
        expression = trimCharacters(expression, tabChar[idx]);
    }

    //return parseNextExpression(parsedExpression, expression);
    IXTExprsnNode* node = parseNextExpression(expression);
    //node->print();
    return node;
}


IXTExprsnNode* XTExprsnParser::parseNextExpression(std::string& expression)
{
    // Set null return value
    IXTExprsnNode* pNode = 0;

    std::string localString = expression;

    int  stringLen  = localString.length();

    // Look for the next operator outside of enclosing parenthesis
    int strnLength  = expression.length();
    int startPos    = 0;

    // Search for the next delimiter
    DelimPair   fndDelimPr = findNextDelimiter(expression, startPos);
    std::string fndDelim   = fndDelimPr.first;

    // Delimiter found
    if (startPos < strnLength && fndDelim != "")
    {
        IXTExprsnNode* pNodeL = 0;
        IXTExprsnNode* pNodeR = 0;

        // If the delimiter is not the first character then process 
        // the string to the left of the delimiter
        if (startPos > 0)
        {
            // Substring to left of delimiter
            std::string temp = localString.substr(0, startPos);

            // Parse it
            pNodeL = parseNextExpression(temp);
        }
        
        // If the delimiter (e.g. a ")") is not the last character
        // then process to the right of the delimiter
        if (startPos < stringLen)
        {
            // Substring to right of delimiter
            std::string temp = localString.erase(0, startPos+fndDelim.length());

            // Parse it
            pNodeR = parseNextExpression(temp);
        }

        // If only a left or right the return
        if      (pNodeL != 0 && pNodeR == 0) pNode = pNodeL;
        else if (pNodeL == 0 && pNodeR != 0) pNode = pNodeR;
        // Otherwise create new Expression Node
        else
        {
            // An expression node operates on a right and left node to produce 
            // some output. The type of output will depend on the operator, the 
            // type of input is extracted from the nodes (where it is assumed
            // that both nodes have the same type)
            std::string& opType = m_delimMap[fndDelim];

            std::string inputType = pNodeL->getTypeId();

            // Exact type of node depends on the operator
            // Logical operation node, input type: bool or double, output type: bool
            if (opType == "B")
            {
                int pos = inputType.find("bool",0);

                if ( pos > -1 )
                {
                    pNode = new XTExprsnNode<bool,bool>(fndDelim, *pNodeL, *pNodeR);
                }
                else if ( (pos = inputType.find("string",0)) > -1)
                {
                    pNode = new XTExprsnNode<bool,std::string>(fndDelim, *pNodeL, *pNodeR);
                }
                else
                {
                    pNode = new XTExprsnNode<bool,double>(fndDelim, *pNodeL, *pNodeR);
                }
            }
            // Arithmetic operation node, input type: double, output type: double
            else 
            {
                pNode = new XTExprsnNode<double,double>(fndDelim, *pNodeL, *pNodeR);
            }
        }
    }
    // Look for special case of string enclosed by parens
    else if (localString[0] == '(' && localString[strnLength-1] == ')')
    {
        std::string temp = localString.substr(1,strnLength-2);

        pNode = parseNextExpression(temp);
    }
    // Must be a value of some sort
    else 
    {
        pNode = parseValue(localString);
    }

    return pNode;
}

IXTExprsnNode* XTExprsnParser::parseValue(std::string& expression)
{
    IXTExprsnNode* pNode = 0;

    // First case: expression is a constant value to convert to a double
    try
    {
        double  value  = facilities::Util::stringToDouble(expression);
        double* pValue = new double;

        *pValue = value;
        
        pNode = new XTExprsnValue<double>(expression, pValue);
    }
    // Not a constant, ends up here
    catch(facilities::WrongType&)
    {
        // Look up possible functions
        if (!(pNode = parseFunction(expression)))
        {
            // Last posibility: an ntuple variable
            //XTcolumnVal<double>*  tupleVal = m_tuple.addNewDataItem(expression);
            XTcolumnValBase* tupleVal = 0;
            
            XTtupleMap::iterator dataIter = m_tuple.find(expression);

            if (dataIter != m_tuple.end())
            {
                tupleVal = dataIter->second;
            }
            else
            {
                // Ok, last last possibility is that it is a string constant...
                int stringLen  = expression.length();
                int firstQuote = expression.find("\"", 0);
                int secndQuote = expression.find("\"", firstQuote+1);

                if (firstQuote == 0 && secndQuote == stringLen-1)
                {
                    std::string dataValue = expression.substr(firstQuote+1,secndQuote-1);
                    XTcolumnVal<std::string>* newValue = new XTcolumnVal<std::string>(dataValue,"categorical");
                    newValue->setDataValue(dataValue);
                    tupleVal = newValue;
                }
                else
                {
                    tupleVal = new XTcolumnVal<double>(expression);
                }

                m_tuple[expression] = tupleVal;
            }

            if (tupleVal->getType() == "continuous")
                pNode = new XTExprsnTupleVal<XTcolumnVal<double> >(expression, dynamic_cast<XTcolumnVal<double>*>(tupleVal));
            else
                pNode = new XTExprsnTupleVal<XTcolumnVal<std::string> >(expression, dynamic_cast<XTcolumnVal<std::string>*>(tupleVal));
        }
    }

    return pNode;
}
    
IXTExprsnNode* XTExprsnParser::parseFunction(std::string& expression)
{
    IXTExprsnNode* pNode = 0;

    int leftPos = 0;
    int rightPos = expression.length();

    std::string operand = findEnclosingParens(expression, leftPos, rightPos);

    // Is there a possible function here?
    if (leftPos > 0)
    {
        std::string funcCand = expression.substr(0,leftPos);

        // Special case of IM asking to retrieve data from tuple
        if (funcCand == "get" || funcCand == "\"Pr" || funcCand == "Pr")
        {
            expression = operand;
            pNode = parseFunction(expression);
        }
        // Special case of an IM "ifelse" clause
        else if (funcCand == "ifelse")
        {
            // operand will contain the conditional expression and the two results
            int firstPos   = 1;
            int operandLen = operand.length();
            int firstComma = operand.find(",", firstPos);
            int secndComma = operand.find(",", firstComma+1);

            std::string conExpression = operand.substr(0, firstComma);
            std::string ifResult      = operand.substr(firstComma+1,secndComma-firstComma-1);
            std::string elseResult    = operand.substr(secndComma+1,operandLen-secndComma-1);

            IXTExprsnNode* condNode = parseNextExpression(conExpression);
            IXTExprsnNode* ifNode   = parseNextExpression(ifResult);
            IXTExprsnNode* elseNode = parseNextExpression(elseResult);

            // Output type of if else presumed same whether "if" or "else" taken
            std::string inputType = ifNode->getTypeId();

            int stringPos = inputType.find("string",0);

            if (stringPos > -1)
                pNode = new XTIfElseNode<std::string>(funcCand,*condNode,*ifNode,*elseNode);
            else
                pNode = new XTIfElseNode<double>(funcCand,*condNode,*ifNode,*elseNode);
        }
        // Otherwise, must be a function
        else
        {
            IXTExprsnNode* operandNode1 = 0;
            IXTExprsnNode* operandNode2 = 0;

            // First check to see if the operand is itself a function
            // (e.g. acos(min(1.,angle))
            if (operandNode1 = parseFunction(operand))
            {
                operandNode2 = new XTExprsnValue<double>("",0);
            }
            // Ok, standard argument
            else
            {
                // Check to see if this function has more than one operand
                int firstPos = 0;
                int commaPos = operand.find(",",firstPos);

             if (commaPos > 0)
                {
                    std::string operand1 = operand.substr(firstPos, commaPos);
                    operandNode1 = parseNextExpression(operand1);

                    std::string operand2 = operand.substr(commaPos+1, operand.length()-commaPos);
                    operandNode2 = parseNextExpression(operand2);
                }
                else
                {
                    operandNode1 = parseNextExpression(operand);
                    operandNode2 = new XTExprsnValue<double>("",0);
                }
            }

            try
            {
                pNode = new XTExprsnFunction<double>(funcCand, *operandNode1, *operandNode2);
            }
            catch(...)
            {
                pNode = operandNode1;
                delete operandNode2;
            }
        }
    }

    return pNode;
}

XTExprsnParser::DelimPair XTExprsnParser::findNextDelimiter(const std::string& inString, int& startPos)
{
    // null string in case we don't find anything
    DelimPair fndDelim("",0);
    int       stringLen = inString.length();

    // Check to find a matching paren pair
    int leftPos  = startPos;
    int rightPos = stringLen;
    std::string parenArg = findEnclosingParens(inString, leftPos, rightPos);

    // if no parens then look for an embedded string value (a "categorical" variable value)
    if (parenArg == "")
    {
        leftPos  = startPos;
        rightPos = stringLen;
        parenArg = findCategoricalVal(inString, leftPos, rightPos);
    }

    // If enclosing parenthesis then check to the left and right
    if (parenArg != "")
    {
        DelimPair leftOp  = DelimPair("",0);
        int       lDelim  = 0;
        DelimPair rightOp = DelimPair("",0);
        int       rDelim  = 0;

        std::string test1 = inString.substr(leftPos,1);
        std::string test2 = inString.substr(rightPos,1);

        // Check if parens are to the right of the next delimiter
        if (leftPos > startPos)
        {
            std::string tryString = inString.substr(0, leftPos);

            leftOp = findNextDelimiter(tryString, lDelim);
        }

        // Now check if to the left
        if (rightPos < stringLen)
        {
            std::string tryString = inString.substr(rightPos+1, stringLen - rightPos);

            rightOp = findNextDelimiter(tryString, rDelim);
        }

        // Which do we take?
        fndDelim = leftOp;
        startPos = lDelim;

        if (rightOp.second > leftOp.second)
        {
            fndDelim = rightOp;
            startPos = rightPos + rDelim + 1;
        }
    }

    // If here we either have a valid delimiter about to appear OR we have hit 
    // the end of the string... 
    else if (startPos < stringLen)
    {
        for(DelimVector::iterator delIter = m_delimiters.begin(); delIter != m_delimiters.end(); delIter++)
        {
            const DelimPair& delimOp = *delIter;

            int subStrPos = inString.find(delimOp.first, startPos);

            // position in string > -1 if we have a match
            if (subStrPos > -1)
            {
                // Attempt to catch special case of unary + or - operator
                if (subStrPos == 0 && (delimOp.first == "-" || delimOp.first == "+")) continue;

                // Ugliness to check for exponential notation
                if (subStrPos > 0 && (delimOp.first == "-" || delimOp.first == "+"))
                {
                    // Ok, this will fix the problem for when we have just an exponentiated number
                    // lying around. It does not fix the problem when we have arithmetic involving
                    // exponentiated numbers... So, here is a place that needs more work.
                    int expPos = inString.find("e", subStrPos-1);

                    // Possible exponentiated number
                    if (expPos >= 0)
                    {
                        try
                        {
                            double value = facilities::Util::stringToDouble(inString);

                            // Was able to read number so "continue" past this delimiter
                            continue;
                        }
                        // Was not a valid number, catch here and go to next step
                        catch(facilities::WrongType&) {}
                    }
                }
                
                // Valid delimiter! 
                fndDelim = delimOp;
                startPos = subStrPos;
                break;
            }
        }
    }

    return fndDelim;
}

std::string XTExprsnParser::findEnclosingParens(const std::string& expression, int& startPos, int& endPos)
{
    // This method aims to return the location of matching parenthesis
    // contained within the string "expression". 
    // null string in case we don't find anything
    std::string subStr = "";
    
    // Initialize end... 
    endPos = expression.length();

    // Special section to handle parenthesis... (not pretty...)
    startPos = expression.find("(", startPos);

    if (startPos > -1)
    {
        // Ugly code starts here for looking for matching right paren
        int skipNparens = 1;
        int stringLen   = endPos;

        endPos = startPos + 1;

        // Loop through characters and match this "(" to its partner ")"
        while(endPos < stringLen)
        {
            std::string nextChar = expression.substr(endPos,1);

            // Beginning parenthesis so we need to skip
            if (nextChar == "(") skipNparens++;

            // Found ending parens, decrement paren depth count
            if (nextChar == ")") 
            {
                skipNparens--;
                if (skipNparens == 0) break;
            }
        
            // Increment counter
            endPos++;
        }

        // Set the return string minus enclosing parens
        subStr = expression.substr(startPos+1, endPos-startPos-1);
    }

    return subStr;
}

std::string XTExprsnParser::findCategoricalVal(const std::string& expression, int& startPos, int& endPos)
{
    // This method aims to return the location of a "categorical" value 
    // contained within the string "expression". 
    // null string in case we don't find anything
    std::string subStr = "";
    
    // Initialize end... 
    endPos = expression.length();

    // Special section to handle parenthesis... (not pretty...)
    startPos = expression.find("\"", startPos);

    // if we find a left " then look for the next one following it
    if (startPos > -1)
    {
        endPos = expression.find("\"",startPos+1);

        // Set the return string including the surrounding quotes
        subStr = expression.substr(startPos, endPos-startPos+1);
    }

    return subStr;
}


// Remove all blank spaces from a given string
//std::string XTExprsnParser::trimCharacters(std::string& expression, const std::string& charToTrim)
std::string XTExprsnParser::trimCharacters(std::string& expression, const char& charToTrim)
{
    std::string trimmedString = expression;
    int blankPos = trimmedString.find(charToTrim,0);

    if (blankPos > -1) 
    {
        std::string temp = trimmedString.erase(blankPos, 1);

        trimmedString = trimCharacters(trimmedString, charToTrim);
    }

    return trimmedString;
}
    
std::string XTExprsnParser::trimTrailing(std::string& expression)
{
    std::string temp = expression;
    facilities::Util::trimTrailing(&temp);
    return temp;
}
