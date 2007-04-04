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
    m_delimiters.push_back(DelimPair("!",3));
    m_delimiters.push_back(DelimPair("+",3));
    m_delimiters.push_back(DelimPair("-",3));
    m_delimiters.push_back(DelimPair("*",2));
    m_delimiters.push_back(DelimPair("/",2));
    m_delimiters.push_back(DelimPair("^",1));

    m_delimMap["("]  = " ";
    m_delimMap[")"]  = " ";
    m_delimMap["!"]  = "B";
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
    // Declare the pointer we will create
    IXTExprsnNode* pNode = 0;

    // Go through the list of of possibilities for this expression
    // First is that it is an expression with an operator
    if (pNode = parseOperator(expression)) return pNode;

    // Second is that there are parens
    int stringLen  = expression.length();
    if (expression[0] == '(' && expression[stringLen-1] == ')')
    {
        int         stringLen  = expression.length();
        std::string temp = expression.substr(1,stringLen-2);

        return parseNextExpression(temp);
    }

    // Third is that we have a function
    if (pNode = parseFunction(expression)) return pNode;

    // Fourth is that it is a constant value
    if (pNode = parseValue(expression)) return pNode;

    // Fifth is that it is a variable
    if (pNode = parseVariable(expression)) return pNode;

    // If here then we return a node which is set to zero
    pNode = new XTExprsnValue<REALNUM>("",0);

    return pNode;
}

IXTExprsnNode* XTExprsnParser::parseOperator(std::string& expression)
{
    // Declare the pointer we will create
    IXTExprsnNode* pNode = 0;

    std::string localString = expression;

    int  stringLen = localString.length();

    // Look for the next operator outside of enclosing parenthesis
    int strnLength = expression.length();
    int delimPos   = 0;

    // Search for the next delimiter
    DelimPair   fndDelimPr = findNextDelimiter(expression, delimPos);
    std::string fndDelim   = fndDelimPr.first;

    // Delimiter found, process the expression
    if (delimPos < strnLength && fndDelim != "")
    {
        IXTExprsnNode* pNodeL    = 0;
        IXTExprsnNode* pNodeR    = 0;
        std::string    inputType = "";

        // If the delimiter is not the first character then process 
        // the string to the left of the delimiter
        if (delimPos > 0)
        {
            // Substring to left of delimiter
            std::string temp = localString.substr(0, delimPos);

            // Parse it
            pNodeL = parseNextExpression(temp);

            inputType = pNodeL->getTypeId();
        }
        //else pNodeL = new XTExprsnValue<REALNUM>("",0);
        
        // If the delimiter (e.g. a ")") is not the last character
        // then process to the right of the delimiter
        if (delimPos < stringLen)
        {
            // Substring to right of delimiter
            std::string temp = localString.erase(0, delimPos+fndDelim.length());

            // Parse it
            pNodeR = parseNextExpression(temp);

            inputType = pNodeR->getTypeId();
        }
        //else pNodeR = new XTExprsnValue<REALNUM>("",0);

        // An expression node operates on a right and left node to produce 
        // some output. The type of output will depend on the operator, the 
        // type of input is extracted from the nodes (where it is assumed
        // that both nodes have the same type)
        std::string& opType = m_delimMap[fndDelim];

        // Exact type of node depends on the operator
        // Logical operation node, input type: bool or double, output type: bool
        if (opType == "B")
        {
            int pos = inputType.find("bool",0);

            if ( pos > -1 )
            {
                pNode = new XTExprsnNode<bool,bool>(fndDelim, pNodeL, pNodeR);
            }
            else if ( (pos = inputType.find("string",0)) > -1)
            {
                pNode = new XTExprsnNode<bool,std::string>(fndDelim, pNodeL, pNodeR);
            }
            else
            {
                pNode = new XTExprsnNode<bool,REALNUM>(fndDelim, pNodeL, pNodeR);
            }
        }
        // Arithmetic operation node, input type: double, output type: double
        else 
        {
            pNode = new XTExprsnNode<REALNUM,REALNUM>(fndDelim, pNodeL, pNodeR);
        }
    }

    return pNode;
}

IXTExprsnNode* XTExprsnParser::parseValue(std::string& expression)
{
    IXTExprsnNode* pNode = 0;

    // First case: expression is a constant value to convert to a double
    try
    {
        float    value  = facilities::Util::stringToDouble(expression);
        REALNUM* pValue = new REALNUM;

        *pValue = value;
        
        pNode = new XTExprsnValue<REALNUM>(expression, pValue);
    }
    // Not a constant, ends up here
    catch(facilities::WrongType&)
    {
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
    if (leftPos > 0 && rightPos > expression.size()-3)
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
            int startPos = 0;
            int endPos   = operand.length();

            std::string conExpression = findFuncArgument(operand,startPos,endPos);

            startPos = endPos + 1;
            endPos   = operand.length();
            std::string ifResult      = findFuncArgument(operand,startPos,endPos);

            startPos = endPos + 1;
            endPos   = operand.length();
            std::string elseResult    = findFuncArgument(operand,startPos,endPos);

            IXTExprsnNode* condNode = parseNextExpression(conExpression);
            IXTExprsnNode* ifNode   = parseNextExpression(ifResult);
            IXTExprsnNode* elseNode = parseNextExpression(elseResult);

            // Output type of if else presumed same whether "if" or "else" taken
            std::string inputType = ifNode->getTypeId();

            int stringPos = inputType.find("string",0);

            if (stringPos > -1)
                pNode = new XTIfElseNode<std::string>(funcCand,*condNode,*ifNode,*elseNode);
            else
                pNode = new XTIfElseNode<REALNUM>(funcCand,*condNode,*ifNode,*elseNode);
        }
        // Otherwise, must be a function
        else
        {
            IXTExprsnNode* operandNode1 = 0;
            IXTExprsnNode* operandNode2 = 0;

            // Check to see if this function has more than one operand
            int startPos = 0;
            int lastPos  = operand.length();
            int endPos   = lastPos;

            std::string argument = findFuncArgument(operand, startPos, endPos);

            // First check to see if the operand is itself a function
            // (e.g. acos(min(1.,angle))
            if (!(operandNode1 = parseFunction(argument)))
            {
                // Otherwise it must be an expression
                operandNode1 = parseNextExpression(argument);
            }

            // Is there a second argument to deal with?
            if (endPos < lastPos)
            {
                startPos     = endPos + 1;
                endPos       = lastPos;
                argument     = findFuncArgument(operand, startPos, endPos);

                // Is it a function?
                if (!(operandNode2 = parseFunction(argument)))
                {
                    // Otherwise it is an expression
                    operandNode2 = parseNextExpression(argument);
                }
            }
            else operandNode2 = new XTExprsnValue<REALNUM>("",0);

            try
            {
                pNode = new XTExprsnFunction<REALNUM>(funcCand, *operandNode1, *operandNode2);
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

IXTExprsnNode* XTExprsnParser::parseVariable(std::string& expression)
{
    // Assumption is that this is an ntuple variable which may or may not have been declared
    XTcolumnValBase* tupleVal = 0;
    
    XTtupleMap::iterator dataIter = m_tuple.find(expression);

    // Was it already declared?
    if (dataIter != m_tuple.end())
    {
        tupleVal = dataIter->second;
    }
    // We declare it here
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
            tupleVal = new XTcolumnVal<REALNUM>(expression);
        }

        m_tuple[expression] = tupleVal;
    }

    IXTExprsnNode* pNode = 0;

    if (tupleVal->getType() == "continuous")
        pNode = new XTExprsnTupleVal<XTcolumnVal<REALNUM> >(expression, dynamic_cast<XTcolumnVal<REALNUM>*>(tupleVal));
    else
        pNode = new XTExprsnTupleVal<XTcolumnVal<std::string> >(expression, dynamic_cast<XTcolumnVal<std::string>*>(tupleVal));

    return pNode;
}

XTExprsnParser::DelimPair XTExprsnParser::findNextDelimiter(const std::string& inString, int& startPos, bool checkUnary)
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

            rightOp = findNextDelimiter(tryString, rDelim, false);
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
                if (checkUnary && subStrPos == 0 && (delimOp.first == "-" || delimOp.first == "+")) continue;

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

std::string XTExprsnParser::findFuncArgument(const std::string& expression, int& startPos, int& endPos)
{
    // Look for enclosing parenthesis as a sign of an embedded function in our expression
    // (which may have an embedded comma which can screw us up)
    int leftParen  = startPos;
    int rghtParen  = endPos;
    findEnclosingParens(expression, leftParen, rghtParen);

    // Where is the comma in the expression?
    int firstComma = expression.find(",", startPos);

    // Use this loop to try to make sure we are not splitting at an embedded function
    while (leftParen > -1 && firstComma > leftParen && firstComma < rghtParen)
    {
        firstComma = expression.find(",", rghtParen);
        leftParen  = rghtParen;
        findEnclosingParens(expression, leftParen, rghtParen);
    }

    // If we have comma then set that as the end point
    if (firstComma > -1) endPos = firstComma;

    std::string result = expression.substr(startPos, endPos-startPos);

    return result;
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
