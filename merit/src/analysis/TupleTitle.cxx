// $Header$
#include "analysis/TupleTitle.h"


TupleTitle::TupleTitle(std::string head)
: m_head(head)
{
}

TupleTitle::~TupleTitle()
{
    iterator it=begin();
    while(it != end() ) delete *it++;
}

void TupleTitle::add(TitleClientBase* me){ push_back(me);}

std::string TupleTitle::operator()()const {
    std::string result(m_head);
    for( const_iterator it = begin(); it != end(); ++it)
	result += (*it)->addToTitle();
    return result;
}
