// $Header$

#include "rdbModel/Tables/Index.h"

namespace rdbModel {

  Visitor::VisitorState  Index::accept(Visitor * v) {
    Visitor::VisitorState state = v->visitIndex(this);
    if (state == Visitor::VBRANCHDONE) return Visitor::VCONTINUE;
    return state;
  }

  const std::vector<std::string>& Index::getColumnNames() {
    return m_indexCols;
  }


}
