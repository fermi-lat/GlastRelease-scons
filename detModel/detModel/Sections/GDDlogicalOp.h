#ifndef GDDLOGICALOP_H
#define GDDLOGICALOP_H
#include "detModel/Sections/GDDvolume.h"
class GDDanySinglePosition;
enum logicaltype{unions,intersection,subtraction}
//union is inusable (reserved word)
class GDDlogicalOp : GDDvolume {
private:    
    logicaltype operation;
};
#endif //GDDLOGICALOP_H


/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDlogicalOp : GDDvolume {
public:

    GDDlogicalOp(logicaltype plogictype,typeVolume ptype,std::string pname);

    void addPosition(GDDanySinglePosition* ppos);
private:    
    logicaltype logictype;

    /** @link aggregation
     *  @supplierCardinality 1..n */
    vector < GDDanySinglePosition * > anySinglePosition;
};


