#ifndef G4GenErrorRecord_H
#define G4GenErrorRecord_H

class MsgStream ;
#include <vector>
#include <string>

class G4GenErrorRecord 
{
public:
    G4GenErrorRecord(int run, 
                     int event, 
                     double lastTime,
                     const std::string &catcherName,
                     const std::string & comment )
                     : m_run(run), 
                       m_event(event), 
                       m_lastTime(lastTime),
                       m_catcherName(catcherName),
                       m_comment(comment) {}
    ~G4GenErrorRecord() {}
      
    int                getRun()         {return m_run;         }
    int                getEvent()       {return m_event;       }
    double             getLastTime()    {return m_lastTime;    }
    const std::string& getCatcherName() {return m_catcherName; }
    const std::string& getComment()     {return m_comment;     }
     
    void print( MsgStream & ) const ;

private:
    int         m_run;
    int         m_event;
    double      m_lastTime;
    std::string m_catcherName;
    std::string m_comment;
} ;

#endif
