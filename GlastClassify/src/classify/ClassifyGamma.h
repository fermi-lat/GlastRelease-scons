/**@file ClassifyGamma.h
@brief 

$Header$

*/
#include "GlastClassify.h"

#include <cmath>

class ClassifyGamma : public GlastClassify
{
public:
    ClassifyGamma(const std::string& info_path)
        : GlastClassify(info_path, false) // flag that not mixed.

        // these initializers are for local filters

        ,CalEnergyRaw     ( "CalEnergyRaw") 
        ,AcdActiveDist    ( "AcdActiveDist")
        ,AcdDoca          ( "AcdDoca")         
        ,AcdRibbonActDist ( "AcdRibbonActDist")
        ,CalTrackAngle    ( "CalTrackAngle")    
        ,CalTrackDoca     ( "CalTrackDoca")     
        ,CalXtalRatio     ( "CalXtalRatio")     
        ,EvtECalTransRms  ( "EvtECalTransRms")  
        ,Tkr1FirstChisq   ( "Tkr1FirstChisq")   
        ,Tkr1ToTTrAve     ( "Tkr1ToTTrAve") 
        ,CalTotRLn        ( "CalTotRLn")
        ,FilterStatus_HI  ( "FilterStatus_HI")
        ,CTvertex         ( "CTvertex")
        ,Tkr1ZDir         ( "Tkr1ZDir")
        ,Tkr1FirstLayer   ( "Tkr1FirstLayer")
#if 0 //these are not in the tuple, but Bill wants them
        ,AcdLowerTileCount( "AcdLowerTileCount")
        ,AcdUpperTileCount( "AcdUpperTileCount")
        ,CalMaxXtalRatio  ( "CalMaxXtalRatio")  
#endif
    {

        if(      info_path =="gamma/vertex/highcal") m_case = 0;
        else if( info_path =="gamma/vertex/medcal")  m_case = 1;
        else if( info_path =="gamma/vertex/thin")    m_case = 2;
        else if( info_path =="gamma/vertex/thick")   m_case = 3;
        else if( info_path =="gamma/track/highcal")  m_case = 4;
        else if( info_path =="gamma/track/medcal")   m_case = 5;
        else if( info_path =="gamma/track/thin")     m_case = 6;
        else if( info_path =="gamma/track/thick")    m_case = 7;
        else{
            throw std::invalid_argument("unrecognized gamma case: " +info_path);
        }

    }

    //acceptance function, applied to background and signal

    virtual bool accept()
    {
        bool firstcut=
               CalEnergyRaw > 5.0 
            && CalTotRLn    > 4.0
            && FilterStatus_HI==0
            && Tkr1ZDir     < -0.3;
        if( !firstcut ) return false;

        bool  keep=false;
        switch  (m_case)
        {
        case 0: keep = vertex()  && highcal(); break;
        case 1: keep = vertex()  && medcal(); break;
        case 2: keep = vertex()  && lowcal() && thin(); break;
        case 3: keep = vertex()  && lowcal() && !thin(); break; 
        case 4: keep = !vertex() && highcal(); break;
        case 5: keep = !vertex() && medcal(); break;
        case 6: keep = !vertex() && lowcal() && thin(); break;
        case 7: keep = !vertex() && lowcal() && !thin(); break;
        };
        return keep;
    }

private:
    bool vertex()const { return CTvertex > 0.5; }
    bool highcal()const{ return CalEnergyRaw > 3500.;}
    bool medcal()const { return !highcal() && CalEnergyRaw > 350.;}
    bool lowcal()const { return CalEnergyRaw<350.; }
    bool thin()const   { return Tkr1FirstLayer >5;}

    GlastClassify::Entry CalEnergyRaw;

    GlastClassify::Entry AcdActiveDist    ;
    GlastClassify::Entry AcdDoca          ;
    GlastClassify::Entry AcdRibbonActDist ;
#if 0 // not in tuple
    GlastClassify::Entry AcdLowerTileCount;
    GlastClassify::Entry AcdUpperTileCount;
    GlastClassify::Entry CalMaxXtalRatio  ;
#endif
    GlastClassify::Entry CalTrackAngle    ;
    GlastClassify::Entry CalTrackDoca     ;
    GlastClassify::Entry CalXtalRatio     ;
    GlastClassify::Entry EvtECalTransRms  ;
    GlastClassify::Entry Tkr1FirstChisq   ;
    GlastClassify::Entry Tkr1ToTTrAve     ;

    GlastClassify::Entry CalTotRLn        ;
    GlastClassify::Entry FilterStatus_HI  ;
    GlastClassify::Entry CTvertex         ;
    GlastClassify::Entry Tkr1ZDir         ;
    GlastClassify::Entry Tkr1FirstLayer   ;

    int m_case;
};

