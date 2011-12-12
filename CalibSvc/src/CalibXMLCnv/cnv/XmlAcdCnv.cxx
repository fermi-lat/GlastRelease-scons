#include "XmlAcdCnv.h"

#include "GaudiKernel/CnvFactory.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"
#include "CalibData/Acd/AcdRibbon.h"
#include "CalibData/Acd/AcdHighPed.h"
#include "CalibData/Acd/AcdCarbon.h"
#include "CalibData/Acd/AcdVetoFit.h"
#include "CalibData/Acd/AcdCnoFit.h"
#include "CalibData/Acd/AcdPE.h"

typedef XmlAcdCnv<CalibData::AcdPed> XmlAcdPedCnv;
typedef XmlAcdCnv<CalibData::AcdGain> XmlAcdGainCnv;
typedef XmlAcdCnv<CalibData::AcdVeto> XmlAcdVetoCnv;
typedef XmlAcdCnv<CalibData::AcdCno> XmlAcdCnoCnv;
typedef XmlAcdCnv<CalibData::AcdRange> XmlAcdRangeCnv;
typedef XmlAcdCnv<CalibData::AcdHighRange> XmlAcdHighRangeCnv;
typedef XmlAcdCnv<CalibData::AcdCoherentNoise> XmlAcdCoherentNoiseCnv;
typedef XmlAcdCnv<CalibData::AcdRibbon> XmlAcdRibbonCnv;
typedef XmlAcdCnv<CalibData::AcdHighPed> XmlAcdHighPedCnv;
typedef XmlAcdCnv<CalibData::AcdCarbon> XmlAcdCarbonCnv;
typedef XmlAcdCnv<CalibData::AcdVetoFit> XmlAcdVetoFitCnv;
typedef XmlAcdCnv<CalibData::AcdCnoFit> XmlAcdCnoFitCnv;
typedef XmlAcdCnv<CalibData::AcdPE> XmlAcdPECnv;

//static CnvFactory< XmlAcdPedCnv > s_PedFactory;
//const  ICnvFactory& XmlAcdPedCnvFactory = s_PedFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdPedCnv);

//static CnvFactory< XmlAcdGainCnv > s_GainFactory;
//const  ICnvFactory& XmlAcdGainCnvFactory = s_GainFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdGainCnv);

//static CnvFactory< XmlAcdVetoCnv > s_VetoFactory;
//const  ICnvFactory& XmlAcdVetoCnvFactory = s_VetoFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdVetoCnv);

//static CnvFactory< XmlAcdCnoCnv > s_CnoFactory;
//const  ICnvFactory& XmlAcdCnoCnvFactory = s_CnoFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdCnoCnv);

//static CnvFactory< XmlAcdRangeCnv > s_RangeFactory;
//const  ICnvFactory& XmlAcdRangeCnvFactory = s_RangeFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdRangeCnv);

//static CnvFactory< XmlAcdHighRangeCnv > s_HighRangeFactory;
//const  ICnvFactory& XmlAcdHighRangeCnvFactory = s_HighRangeFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdHighRangeCnv);

//static CnvFactory< XmlAcdCoherentNoiseCnv > s_CoherentNoiseFactory;
//const  ICnvFactory& XmlAcdCoherentNoiseCnvFactory = s_CoherentNoiseFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdCoherentNoiseCnv);

//static CnvFactory< XmlAcdRibbonCnv > s_RibbonFactory;
//const  ICnvFactory& XmlAcdRibbonCnvFactory = s_RibbonFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdRibbonCnv);


//static CnvFactory< XmlAcdHighPedCnv > s_HighPedFactory;
//const  ICnvFactory& XmlAcdHighPedCnvFactory = s_HighPedFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdHighPedCnv);


//static CnvFactory< XmlAcdCarbonCnv > s_CarbonFactory;
//const  ICnvFactory& XmlAcdCarbonCnvFactory = s_CarbonFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdCarbonCnv);


//static CnvFactory< XmlAcdVetoFitCnv > s_VetoFitFactory;
//const  ICnvFactory& XmlAcdVetoFitCnvFactory = s_VetoFitFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdVetoFitCnv);


//static CnvFactory< XmlAcdCnoFitCnv > s_CnoFitFactory;
//const  ICnvFactory& XmlAcdCnoFitCnvFactory = s_CnoFitFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdCnoFitCnv);


//static CnvFactory< XmlAcdPECnv > s_PEFactory;
//const  ICnvFactory& XmlAcdPECnvFactory = s_PEFactory;
DECLARE_CONVERTER_FACTORY(XmlAcdPECnv);

