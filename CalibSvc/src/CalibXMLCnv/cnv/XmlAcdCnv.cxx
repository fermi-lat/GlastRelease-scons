#include "XmlAcdCnv.h"

#include "GaudiKernel/CnvFactory.h"

#include "CalibData/Acd/AcdPed.h"
#include "CalibData/Acd/AcdGain.h"
#include "CalibData/Acd/AcdVeto.h"
#include "CalibData/Acd/AcdCno.h"
#include "CalibData/Acd/AcdRange.h"
#include "CalibData/Acd/AcdHighRange.h"
#include "CalibData/Acd/AcdCoherentNoise.h"

typedef XmlAcdCnv<CalibData::AcdPed> XmlAcdPedCnv;
typedef XmlAcdCnv<CalibData::AcdGain> XmlAcdGainCnv;
typedef XmlAcdCnv<CalibData::AcdVeto> XmlAcdVetoCnv;
typedef XmlAcdCnv<CalibData::AcdCno> XmlAcdCnoCnv;
typedef XmlAcdCnv<CalibData::AcdRange> XmlAcdRangeCnv;
typedef XmlAcdCnv<CalibData::AcdHighRange> XmlAcdHighRangeCnv;
typedef XmlAcdCnv<CalibData::AcdCoherentNoise> XmlAcdCoherentNoiseCnv;

static CnvFactory< XmlAcdPedCnv > s_PedFactory;
const  ICnvFactory& XmlAcdPedCnvFactory = s_PedFactory;

static CnvFactory< XmlAcdGainCnv > s_GainFactory;
const  ICnvFactory& XmlAcdGainCnvFactory = s_GainFactory;

static CnvFactory< XmlAcdVetoCnv > s_VetoFactory;
const  ICnvFactory& XmlAcdVetoCnvFactory = s_VetoFactory;

static CnvFactory< XmlAcdCnoCnv > s_CnoFactory;
const  ICnvFactory& XmlAcdCnoCnvFactory = s_CnoFactory;

static CnvFactory< XmlAcdRangeCnv > s_RangeFactory;
const  ICnvFactory& XmlAcdRangeCnvFactory = s_RangeFactory;

static CnvFactory< XmlAcdHighRangeCnv > s_HighRangeFactory;
const  ICnvFactory& XmlAcdHighRangeCnvFactory = s_HighRangeFactory;

static CnvFactory< XmlAcdCoherentNoiseCnv> s_CoherentNoiseFactory;
const  ICnvFactory& XmlAcdCoherentNoiseCnvFactory = s_CoherentNoiseFactory;

