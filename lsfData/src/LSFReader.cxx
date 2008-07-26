#include "eventFile/LSE_Context.h"
#include "eventFile/EBF_Data.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/LPA_Handler.h"
#include "eventFile/LSE_Keys.h"

#include "lsfData/LSFReader.h"
#include "lsfData/LsfCcsds.h"
#include "lsfData/LsfDatagramInfo.h"
#include "lsfData/LsfGemScalers.h"
#include "lsfData/LsfGemTime.h"
#include "lsfData/LsfRunInfo.h"
#include "lsfData/LsfTime.h"
#include "lsfData/LsfTimeTone.h"
#include "lsfData/LsfMetaEvent.h"

namespace lsfData {
  
  bool LSFReader::read( LsfCcsds& lccsds, MetaEvent& lmeta, eventFile::EBF_Data& ebf )
  {
    // create LSE objects to hold the retrieved values
    eventFile::LSE_Context        ctx;
    eventFile::LSE_Info::InfoType infotype;
    eventFile::LPA_Info           pinfo;
    eventFile::LCI_ACD_Info       ainfo;
    eventFile::LCI_CAL_Info       cinfo;
    eventFile::LCI_TKR_Info       tinfo;
    eventFile::LSE_Keys::KeysType ktype;
    eventFile::LPA_Keys           pakeys;
    eventFile::LCI_Keys           cikeys;

    // read the native objects
    if ( !eventFile::LSEReader::read( ctx, ebf, infotype, pinfo, ainfo, cinfo, tinfo, ktype, pakeys, cikeys ) ) {
      return false;
    }

    // transfer the CCSDS information
    transferCcsds( ctx, lccsds );

    // transfer the context information
    transferContext( ctx, lmeta );

    // transfer the type-specific meta-information
    switch ( infotype ) {
    case eventFile::LSE_Info::LPA:
      transferInfo( ctx, pinfo, lmeta );
      transferKeys( pakeys, lmeta );
      break;
    case eventFile::LSE_Info::LCI_ACD:
      transferInfo( ctx, ainfo, lmeta );
      transferKeys( cikeys, lmeta );
      break;
    case eventFile::LSE_Info::LCI_CAL:
      transferInfo( ctx, cinfo, lmeta );
      transferKeys( cikeys, lmeta );
      break;
    case eventFile::LSE_Info::LCI_TKR:
      transferInfo( ctx, tinfo, lmeta );
      transferKeys( cikeys, lmeta );
      break;
    default:
      break;
    }

    return true;
  }

  void LSFReader::transferCcsds( const eventFile::LSE_Context& ctx, LsfCcsds& lccsds )
  {
    lccsds.initialize( ctx.ccsds.scid, ctx.ccsds.apid, ctx.ccsds.utc );
  }
  
  void LSFReader::transferContext( const eventFile::LSE_Context& ctx, MetaEvent& lsfmeta )
  {
    // set the datagram information
    enums::Lsf::Open::Action  ao = static_cast< enums::Lsf::Open::Action  >( ctx.open.action );
    enums::Lsf::Open::Reason  ro = static_cast< enums::Lsf::Open::Reason  >( ctx.open.reason );
    enums::Lsf::Crate         cr = static_cast< enums::Lsf::Crate         >( ctx.open.crate );
    enums::Lsf::Mode          md = static_cast< enums::Lsf::Mode          >( ctx.open.mode );
    enums::Lsf::Close::Action ac = static_cast< enums::Lsf::Close::Action >( ctx.close.action );
    enums::Lsf::Close::Reason rc = static_cast< enums::Lsf::Close::Reason >( ctx.close.reason );
    DatagramInfo dgm( ao, ro, cr, md, ac, rc,
			       ctx.open.datagrams, ctx.open.modeChanges 
			       );
    lsfmeta.setDatagram( dgm );

    // set the run information
    enums::Lsf::Platform   pl = static_cast< enums::Lsf::Platform   >( ctx.run.platform );
    enums::Lsf::DataOrigin od = static_cast< enums::Lsf::DataOrigin >( ctx.run.origin );
    RunInfo  run( pl, od,
			   ctx.run.groundId, ctx.run.startedAt, runid()
			   );
    lsfmeta.setRun( run );

    // set the GEM scalers
    GemScalers sca( ctx.scalers.elapsed, ctx.scalers.livetime,
			  ctx.scalers.prescaled, ctx.scalers.discarded,
			  ctx.scalers.sequence, ctx.scalers.deadzone
			  );
    lsfmeta.setScalers( sca );
  }
  
  void LSFReader::transferTime( const eventFile::LSE_Context& ctx, const eventFile::LSE_Info& info, MetaEvent& lsfmeta )
  {
    // get the current timetone info
    unsigned char flags = 0;
    flags |= ( ctx.current.missingTimeTone ) ? enums::Lsf::TimeTone::MISSING_TIMETONE_MASK : 0x0;
    flags |= ( ctx.current.missingLatPps )   ? enums::Lsf::TimeTone::MISSING_LAT_MASK      : 0x0;
    flags |= ( ctx.current.missingCpuPps )   ? enums::Lsf::TimeTone::MISSING_CPU_MASK      : 0x0;
    flags |= ( ctx.current.earlyEvent )      ? enums::Lsf::TimeTone::EARLY_EVENT_MASK      : 0x0;
    flags |= ( ctx.current.sourceGps )      ? enums::Lsf::TimeTone::SOURCE_GPS_MASK      : 0x0;
    TimeTone curr( ctx.current.incomplete, ctx.current.timeSecs, 
			    ctx.current.flywheeling, flags,
			    GemTime( ctx.current.timeHack.hacks, ctx.current.timeHack.tics )
			    );

    // get the previous timetone info
    flags = 0;
    flags |= ( ctx.previous.missingTimeTone ) ? enums::Lsf::TimeTone::MISSING_TIMETONE_MASK : 0x0;
    flags |= ( ctx.previous.missingLatPps )   ? enums::Lsf::TimeTone::MISSING_LAT_MASK      : 0x0;
    flags |= ( ctx.previous.missingCpuPps )   ? enums::Lsf::TimeTone::MISSING_CPU_MASK      : 0x0;
    flags |= ( ctx.previous.earlyEvent )      ? enums::Lsf::TimeTone::EARLY_EVENT_MASK      : 0x0;
    flags |= ( ctx.previous.sourceGps )      ? enums::Lsf::TimeTone::SOURCE_GPS_MASK      : 0x0;
    TimeTone prev( ctx.previous.incomplete, ctx.previous.timeSecs, 
			    ctx.previous.flywheeling, flags,
			    GemTime( ctx.previous.timeHack.hacks, ctx.previous.timeHack.tics )
			    );

    // set the current, previous, and event-time fields
    lsfmeta.setTime( Time( curr, prev,
				    GemTime( info.timeHack.hacks, info.timeHack.tics ),
				    info.timeTics )
		     );
  }

  void LSFReader::transferLciCfg( const eventFile::LCI_Info& info, LciConfiguration& lcfg )
  {
    unsigned char flags(0);
    flags |= ( info.autoRange      ) ? enums::Lsf::AUTORANGE        : 0x0;
    flags |= ( info.zeroSupression ) ? enums::Lsf::ZERO_SUPPRESSION : 0x0;
    flags |= ( info.strobe         ) ? enums::Lsf::STROBE           : 0x0;
    lcfg.set( info.softwareKey, 
	      info.writeCfg,
	      info.readCfg,
	      info.periodicPrescale,
	      flags
	      );
  }

  void LSFReader::transferKeys( const eventFile::LPA_Keys& pakeys, MetaEvent& lmeta )
  {
    // create & populate a local LpaKeys object
    LpaKeys lkeys( pakeys.LATC_master, pakeys.LATC_ignore, pakeys.SBS,
                   pakeys.LPA_db );

    // install the keys object into the MetaEvent
    lmeta.setKeys( lkeys );
  }

  void LSFReader::transferInfo( const eventFile::LSE_Context& ctx, const eventFile::LPA_Info& info, MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local Configuration object
    LpaConfiguration lcfg( info.hardwareKey, info.softwareKey );

    // install the configuration object into the MetaEvent
    lmeta.setConfiguration( lcfg );

    std::vector<eventFile::LPA_Handler>::const_iterator handlerIt;
    for (handlerIt = info.handlers.begin(); handlerIt != info.handlers.end(); handlerIt++) {
    const eventFile::PassthruHandlerRsdV0* evtPass( handlerIt->passthruRsdV0() );
        if (((enums::Lsf::HandlerId)(handlerIt->id)) == enums::Lsf::DGN) {
            lsfData::DgnHandler dgn;
            dgn.set(handlerIt->masterKey,handlerIt->cfgKey,handlerIt->cfgId,
                (enums::Lsf::RsdState)(handlerIt->state),
                (enums::Lsf::LeakedPrescaler)(handlerIt->prescaler),
                handlerIt->version,
                (enums::Lsf::HandlerId)(handlerIt->id),handlerIt->has);
            const eventFile::DgnHandlerRsdV0*  evtDgn( handlerIt->dgnRsdV0() );
            if (evtDgn) {
                dgn.setStatus(evtDgn->status);
            }
            lmeta.addDgnHandler(dgn);
        } else if (((enums::Lsf::HandlerId)(handlerIt->id)) == enums::Lsf::GAMMA) {
            lsfData::GammaHandler gam;
            gam.set(handlerIt->masterKey,handlerIt->cfgKey,handlerIt->cfgId,
                (enums::Lsf::RsdState)(handlerIt->state),
                (enums::Lsf::LeakedPrescaler)(handlerIt->prescaler),
                handlerIt->version,
                (enums::Lsf::HandlerId)(handlerIt->id),handlerIt->has);
            const eventFile::GammaHandlerRsdV0*  evtGamma( handlerIt->gammaRsdV0() );
            if (evtGamma)
                gam.setStatus(evtGamma->status,evtGamma->stage(),
                evtGamma->energyValid, evtGamma->energyInLeus);
            lmeta.addGammaHandler(gam);
        } else if (((enums::Lsf::HandlerId)(handlerIt->id)) == enums::Lsf::HIP) {
            lsfData::HipHandler hip;
            hip.set(handlerIt->masterKey,handlerIt->cfgKey,handlerIt->cfgId,
                (enums::Lsf::RsdState)(handlerIt->state),
                (enums::Lsf::LeakedPrescaler)(handlerIt->prescaler),
                handlerIt->version,
                (enums::Lsf::HandlerId)(handlerIt->id),handlerIt->has);
            const eventFile::HipHandlerRsdV0*  evtHip( handlerIt->hipRsdV0() );
            if (evtHip)
                hip.setStatus(evtHip->status);
            lmeta.addHipHandler(hip);
        } else if (((enums::Lsf::HandlerId)(handlerIt->id)) == enums::Lsf::MIP) {
            lsfData::MipHandler mip;
            mip.set(handlerIt->masterKey,handlerIt->cfgKey,handlerIt->cfgId,
                (enums::Lsf::RsdState)(handlerIt->state),
                (enums::Lsf::LeakedPrescaler)(handlerIt->prescaler),
                handlerIt->version,
                (enums::Lsf::HandlerId)(handlerIt->id),handlerIt->has);
            const eventFile::MipHandlerRsdV0* evtMip( handlerIt->mipRsdV0() );
            if (evtMip)
                mip.setStatus(evtMip->status);
            lmeta.addMipHandler(mip);
        } else if (((enums::Lsf::HandlerId)(handlerIt->id)) == enums::Lsf::PASS_THRU) {
            lsfData::PassthruHandler pass;
            pass.set(handlerIt->masterKey,handlerIt->cfgKey,handlerIt->cfgId,
                (enums::Lsf::RsdState)(handlerIt->state),
                (enums::Lsf::LeakedPrescaler)(handlerIt->prescaler),
                handlerIt->version,
                (enums::Lsf::HandlerId)(handlerIt->id),handlerIt->has);
            const eventFile::PassthruHandlerRsdV0* evtPass( handlerIt->passthruRsdV0() );
            if (evtPass) 
                pass.setStatus(evtPass->status);
            lmeta.addPassthruHandler(pass);
        } else {
          /* do nothing for now
            lsfData::LpaHandler lpa;
            lpa.set(handlerIt->masterKey,handlerIt->cfgKey,handlerIt->cfgId,
                    (enums::Lsf::RsdState)(handlerIt->state),
                    (enums::Lsf::LeakedPrescaler)(handlerIt->prescaler),
                    handlerIt->version,
                    (enums::Lsf::HandlerId)(handlerIt->id),handlerIt->has);
            //lmeta.addLpaHandler((enums::Lsf::HandlerId)handlerIt->id, lpa);
            lmeta.addLpaHandler(lpa);
          */
           
        }
    }
  }

  void LSFReader::transferKeys( const eventFile::LCI_Keys& cikeys, MetaEvent& lmeta )
  {
    // create & populate a local LciKeys object
    LciKeys lkeys( cikeys.LATC_master, cikeys.LATC_ignore, cikeys.LCI_script );
    
    // install the keys object into the MetaEvent
    lmeta.setKeys( lkeys );
  }

  void LSFReader::transferInfo( const eventFile::LSE_Context& ctx, const eventFile::LCI_ACD_Info& info, MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local Configuration object
    LciAcdConfiguration lcfg( info.injected, 
				       info.threshold, 
				       info.biasDac, 
				       info.holdDelay,
                                       info.hitmapDelay,
                                       info.range,
				       LciAcdConfiguration::AcdTrigger( info.trigger.veto, 
										 info.trigger.vetoVernier,
										 info.trigger.hld ),
				       Channel( info.channel.single, 
							 info.channel.all, 
							 info.channel.latc )
				       );
    transferLciCfg( info, lcfg );

    // install the configuration object into the MetaEvent
    lmeta.setConfiguration( lcfg );

  }

  void LSFReader::transferInfo( const eventFile::LSE_Context& ctx, const eventFile::LCI_CAL_Info& info, MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local Configuration object
    LciCalConfiguration lcfg( info.uld,
				       info.injected,
				       info.delay,
                                       info.firstRange,
				       info.threshold,
                                       info.calibGain,
                                       info.highCalEna,
                                       info.highRngEna,
                                       info.highGain,
                                       info.lowCalEna,
                                       info.lowRngEna,
                                       info.highGain,
				       LciCalConfiguration::CalTrigger( info.trigger.le, info.trigger.lowTrgEna,
										 info.trigger.he, info.trigger.highTrgEna ),
				       Channel( info.channel.single, 
							 info.channel.all, 
							 info.channel.latc )
				       );
    transferLciCfg( info, lcfg );

    // install the configuration object into the MetaEvent
    lmeta.setConfiguration( lcfg );
  }

  void LSFReader::transferInfo( const eventFile::LSE_Context& ctx, const eventFile::LCI_TKR_Info& info, MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local Configuration object
    LciTkrConfiguration lcfg( info.injected,
				       info.delay,
				       info.threshold,
                                       info.splitLow,
                                       info.splitHigh,
				       Channel( info.channel.single, 
							 info.channel.all, 
							 info.channel.latc )
				       );
    transferLciCfg( info, lcfg );

    // install the configuration object into the MetaEvent
    lmeta.setConfiguration( lcfg );
  }

}
