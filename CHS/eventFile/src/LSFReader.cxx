#include "eventFile/LSFReader.h"
#include "eventFile/LSE_Context.h"
#include "eventFile/EBF_Data.h"
#include "eventFile/LSE_Info.h"

#include "lsfData/LsfCcsds.h"
#include "lsfData/LsfDatagramInfo.h"
#include "lsfData/LsfGemScalers.h"
#include "lsfData/LsfGemTime.h"
#include "lsfData/LsfRunInfo.h"
#include "lsfData/LsfTime.h"
#include "lsfData/LsfTimeTone.h"
#include "lsfData/LsfMetaEvent.h"

namespace eventFile {
  
  bool LSFReader::read( lsfData::LsfCcsds& lccsds, lsfData::MetaEvent& lmeta, EBF_Data& ebf )
  {
    // create LSE objects to hold the retrieved values
    LSE_Context        ctx;
    LSE_Info::InfoType infotype;
    LPA_Info           pinfo;
    LCI_ACD_Info       ainfo;
    LCI_CAL_Info       cinfo;
    LCI_TKR_Info       tinfo;

    // read the native objects
    if ( !LSEReader::read( ctx, ebf, infotype, pinfo, ainfo, cinfo, tinfo ) ) {
      return false;
    }

    // transfer the CCSDS information
    transferCcsds( ctx, lccsds );

    // transfer the context information
    transferContext( ctx, lmeta );

    // transfer the type-specific meta-information
    switch ( infotype ) {
    case LSE_Info::LPA:
      transferInfo( ctx, pinfo, lmeta );
      break;
    case LSE_Info::LCI_ACD:
      transferInfo( ctx, ainfo, lmeta );
      break;
    case LSE_Info::LCI_CAL:
      transferInfo( ctx, cinfo, lmeta );
      break;
    case LSE_Info::LCI_TKR:
      transferInfo( ctx, tinfo, lmeta );
      break;
    }

    return true;
  }

  void LSFReader::transferCcsds( const LSE_Context& ctx, lsfData::LsfCcsds& lccsds )
  {
    lccsds.initialize( ctx.ccsds.scid, ctx.ccsds.apid, ctx.ccsds.utc );
  }
  
  void LSFReader::transferContext( const LSE_Context& ctx, lsfData::MetaEvent& lsfmeta )
  {
    // set the datagram information
    enums::Lsf::Open::Action  ao = static_cast< enums::Lsf::Open::Action  >( ctx.open.action );
    enums::Lsf::Open::Reason  ro = static_cast< enums::Lsf::Open::Reason  >( ctx.open.reason );
    enums::Lsf::Crate         cr = static_cast< enums::Lsf::Crate         >( ctx.open.crate );
    enums::Lsf::Mode          md = static_cast< enums::Lsf::Mode          >( ctx.open.mode );
    enums::Lsf::Close::Action ac = static_cast< enums::Lsf::Close::Action >( ctx.close.action );
    enums::Lsf::Close::Reason rc = static_cast< enums::Lsf::Close::Reason >( ctx.close.reason );
    lsfData::DatagramInfo dgm( ao, ro, cr, md, ac, rc,
			       ctx.open.datagrams, ctx.open.modeChanges 
			       );
    lsfmeta.setDatagram( dgm );

    // set the run information
    enums::Lsf::Platform   pl = static_cast< enums::Lsf::Platform   >( ctx.run.platform );
    enums::Lsf::DataOrigin od = static_cast< enums::Lsf::DataOrigin >( ctx.run.origin );
    lsfData::RunInfo  run( pl, od,
			   ctx.run.groundId, ctx.run.startedAt 
			   );
    lsfmeta.setRun( run );

    // set the GEM scalers
    lsfData::GemScalers sca( ctx.scalers.elapsed, ctx.scalers.livetime,
			  ctx.scalers.prescaled, ctx.scalers.discarded,
			  ctx.scalers.sequence, ctx.scalers.deadzone
			  );
    lsfmeta.setScalers( sca );
  }
  
  void LSFReader::transferTime( const LSE_Context& ctx, const LSE_Info& info, lsfData::MetaEvent& lsfmeta )
  {
    // get the current timetone info
    unsigned char flags = 0;
    flags |= ( ctx.current.missingTimeTone ) ? lsfData::TimeTone::MISSING_TIMETONE_MASK : 0x0;
    flags |= ( ctx.current.missingLatPps )   ? lsfData::TimeTone::MISSING_LAT_MASK      : 0x0;
    flags |= ( ctx.current.missingCpuPps )   ? lsfData::TimeTone::MISSING_CPU_MASK      : 0x0;
    flags |= ( ctx.current.missingGps )      ? lsfData::TimeTone::MISSING_GPS_MASK      : 0x0;
    lsfData::TimeTone curr( ctx.current.incomplete, ctx.current.timeSecs, 
			    ctx.current.flywheeling, flags,
			    lsfData::GemTime( ctx.current.timeHack.hacks, ctx.current.timeHack.tics )
			    );

    // get the previous timetone info
    flags = 0;
    flags |= ( ctx.previous.missingTimeTone ) ? lsfData::TimeTone::MISSING_TIMETONE_MASK : 0x0;
    flags |= ( ctx.previous.missingLatPps )   ? lsfData::TimeTone::MISSING_LAT_MASK      : 0x0;
    flags |= ( ctx.previous.missingCpuPps )   ? lsfData::TimeTone::MISSING_CPU_MASK      : 0x0;
    flags |= ( ctx.previous.missingGps )      ? lsfData::TimeTone::MISSING_GPS_MASK      : 0x0;
    lsfData::TimeTone prev( ctx.previous.incomplete, ctx.previous.timeSecs, 
			    ctx.previous.flywheeling, flags,
			    lsfData::GemTime( ctx.previous.timeHack.hacks, ctx.previous.timeHack.tics )
			    );

    // set the current, previous, and event-time fields
    lsfmeta.setTime( lsfData::Time( curr, prev,
				    lsfData::GemTime( info.timeHack.hacks, info.timeHack.tics ),
				    info.timeTics )
		     );
  }

  void LSFReader::transferLciCfg( const LCI_Info& info, lsfData::LciConfiguration& lcfg )
  {
    unsigned char flags(0);
    flags |= ( info.autoRange      ) ? enums::Lsf::AUTORANGE        : 0x0;
    flags |= ( info.zeroSupression ) ? enums::Lsf::ZERO_SUPPRESSION : 0x0;
    lcfg.set( info.softwareKey, 
	      info.writeCfg,
	      info.readCfg,
	      info.periodicPrescale,
	      flags
	      );
  }

  void LSFReader::transferInfo( const LSE_Context& ctx, const LPA_Info& info, lsfData::MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local lsfData::Configuration object
    lsfData::LpaConfiguration lcfg( info.hardwareKey, info.softwareKey );

    // install the configuration object into the lsfData::MetaEvent
    lmeta.setConfiguration( lcfg );
  }

  void LSFReader::transferInfo( const LSE_Context& ctx, const LCI_ACD_Info& info, lsfData::MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local lsfData::Configuration object
    lsfData::LciAcdConfiguration lcfg( info.injected, 
				       info.threshold, 
				       info.biasDac, 
				       info.holdDelay,
				       lsfData::LciAcdConfiguration::AcdTrigger( info.trigger.veto, 
										 info.trigger.vetoVernier,
										 info.trigger.hld ),
				       lsfData::Channel( info.channel.single, 
							 info.channel.all, 
							 info.channel.latc )
				       );
    transferLciCfg( info, lcfg );

    // install the configuration object into the lsfData::MetaEvent
    lmeta.setConfiguration( lcfg );
  }

  void LSFReader::transferInfo( const LSE_Context& ctx, const LCI_CAL_Info& info, lsfData::MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local lsfData::Configuration object
    lsfData::LciCalConfiguration lcfg( info.uld,
				       info.injected,
				       info.delay,
				       info.threshold,
				       lsfData::LciCalConfiguration::CalTrigger( info.trigger.le,
										 info.trigger.he ),
				       lsfData::Channel( info.channel.single, 
							 info.channel.all, 
							 info.channel.latc )
				       );
    transferLciCfg( info, lcfg );

    // install the configuration object into the lsfData::MetaEvent
    lmeta.setConfiguration( lcfg );
  }

  void LSFReader::transferInfo( const LSE_Context& ctx, const LCI_TKR_Info& info, lsfData::MetaEvent& lmeta )
  {
    // set the timing information
    transferTime( ctx, info, lmeta );

    // create & populate a local lsfData::Configuration object
    lsfData::LciTkrConfiguration lcfg( info.injected,
				       info.delay,
				       info.threshold,
				       lsfData::Channel( info.channel.single, 
							 info.channel.all, 
							 info.channel.latc )
				       );
    transferLciCfg( info, lcfg );

    // install the configuration object into the lsfData::MetaEvent
    lmeta.setConfiguration( lcfg );
  }

}
