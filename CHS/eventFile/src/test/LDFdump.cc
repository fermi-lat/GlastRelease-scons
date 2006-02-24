/*
** ++
**
**  Abstract:
*+      $Id$
**
**      Example program showing how to use the LDF iterators to dump LDF data.
*-
**  Author:
**      R. Claus
**
**  $Revision$
**      February 21, 2003 -- Created
**
**  Copyright:
**                                Copyright 2003
**                                      by
**                         The Board of Trustees of the
**                       Leland Stanford Junior University.
**                              All rights reserved.
**
** --
*/

#include "copyright_SLAC.h"

extern "C" {
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
}

#include "DFC/DFC_endianness.h"
#include "EBF_swap.h"

#include "LATp.h"
#include "eventSummary.h"
#include "errorSummary.h"
#include "AEMheader.h"
#include "ACDpha.h"
#include "CALlog.h"
#include "CALlogEnd.h"
#include "TKRstrip.h"
#include "EBFevent.h"
#include "EBFcontribution.h"
#include "GEMcontribution.h"
#include "GLTcontribution.h"
#include "AEMcontribution.h"
#include "AEMcontributionIterator.h"
#include "CALcontribution.h"
#include "TKRcontribution.h"
#include "CALcontributionIterator.h"
#include "TKRcontributionIterator.h"
#include "DIAGcontributionIterator.h"
#include "ERRcontributionIterator.h"
#include "LATcomponentIterator.h"
#include "EBFeventIterator.h"
#include "LATdatagram.h"
#include "LATdatagramIterator.h"
#include "LATcontributionIterator.h"
#include "GXXXerror.h"
#include "CALdiagnostic.h"
#include "OSWcontribution.h"
#include "OSWtimeContribution.h"
#include "OSWcontributionIterator.h"
#include "OSWtimeBase.h"
//#include "ASCtileContributionIterator.h"
//#include "ASCcontributionIterator.h"

class MyLATPcellHeader : public LATPcellHeader
{
public:
  MyLATPcellHeader()  {}
  ~MyLATPcellHeader() {}

  static void dump(const unsigned header, const char* prefix = "");
};

void MyLATPcellHeader::dump(const unsigned header, const char* pfx)
{
  printf("%s     source = 0x%02x\n", pfx, LATPcellHeader::source     (header));
  printf("%sdestination = 0x%02x\n", pfx, LATPcellHeader::destination(header));
  printf("%s    respond = %d\n"    , pfx, LATPcellHeader::respond    (header));
  printf("%s   protocol = %d\n"    , pfx, LATPcellHeader::protocol   (header));
  printf("%s     parity = %d\n"    , pfx, LATPcellHeader::parity     (header));
}

class MyEventSummary : public EventSummary
{
public:
  MyEventSummary()  {}
  ~MyEventSummary() {}

  static void dump(const unsigned summary, const char* prefix = "");
};

void MyEventSummary::dump(const unsigned summary, const char* pfx)
{
  unsigned evtNo = EventSummary::eventNumber(summary);
  unsigned tag   = EventSummary::tag        (summary);
  unsigned seq   = (evtNo << 2) | tag;

  printf("%s eventSequence = %d = 0x%08x = ((eventNumber << 2) | tag)\n",
         pfx, seq, seq);
  printf("%s   eventNumber = %d\n", pfx, evtNo);
  printf("%s           tag = %d\n", pfx, tag);
  printf("%s     calStrobe = %d\n", pfx, EventSummary::calStrobe     (summary));
  printf("%s          TACK = %d\n", pfx, EventSummary::TACK          (summary));
  printf("%s      readout4 = %d\n", pfx, EventSummary::readout4      (summary));
  printf("%s  zeroSuppress = %d\n", pfx, EventSummary::zeroSuppress  (summary));
  printf("%s        marker = %d\n", pfx, EventSummary::marker        (summary));
  printf("%s         error = %d\n", pfx, EventSummary::error         (summary));
  printf("%s    diagnostic = %d\n", pfx, EventSummary::diagnostic    (summary));
  printf("%strgParityError = %d\n", pfx, EventSummary::trgParityError(summary));
}

class MyEBFevent : public EBFevent
{
public:
  MyEBFevent()  {}
  ~MyEBFevent() {}

  void dump(const char* prefix);
};

void MyEBFevent::dump(const char* prefix)
{
  printf("%sEBF format identity     = 0x%08x\n",
         prefix, identity().value());
  printf("%sEvent status            = 0x%04x = %d\n",
         prefix, status(), status());
  printf("%sEvent length            = 0x%04x = %d Bytes\n",
         prefix, length(), length());
  printf("%sEvent summary           = 0x%08x\n",
         prefix, summary());
  printf("\n");
  //  printf("%x",this);
}


class MyEBFcontribution : public EBFcontribution
{
public:
  MyEBFcontribution()  {}
  ~MyEBFcontribution() {}

  void dumpCommon(EBFevent* event, const char* prefix = "");
  void dump      (EBFevent* event, const char* prefix = "");
};

void MyEBFcontribution::dumpCommon(EBFevent* event, const char* prefix)
{
  char pfx[80];
  //  printf("%sContribution base         = %x \n",    prefix, this);
  printf("%sContribution length     = %d Bytes\n", prefix, length());
  printf("%sContribution error      = 0x%1x\n",    prefix, packetError());
  printf("%sContribution sequence   = %d\n",       prefix, sequence());
  printf("%sLCB Header              = 0x%08x\n",   prefix, header());
  sprintf(pfx, "%s            ", prefix);
  MyLATPcellHeader::dump(header(), pfx);
  printf("%sSummary                 = 0x%08x\n",   prefix, summary());
  sprintf(pfx, "%s         ", prefix);
  MyEventSummary::dump(summary(), pfx);
}

void MyEBFcontribution::dump(EBFevent* event, const char* prefix)
{
  dumpCommon(event, prefix);

  unsigned* buffer = (unsigned*)payload();
  unsigned  count  = (length() >> 2) - 2;
  unsigned  rows   = count >> 2;
  printf("%sData:\n", prefix);
  while (rows--)
  {
    printf("%s  %08x, %08x, %08x, %08x\n", prefix, buffer[0], buffer[1],
                                                   buffer[2], buffer[3]);
    buffer += 4;
  }

  printf("%s  ", prefix);
  unsigned l = (count & 3) - 1;
  while (l--)
  {
    printf("%08x, ", *buffer++);
  }
  printf("%08x\n", *buffer);
}


class MyGLTcontribution : public GLTcontribution
{
public:
  MyGLTcontribution()  {}
  ~MyGLTcontribution() {}

  void dump(EBFevent* event, const char* prefix = "");
};

void MyGLTcontribution::dump(EBFevent* /*event*/, const char* prefix)
{
  const struct timespec* ts = timeStamp();
  const OSWtimeBase*     tb = timebase();

  printf("%sGLT:\n", prefix);
  printf("%s  Event GMT timestamp   = %ld.%09ld seconds after 1/1/1970\n",
         prefix, ts->tv_sec, ts->tv_nsec);
  if (payloadSize() > 8)                // Timebase does not exist in old data
    printf("%s  PPC timebase          = 0x%08x%08x\n", prefix, tb->upper(), tb->lower());
}


class MyGEMcontribution : public GEMcontribution
{
public:
  MyGEMcontribution()  {}
  ~MyGEMcontribution() {}

  void dump(EBFevent* event, const char* prefix = "");
};

void MyGEMcontribution::dump(EBFevent* /*event*/, const char* prefix)
{
  const GEMtileList* tl  = tileList();
  GEMonePPStime      opt = onePPStime();

  printf("%sGEM:\n", prefix);
  printf("%s  ROI vector        = 0x%04x\n", prefix, roiVector());
  printf("%s  TKR vector        = 0x%04x\n", prefix, tkrVector());
  printf("%s  CAL HE vector     = 0x%04x\n", prefix, calHEvector());
  printf("%s  CAL LE vector     = 0x%04x\n", prefix, calLEvector());
  printf("%s  Condition summary = 0x%02x\n", prefix, conditionSummary());
  printf("%s  Missed (deadZone) = 0x%02x\n", prefix, missed()); // or deadZone()
  printf("%s  CNO vector        = 0x%04x\n", prefix, cnoVector());
  printf("%s  tile list:\n", prefix);
  printf("%s    XZM, XZP        = 0x%04x, 0x%04x\n", prefix, tl->XZM(), tl->XZP());
  printf("%s    YZM, YZP        = 0x%04x, 0x%04x\n", prefix, tl->YZM(), tl->YZP());
  printf("%s    XY              = 0x%08x\n",       prefix, tl->XY());
  printf("%s    RBN, NA         = 0x%04x, 0x%04x\n", prefix, tl->RBN(), tl->NA());
  printf("%s  Live time         = 0x%08x = %d\n", prefix, liveTime(), liveTime());
  printf("%s  Prescaled         = 0x%08x = %d\n", prefix, prescaled(), prescaled());
  printf("%s  Discarded         = 0x%08x = %d\n", prefix, discarded(), discarded());
  // printf("%s  Sent              = 0x%08x = %d\n", prefix, sent(), sent());
  printf("%s  Condition arrival:\n", prefix);
  printf("%s    raw             = 0x%08x     \n", prefix, condArrTime().datum());
  printf("%s    external        = 0x%02x = %d\n", prefix, condArrTime().external(), condArrTime().external() );
  printf("%s    cno             = 0x%02x = %d\n", prefix, condArrTime().cno()     , condArrTime().cno()      );
  printf("%s    calHE           = 0x%02x = %d\n", prefix, condArrTime().calHE()   , condArrTime().calHE()    );
  printf("%s    calLE           = 0x%02x = %d\n", prefix, condArrTime().calLE()   , condArrTime().calLE()    );
  printf("%s    tkr             = 0x%02x = %d\n", prefix, condArrTime().tkr()     , condArrTime().tkr()      );
  printf("%s    roi             = 0x%02x = %d\n", prefix, condArrTime().roi()     , condArrTime().roi()      );
  printf("%s  Trigger time      = 0x%08x = %d\n", prefix, triggerTime(), triggerTime());
  printf("%s  One PPS time:\n", prefix);
  printf("%s    seconds, timebase = 0x%02x, 0x%07x = %d, %d\n",
         prefix, opt.seconds(), opt.timebase(),
         opt.seconds(), opt.timebase());
  printf("%s  Delta window open time  = 0x%04x = %d\n", prefix, deltaWindowOpenTime(), deltaWindowOpenTime());
  printf("%s  Delta event time  = 0x%04x = %d\n", prefix, deltaEventTime(), deltaEventTime());
}


class MyOSWiterator : public OSWcontributionIterator
{
public:
  MyOSWiterator(EBFevent*        event,
                OSWcontribution* contribution,
                const char*      prefix) :
    OSWcontributionIterator(event, contribution),
    _prefix(prefix) {};
  virtual ~MyOSWiterator() {}

  virtual int OSW_time(EBFevent*, OSWtimeContribution*) const;
  virtual int OSW_UDF(EBFevent*, OSWcontribution*) const;
  virtual int handleError(OSWcontribution* contribution, unsigned code, unsigned p1=0, unsigned p2=0) const;

  void dump();

private:
  const char* _prefix;
};

void MyOSWiterator::dump()
{
  iterate();
}

int MyOSWiterator::OSW_time(EBFevent*            event,
                            OSWtimeContribution* contribution) const
{
  const struct timespec* ts = contribution->timeStamp();
  const OSWtimeBase*     tb = contribution->timebase();
  printf("%s  OSWtime:\n", _prefix);
  printf("%s    Event GMT timestamp = %ld.%09ld seconds after 1/1/1970\n",
         _prefix, ts->tv_sec, ts->tv_nsec);
  printf("%s    PPC timebase        = 0x%08x%08x\n", _prefix,
         tb->upper(), tb->lower());
  printf("%s    Extended event no.  = 0x%08x = %d\n", _prefix,
         contribution->evtSequence(), contribution->evtSequence());

  return 0;
}

int MyOSWiterator::OSW_UDF(EBFevent*        /* event */,
                           OSWcontribution* contribution) const
{
  printf("%sUndefined OSW contribution encountered with identity %08x\n",
         _prefix, contribution->identity().value());
  return 0;
}

int MyOSWiterator::handleError(OSWcontribution* contribution,
                               unsigned code, unsigned p1, unsigned p2) const
{
  fprintf(stderr, "%sMyOSWiterator::handleError:  Somehow an error occured. \n",
         _prefix);
  fprintf(stderr, "%s  code=%d, p1=%d, p2=%d\n", _prefix, code, p1, p2);
  return 0;
}

class MyACDiterator : public AEMcontributionIterator
{
public:
  MyACDiterator(EBFevent*        event,
                AEMcontribution* contribution,
                const char*      prefix);
  virtual ~MyACDiterator() {};
  void dump();

  virtual void header(unsigned cable, AEMheader hdr);
  virtual void pha(unsigned cable, unsigned channel, ACDpha p);
  virtual int handleError(AEMcontribution* contribution, unsigned code, unsigned p1=0, unsigned p2=0) const;
private:
  const char* _prefix;
};

MyACDiterator::MyACDiterator(EBFevent*        event,
                             AEMcontribution* contribution,
                             const char*      prefix):
  AEMcontributionIterator(event, contribution),
  _prefix(prefix)
{
}

void MyACDiterator::dump()
{
  iterate();
}

void MyACDiterator::header(unsigned cable, AEMheader hdr)
{
  printf("%s  Cable %d = FREE board %s header:\n", _prefix, cable,
         map()->freeName(event()->identity(), cable));
  printf("%s    Start-bit           = %01x\n",     _prefix, hdr.startBit());
  printf("%s    Hit-map             = 0x%05x\n",   _prefix, hdr.hitMap());
  printf("%s    Accept-map          = 0x%05x\n",   _prefix, hdr.acceptMap());
  printf("%s    PHA vector          = %01x\n",     _prefix, hdr.phaVector());
  printf("%s    Header parity error = %01x\n",     _prefix, hdr.parityError());
  printf("%s    End of Cables       = %01x\n",     _prefix, hdr.endOfCables());
  printf("%s    Cable Number        = %02x\n",     _prefix, hdr.cableNumber());
  if (hdr.phaVector())
  {
    printf("%s                                 Value     Parity\n",      _prefix);
    printf("%s    Channel  Tile  Side Range (dec)  (hex) error  More\n", _prefix);
  }
  else
    printf("%s  Cable %d has no data\n", _prefix, cable);
}

void MyACDiterator::pha(unsigned cable, unsigned channel, ACDpha p)
{
  LATtypeId          id   = event()->identity();
  const ACDtileSide* tbl  = map()->lookup(id, cable, channel);
  char               side = tbl->a() ? 'A' : 'B';
  printf("%s       %2d   %4s     %c    %d   %4d = 0x%03x    %d     %d\n",
         _prefix, channel, tbl->name(), side,
         p.ADCrange(), p.ADCvalue(), p.ADCvalue(), p.parityError(), p.more());
}

int MyACDiterator::handleError(AEMcontribution* contribution,
                               unsigned code, unsigned p1, unsigned p2) const
{
  switch (code)
  {
    case AEMcontributionIterator::ERR_TooManyPhas:
    {
      fprintf(stderr, "%sAEMcontributionIterator::iterate: %s", _prefix,
                      "more PHA values found than bits set in the Accept Map\n");
      break;
    }

    default: break;
  }
  return 0;
}

class MyCALiterator : public CALcontributionIterator
{
public:
  MyCALiterator(EBFevent*        event,
                CALcontribution* contribution,
                const char*      prefix);
  virtual ~MyCALiterator() {}

  virtual void log(unsigned tower, unsigned layer, CALlog theLog);
  virtual int  handleError(CALcontribution*, unsigned code, unsigned p1=0, unsigned p2=0) const;
  void         dump();

private:
  const char* _prefix;
};

MyCALiterator::MyCALiterator(EBFevent*        event,
			     CALcontribution* contribution,
			     const char*      prefix) :
  CALcontributionIterator(event, contribution),
  _prefix(prefix)
{
}

void MyCALiterator::dump()
{
  // Dump the data
  unsigned length = contribution()->numLogAccepts();
  if (length == 0)
  {
    printf("%sCAL: no data - empty contribution\n", _prefix);
    //CALend(3*sizeof(unsigned));
    return;
  }

  printf("%sCAL:\n", _prefix);
  printf("%s  Log Accepts           = 0x%08x\n", _prefix, contribution()->logAccepts());
  printf("%s  Number of Log Accepts = %d\n", _prefix,     contribution()->numLogAccepts());

  printf("%s  Data:\n", _prefix);
  printf("%s                               +ADC      +ADC     -ADC      -ADC\n",  _prefix);
  printf("%s   tower GCCC  layer column    value     range    value     range\n", _prefix);

  iterate();
}

void MyCALiterator::log(unsigned tower, unsigned layer, CALlog theLog)
{
  // dump out one log's worth of data
  const char* layerTag[] = {"x0", "y0", "x1", "y1", "x2", "y2", "x3", "y3"};
  const unsigned  gccc[] = {  0,    3,    0,    3,    2,    1,    2,    1};

  char name[10];
  sprintf(name, "  %1d    %2s ", gccc[layer], layerTag[layer]);

  printf("%s     %2d   %10s   %2d   %4d = 0x%03x   %1d   %4d = 0x%03x   %1d\n",
	 _prefix, tower, name, theLog.column(),
	 theLog.positive().value(), theLog.positive().value(), theLog.positive().range(),
	 theLog.negative().value(), theLog.negative().value(), theLog.negative().range());
}

int MyCALiterator::handleError(CALcontribution* contribution,
                               unsigned code, unsigned p1, unsigned p2) const
{
  fprintf(stderr, "%sMyCALiterator::handleError:  Somehow an error occured. \n",
         _prefix);
  fprintf(stderr, "%s  code=%d, p1=%d, p2=%d\n", _prefix, code, p1, p2);
  return 0;
}


class MyTKRiterator : public TKRcontributionIterator
{
public:
  MyTKRiterator(EBFevent*        event,
                TKRcontribution* contribution,
                const char*      prefix);
  virtual ~MyTKRiterator() {}

  virtual void strip(unsigned tower, unsigned layerEnd, unsigned hit);
  virtual void TOT  (unsigned tower, unsigned layerEnd, unsigned hit);
  virtual int handleError(TKRcontribution*, unsigned code, unsigned p1=0, unsigned p2=0) const;
  void         dump ();
private:
  const char* _prefix;
};


MyTKRiterator::MyTKRiterator(EBFevent*        event,
                             TKRcontribution* contribution,
                             const char*      prefix) :
  TKRcontributionIterator(event, contribution),
  _prefix(prefix)
{
}

void MyTKRiterator::strip(unsigned tower, unsigned layerEnd, unsigned hit)
{
  const char  xy[]  = {'Y', 'X', 'X', 'Y'};
  unsigned    layer = layerEnd >> 1;    // 0 <= layer < 36
  unsigned    l     = layer >> 1;       // 0 <= l     < 18
  const char* sp    = l < 10 ? " " : "";
  printf("%s      %2d    %s%c%d          %4d = 0x%03x\n",
         _prefix, tower, sp, xy[layer & 0x3], l, hit, hit);
}

void MyTKRiterator::TOT(unsigned tower, unsigned layerEnd, unsigned tot)
{
  const char  xy[]  = {'Y', 'X', 'X', 'Y'};
  unsigned    side  = layerEnd & 0x1;   // 0 = low, 1 = high
  unsigned    layer = layerEnd >> 1;    // 0 <= layer < 36
  unsigned    l     = layer >> 1;       // 0 <= l     < 18
  const char* sp    = l < 10 ? " " : "";
  printf("%s      %2d    %s%c%d     %c    %3d  = 0x%03x\n",
         _prefix, tower, sp, xy[layer & 0x3], l, side ? 'H' : 'L', tot, tot);
}

void MyTKRiterator::dump()
{
  // Dump the data
  const EBFevent* theEvent = event();
  unsigned length  = contribution()->numAccepts(theEvent);
  if (!length)
  {
    printf("%sTKR: no data - empty contribution\n", _prefix);
    //    unsigned* acc = (unsigned*)contribution()->accepts(theEvent);
    //    printf("%s     endTots: 0x%08x\n",_prefix,&acc[3]);
    endTots((unsigned)(contribution()->accepts(theEvent) + 3*sizeof(unsigned)
		       - (char*)contribution()) );
    return;
  }
  printf("%sTKR:\n", _prefix);
  unsigned       array[8];              // 8 longwords of space
  unsigned*      accepts = contribution()->acceptsMask(theEvent, array);
  const unsigned gtcc[]  = {6, 3, 7, 2, 5, 0, 4, 1};
  const char*    label[] = {"Left", "Right", "Right", "Left",
                            "Right", "Left", "Left", "Right"};
  const char*    xy[]    = {"x", "x", "x", "x", "y", "y", "y", "y"};
  printf("%s           GTCC   odd even  GTCC\n", _prefix);

  unsigned i;
  for (i = 0; i < 8; i += 2)
  {
    printf("%s    %s %5s  %d  0x%03x 0x%03x  %d  %s %-5s\n",
           _prefix, xy[i],        label[i],  gtcc[i], accepts[i],
           accepts[i+1], gtcc[i+1], xy[i+1], label[i+1]);
  }

  printf("%s  Number of Accepts     = %d\n",
	 _prefix, contribution()->numAccepts(theEvent));
  printf("%s  Data:\n", _prefix);

  // Iterate over strip data
  char pfx[80];
  sprintf(pfx, "%s    ", _prefix);
  //  MyTKRiterator data(event, this, pfx);
  printf("%sTower  Layer        strip = strip\n", pfx);
  iterateStrips();

  // Iterate over TOT data
  printf("\n%sTower  Layer  Side   TOT  =  TOT\n", pfx);
  iterateTOTs();

}

int MyTKRiterator::handleError(TKRcontribution* contribution,
                               unsigned code, unsigned p1, unsigned p2) const
{
  switch (code)
  {
    case TKRcontributionIterator::ERR_WrongOrder:
    {
      fprintf(stderr, "%sTKRiterator.iterateTOTs: TOTs can not be accessed before"
                      "TKRiterator.iterateStrips has executed.\n", _prefix);
      return 0;
      break;
    }
    default: break;
  }
  return 0;
}

class MyDIAGiterator: public DIAGcontributionIterator
{
public:
  MyDIAGiterator(EBFevent*        event,
                 TEMcontribution* contribution,
                 unsigned         dataStart,
                 const char*      prefix);
  virtual ~MyDIAGiterator() {};
  virtual int CALdiag(unsigned tower, unsigned layer, CALdiagnostic diag);
  virtual int TKRdiag(unsigned tower, unsigned gtcc,  TKRdiagnostic diag);
  virtual int     handleError(TEMcontribution*, unsigned code, unsigned p1=0, unsigned p2=0) const;
private:
  const char* _prefix;
};

MyDIAGiterator::MyDIAGiterator(EBFevent*        event,
                               TEMcontribution* contribution,
                               unsigned         dataStart,
                               const char*      prefix) :
  DIAGcontributionIterator(event, contribution),
  _prefix(prefix)
{
  // Set the start of the diagnostic data
  offset(dataStart);
}

int MyDIAGiterator::handleError(TEMcontribution* contribution,
                                unsigned code, unsigned p1, unsigned p2) const
{
  fprintf(stderr, "%sMyDIAGiterator::handleError:  Somehow an error occured. \n",
         _prefix);
  fprintf(stderr, "%s  code=%d, p1=%d, p2=%d\n", _prefix, code, p1, p2);
  return 0;
}

int MyDIAGiterator::CALdiag(unsigned      /* tower */,
                            unsigned      layer,
                            CALdiagnostic diag)
{
  // Diagnostic cont structure is:
  // - 8 32 bit contributions for GCCC
  //     Order: 0(x0,x1) 2(x2,x3) 3(y0,y1) 4(y2,y3)
  const char*    layerLabel[] = {"x0", "x1", "x2", "x3", "y0", "y1", "y2", "y3"};
  const unsigned gccc[]       = {0, 2, 3, 1};

  printf("%s    %1d    %2s    0x%03x     %1d   %1d     0x%03x     %1d   %1d\n",
         _prefix, gccc[layer>>1], layerLabel[layer],
         diag.logAccepts (1),        // negative
         diag.high       (1),        // negative
         diag.low        (1),        // negative
         diag.logAccepts (0),        // positive
         diag.high       (0),        // positive
         diag.low        (0));       // positive
  return 0;
}

int MyDIAGiterator::TKRdiag(unsigned      /* tower */,
                            unsigned      gtcc,
                            TKRdiagnostic diag)
{
  // - 8 16 bit contributions for GTCC
  printf("%s    %1d   0x%03x\n", _prefix, gtcc, diag.datum() & 0xffff);
  return 0;
}

class MyERRiterator: public ERRcontributionIterator
{
public:
  MyERRiterator( EBFevent*      event,
                 TEMcontribution* contribution,
                 unsigned         dataStart,
                 const char*      prefix);
  virtual ~MyERRiterator() {}
  virtual int     handleError(TEMcontribution* contribution, unsigned code, unsigned p1=0, unsigned p2=0) const;
  // definable error handlers
  virtual int gcccError      (unsigned tower, unsigned gccc, GCCCerror err);
  virtual int gtccError      (unsigned tower, unsigned gtcc, GTCCerror err);
  virtual int phaseError     (unsigned tower, unsigned short err);
  virtual int timeoutError   (unsigned tower, unsigned short err);
  // TKR errors
  virtual int gtrcPhaseError (unsigned tower, unsigned gtcc, unsigned gtrc, GTRCerror err);
  virtual int gtfePhaseError (unsigned tower, unsigned gtcc, unsigned gtrc,
                               unsigned short err1, unsigned short err2,
                               unsigned short err3, unsigned short err4, unsigned short err5);
  virtual int gtccFIFOerror  (unsigned tower, unsigned gtcc, unsigned short err);
  virtual int gtccTMOerror        (unsigned tower, unsigned gtcc);
  virtual int gtccHDRParityError  (unsigned tower, unsigned gtcc);
  virtual int gtccWCParityError   (unsigned tower, unsigned gtcc);
  virtual int gtrcSummaryError    (unsigned tower, unsigned gtcc);
  virtual int gtccDataParityError (unsigned tower, unsigned gtcc);
private:
  const char* _prefix;
};

MyERRiterator::MyERRiterator(EBFevent* event, TEMcontribution* contribution,
                             unsigned dataStart, const char* prefix) :
  ERRcontributionIterator(event,contribution),
  _prefix(prefix)
{
  offset(dataStart);
}

int MyERRiterator::handleError(TEMcontribution* contrib,
                               unsigned code, unsigned p1, unsigned p2) const
{
  fprintf(stderr, "%sMyERRiterator::handleError:  Somehow an error occured. \n",
         _prefix);
  fprintf(stderr, "%s  code=%d, p1=%d, p2=%d\n", _prefix, code, p1, p2);
  return 0;
}

int MyERRiterator::gcccError(unsigned /* tower */, unsigned gccc, GCCCerror err)
{
  if (err.type() == 0)
    printf("%s    GCCC %d had missing start bit(s) on sweep %d (raw: 0x%04x):\n",
           _prefix, gccc, err.sweep(), err.param());
  else
    printf("%s    GCCC %d had parity error(s) on sweep %d (raw: 0x%04x):\n",
           _prefix, gccc, err.sweep(), err.param());
  return 0;
}

int MyERRiterator::gtccError(unsigned /* tower */, unsigned gtcc, GTCCerror err)
{
  printf("%s  GTCC %d, GTRC %d (msg: 0x%04x):\n",
         _prefix, gtcc, err.GTRC(), err.raw());
  return 0;
}

int MyERRiterator::phaseError(unsigned /* tower */, unsigned short err)
{
  printf("\n");
  printf("%s  Cable Controller phasing error tags (raw: 0x%04x):\n",
         _prefix, err);

  unsigned cc[]  = {3, 6, 2, 7, 0, 5, 1, 4};
  unsigned tag[] = {0, 0, 0, 0, 0, 0, 0, 0};
  unsigned i;
  for (i=0; i<8; i++)
  {
    tag[cc[i]] = (err >> (i << 1)) & 0x0003;
  }
  printf("%s    GTCC:  0  1  2  3  4  5  6  7\n" , _prefix);
  printf("%s     tag:", _prefix);
  for (i=0; i<8; i++)
  {
    printf("  %d", tag[i]);
  }
  printf("\n");
  return 0;
}

int MyERRiterator::timeoutError(unsigned /* tower */, unsigned short err)
{
  printf("\n");
  printf("%s  Cable Controller timeouts (raw: 0x%03x):\n", _prefix, err);
  unsigned short cTmos = err & 0x000f;
  printf("%s    GCCC:", _prefix);
  unsigned i;
  for (i=0; i<4; i++)
  {
    if (cTmos & (1 << i))  printf(" %d", i);
  }
  unsigned short tTmos = (err >> 4) & 0x00ff;
  printf(", GTCC:");
  for (i=0; i<8; i++)
  {
    if (tTmos & (1 << i))  printf(" %d", i);
  }
  printf("\n");
  return 0;
}

int MyERRiterator::gtrcPhaseError (unsigned  /* tower */,
                                   unsigned  /* gtcc */,
                                   unsigned  /* gtrc */,
                                   GTRCerror err)
{
  printf("\n");
  printf("%s    GTRC phasing error (raw: 0x%04x)\n", _prefix, err.raw());
  printf("%s      Current tag number:           %d",
         _prefix, err.currentTag());
  printf("%s      Layer 0 tag number:           %d\n",
         _prefix, err.layer0Tag());
  printf("%s      Expected layer number:        %d\n",
         _prefix, err.expectedLayer());
  printf("%s      Layer number received by TEM: %d\n",
         _prefix, err.receivedLayer());
  return 0;
}

int MyERRiterator::gtfePhaseError (unsigned       /* tower */,
                                   unsigned       /* gtcc */,
                                   unsigned       /* gtrc */,
                                   unsigned short err1,
                                   unsigned short err2,
                                   unsigned short err3,
                                   unsigned short err4,
                                   unsigned short err5)
{
  printf("\n");
  printf("%s    GTFE phasing error, GTFE tags: ", _prefix);
  printf(" 0x%03x 0x%03x 0x%03x 0x%03x 0x%03x", err1, err2, err3, err4, err5);
  printf("\n");
  return 0;
}

int MyERRiterator::gtccFIFOerror  (unsigned       /* tower */,
                                   unsigned       /* gtcc */,
                                   unsigned short err)
{
  printf("\n");
  printf("%s    FIFO full error, projected word count = 0x%04x\n",
         _prefix, err);
  return 0;
}

int MyERRiterator::gtccTMOerror        (unsigned /* tower */,
                                        unsigned /* gtcc */)
{
  printf("%s    Cable timeout error\n", _prefix);
  return 0;
}

int MyERRiterator::gtccHDRParityError  (unsigned /* tower */,
                                        unsigned /* gtcc */)
{
  printf("%s    Header parity error\n", _prefix);
  return 0;
}

int MyERRiterator::gtccWCParityError   (unsigned /* tower */,
                                        unsigned /* gtcc */)
{
  printf("%s    Word count parity error\n", _prefix);
  return 0;
}

int MyERRiterator::gtrcSummaryError    (unsigned /* tower */,
                                        unsigned /* gtcc */)
{
  printf("%s    GTRC summary error\n", _prefix);
  return 0;
}

int MyERRiterator::gtccDataParityError (unsigned /* tower */,
                                        unsigned /* gtcc */)
{
  printf("%s    Data parity error\n", _prefix);
  return 0;
}

class MyLATcomponentIterator : public LATcomponentIterator
{
public:
  MyLATcomponentIterator() :
    LATcomponentIterator(),
    _calSrc(0xffffffff) {}
  virtual ~MyLATcomponentIterator() {}

  virtual int handleError(EBFcontribution* contribution, unsigned code, unsigned p1=0, unsigned p2=0) const;
  virtual int handleError(EBFevent* event, unsigned code, unsigned p1=0, unsigned p2=0) const;

  virtual int UDFcomponent(EBFevent* event, EBFcontribution* contribution);
  virtual int OSWcomponent(EBFevent* event, OSWcontribution* contribution);
  virtual int GLTcomponent(EBFevent* event, GLTcontribution* contribution);
  virtual int GEMcomponent(EBFevent* event, GEMcontribution* contribution);
  virtual int ACDcomponent(EBFevent* event, AEMcontribution* contribution);
  virtual int CALcomponent(EBFevent* event, CALcontribution* contribution);
  virtual int TKRcomponent(EBFevent* event, TKRcontribution* contribution);
  virtual int diagnostic  (EBFevent* event, TEMcontribution* contribution);
  virtual int error       (EBFevent* event, TEMcontribution* contribution);
  virtual int cleanup     (EBFevent* event, TEMcontribution* contribution);

private:
  unsigned _calSrc;
};

int MyLATcomponentIterator::handleError(EBFcontribution* contrib,
                                        unsigned code, unsigned p1, unsigned p2) const
{
  fprintf(stderr, "MyLATcomponentIterator::handleError:  Somehow an error occured. \n");
  fprintf(stderr, "  code=%d, p1=%d, p2=%d\n", code, p1, p2);
  return 0;
}

int MyLATcomponentIterator::handleError(EBFevent* event,
                                        unsigned code, unsigned p1, unsigned p2) const
{
  switch (code)
  {
    case EBFcontributionIterator::ERR_NumContributions:
    {
      fprintf(stderr, "EBFcontributionIterator::iterate: "
                     "Number of contributions found > %d\n",
                     p1);
      break;
    }
    case EBFcontributionIterator::ERR_PastEnd:
    {
      fprintf(stderr, "EBFcontributionIterator::iterate: "
                     "Iterated past end of event by 0x%0x = %d bytes\n",
                     p1, p1);
      break;
    }
    case EBFcontributionIterator::ERR_ZeroLength:
    {
      fprintf(stderr, "\nEBFcontributionIterator::iterate:\n"
             "  Found a contribution with zero length\n");
      break;
    }
    case EBFcontributionIterator::ERR_PacketError:
    {
      /* ERR_PacketError gets an additional parameter */
      char *type;

      switch (p1)
      {
        case EBFcontribution::Parity:      type = "parity";       break;
        case EBFcontribution::Truncated:   type = "truncated";    break;
        case EBFcontribution::WriteFault:  type = "write fault";  break;
        case EBFcontribution::TimedOut:    type = "timed out";    break;
        default:                           type = "unknown";      break;
      }
      fprintf(stderr, "\nEBFcontributionIterator::iterate:\n"
          "  Skipping contribution with source ID %d having packet %s error\n",
             p2, type);
      break;
    }
    case EBFcontributionIterator::ERR_NoMap:
    {
      fprintf(stderr, "No contribution map exists for EBF version %0x\n",
             p1);
      break;
    }
    default: break;
  }
  return code;
}

int MyLATcomponentIterator::UDFcomponent(EBFevent*        event,
                                         EBFcontribution* contribution)
{
  printf ("\nUndefined EBF component\n");
  ((MyEBFcontribution*)contribution)->dump(event, "  ");
  return 0;
}

int MyLATcomponentIterator::OSWcomponent(EBFevent*        event,
                                         OSWcontribution* contribution)
{
  const char* prefix = "  ";
  printf("\nOSW:\n");
  ((MyEBFcontribution*)contribution)->dump(event, prefix);
  printf("%sOSW:\n", prefix);
  MyOSWiterator iterator(event, contribution, prefix);
  iterator.dump();
  return 0;
}

int MyLATcomponentIterator::GLTcomponent(EBFevent*        event,
                                         GLTcontribution* contribution)
{
  printf("\nGLT:\n");
  ((MyEBFcontribution*)contribution)->dump(event, "  ");
  ((MyGLTcontribution*)contribution)->dump(event, "  ");
  return 0;
}

int MyLATcomponentIterator::GEMcomponent(EBFevent*        event,
                                         GEMcontribution* contribution)
{
  printf("\nGEM:\n");
  ((MyEBFcontribution*)contribution)->dump(event, "  ");
  ((MyGEMcontribution*)contribution)->dump(event, "  ");
  return 0;
}

int MyLATcomponentIterator::ACDcomponent(EBFevent*        event,
                                         AEMcontribution* contribution)
{
  printf("\nAEM:\n");
  const char* prefix = "  ";
  ((MyEBFcontribution*)contribution)->dumpCommon(event, prefix);
  printf("%sACD:\n",prefix);
  MyACDiterator iterator(event, contribution, prefix);
  iterator.dump();
  return 0;
}

int MyLATcomponentIterator::CALcomponent(EBFevent*        event,
                                         CALcontribution* contribution)
{
  const char* prefix = "  ";
  _calSrc = LATPcellHeader::source(contribution->header());
  printf("\nTEM %2d:\n", _calSrc);
  ((MyEBFcontribution*)contribution)->dumpCommon(event, prefix);
  MyCALiterator iterator(event,contribution,prefix);
  iterator.dump();
  CALend(iterator.CALend());
  return 0;
}

int MyLATcomponentIterator::TKRcomponent(EBFevent*        event,
                                         TKRcontribution* contribution)
{
  const char* prefix = "  ";
  if (_calSrc != LATPcellHeader::source(contribution->header()))
  {
    printf("\nTEM: %d\n", LATPcellHeader::source(contribution->header()));
    ((MyEBFcontribution*)contribution)->dumpCommon(event, prefix);
  }
  MyTKRiterator iterator(event,contribution,prefix);
  iterator.dump();
  TKRend(iterator.diagnostic());
  return 0;
}

int MyLATcomponentIterator::diagnostic (EBFevent*        event,
                                        TEMcontribution* contribution)
{
  const char* prefix = "  ";

  if ( EventSummary::diagnostic(contribution->summary()))
  {
    MyDIAGiterator iter(event,contribution,TKRend(),prefix);

    printf( "\n");
    printf( "%sCAL diagnostic data:\n",prefix);
    printf( "%s             Negative Layer-End  Positive Layer-End\n",prefix);
    printf( "%s  GCCC layer logAccept high low  logAccept high low\n",prefix);
    iter.iterateCAL();

    printf( "\n");
    printf( "%sTKR diagnostic data:\n", prefix);
    printf( "%s       Trigger\n", prefix);
    printf( "%s       Request\n", prefix);
    printf( "%s  GTCC  (GTRC)\n", prefix);
    iter.iterateTKR();

    diagnosticEnd(TKRend()+iter.size());  // Reset the end of the diagnostic data
  }
  return 0;
}

int MyLATcomponentIterator::error (EBFevent*        event,
                                   TEMcontribution* contribution)
{
  // parse out errors, and determine if the event ends normally

  const char* prefix = "  ";

  if ( EventSummary::error(contribution->summary()))
  {
    unsigned offset;
    if (0 != diagnosticEnd())
      offset=diagnosticEnd();
    else
      offset=TKRend();
    MyERRiterator iter(event,contribution,offset,prefix);
    ErrorSummary theError = iter.theError();

    printf("\n");
    printf("%sError data:\n",prefix);
    printf("%s   Cable    Cable  \n",prefix);
    printf("%s  timeout  phasing  TKR  CAL\n",prefix);
    printf("%s     %1d        %1d    0x%02x  0x%01x\n",
           prefix,
           theError.tmo(),
           theError.phs(),
           theError.tkr(),
           theError.cal());

    iter.iterate();

    errorEnd(offset+iter.size());
  }
  return 0;
}

int MyLATcomponentIterator::cleanup (EBFevent*        /*event*/,
                                     TEMcontribution* contribution)
{
  // determine whether the remainder of the contribution is
  // a: less than one cell away
  // b: all zeroes
  char* prefix = "  ";

//  printf("CALend: %d, TKRend: %d, diagEnd: %d, errEnd: %d\n",
//	   CALend(),TKRend(),diagnosticEnd(),errorEnd());

  unsigned* contribEnd = (unsigned*)((char*)contribution + contribution->length() );

  if (contribution->length() - errorEnd() > 31)
    printf("%sProblem in parsing:  0x%04x bytes after end of TEM\n",
	   prefix, contribution->length() - errorEnd());

  unsigned* mbz = (unsigned*)((char*)contribution + errorEnd());
  //printf("%sError End: 0x%08x 0x%08x\n",prefix, mbz, contribution);

  while (mbz < (unsigned*)contribEnd)
  {
    if (*mbz)
      printf("%sNon-zero bytes beyond end of TEM: 0x%08x\n",
	     prefix,*mbz);
    mbz++;
  }
  return 0;
}

class MyEBFeventIterator : public EBFeventIterator
{
public:
  MyEBFeventIterator();
  virtual ~MyEBFeventIterator();

  virtual int handleError(EBFevent* evt, unsigned code, unsigned p1=0, unsigned p2=0) const;
  virtual int process(EBFevent*);
private:
  MyLATcomponentIterator* _lci;
};

MyEBFeventIterator::MyEBFeventIterator() : EBFeventIterator()
{
  _lci = new MyLATcomponentIterator();
}

MyEBFeventIterator::~MyEBFeventIterator()
{
  delete _lci;
}

int MyEBFeventIterator::handleError(EBFevent* evt, unsigned code, unsigned p1, unsigned p2) const
{
  switch (code)
  {
    case EBFeventIterator::ERR_NonEBFevent:
    {
      fprintf(stderr, "EBFeventIterator::iterate: "
                      "Encountered a non-EBFevent datagram contribution "
                      "with typeId %08x\n", p1);
      return 1;
      break;
    }
    default: break;
  }
  return 0;
}

int MyEBFeventIterator::process(EBFevent* event)
{
  ((MyEBFevent*)event)->dump("  ");

  if (event->status() != 0)
  {
    printf("\n");
    printf("Bad status %d in event header\n", event->status());
    printf("Aborting processing of event\n");
    printf("\n");
    return event->status();
  }
  else
  {
    // Iterate over the contributions to an EBF event
    _lci->iterate(event);
  }

  // Chose to return the lower level status if both are set
  int status = 0;
  if (_lci->EBFcontributionIterator::status())
  {
    status = _lci->EBFcontributionIterator::status();
    printf("EBFcontributionIterator reported bad status %08x = %d\n",
           status, status);
  }
  if (_lci->TEMcontributionIterator::status())
  {
    status = _lci->TEMcontributionIterator::status();
    printf("TEMcontributionIterator reported bad status %08x = %d\n",
           status, status);
  }

  // RiC commented the following out because the following events should be OK.
  //if (status)  return status; // Return non-zero to abort iteration

  return 0;
}

// class MyASCtileContributionIterator : public ASCtileContributionIterator {
// public:
//   MyASCtileContributionIterator() :
//     ASCtileContributionIterator()  {}
//   virtual ~MyASCtileContributionIterator() {}
//   virtual int handleError(ASCtileContribution* contrib, unsigned code, unsigned p1=0, unsigned p2=0) const;
//   virtual void tile(unsigned channel, const ASCtile& tile);
// };

// int MyASCtileContributionIterator::handleError(ASCtileContribution* contrib, unsigned code, unsigned p1, unsigned p2) const {
//   return 0;
// }

// void MyASCtileContributionIterator::tile(unsigned channel, const ASCtile& tile) {
//   static const char* prefix("  ");
//   const ACDmap* theACDmap = map();
//   const char* tileName = theACDmap->tileNameFromGemIndex(channel);
//   printf("%s%sTile %s: ",prefix,prefix,tileName);
//   for ( unsigned i(0); i < 32; i++ ) {
//     printf("%10i ",tile.ptr()[i]);
//   }
//   printf("\n");
// }

// class MyASCcontributionIterator : public ASCcontributionIterator
// {
// public:
//   MyASCcontributionIterator() :
//     ASCcontributionIterator() {}
//   virtual ~MyASCcontributionIterator() {}

//   virtual int handleError(ASCcontribution* asc, unsigned code, unsigned p1=0, unsigned p2=0) const;
//   virtual int process(ASCcontribution*);

//   /*!
//    * \brief override virtual method for calling the application back with an ACSsummary object
//    * \param summary   - The summary object
//    */
//   virtual void summary(ASCsummary& summary);
//   /*!
//    * \brief override virtual method for calling the application back with an ACDcno object
//    * \param cno     - The cno object
//    */
//   virtual void cno(ASCcno& cno);
//   /*!
//    * \brief override virtual method for calling the application back with an ACDtileContribution object
//    * \param tiles     - The tile contribution object
//    */  
//   virtual void tiles(ASCtileContribution& tiles);

// private:

//   MyASCtileContributionIterator _tci;
// };

// int MyASCcontributionIterator::handleError(ASCcontribution* asc, unsigned code, unsigned p1, unsigned p2) const {
//   return 0;
// }

// int MyASCcontributionIterator::process(ASCcontribution*) {
//   return 0;
// }

// void MyASCcontributionIterator::summary(ASCsummary& summary) {

//   static const char* prefix("  ");
//   printf("\n");
//   printf("ASC\n");
//   printf("%sError Summary\n",prefix);
//   printf("%s%s      GEM           OK      Error    Missing\n",prefix,prefix);
//   printf("%s%sAEM OK      %10i %10i %10i\n",prefix,prefix,
// 	 summary.counts(ASCsummary::AEM_Ok,ASCsummary::GEM_Ok),
// 	 summary.counts(ASCsummary::AEM_Ok,ASCsummary::GEM_InError),
// 	 summary.counts(ASCsummary::AEM_Ok,ASCsummary::GEM_Missing));
//   printf("%s%sAEM Error   %10i %10i %10i\n",prefix,prefix,
// 	 summary.counts(ASCsummary::AEM_InError,ASCsummary::GEM_Ok),
// 	 summary.counts(ASCsummary::AEM_InError,ASCsummary::GEM_InError),
// 	 summary.counts(ASCsummary::AEM_InError,ASCsummary::GEM_Missing));
//   printf("%s%sAEM Missing %10i %10i %10i\n",prefix,prefix,
// 	 summary.counts(ASCsummary::AEM_Missing,ASCsummary::GEM_Ok),
// 	 summary.counts(ASCsummary::AEM_Missing,ASCsummary::GEM_InError),
// 	 summary.counts(ASCsummary::AEM_Missing,ASCsummary::GEM_Missing));
//   printf("\n");
// }

// void MyASCcontributionIterator::cno(ASCcno& cno) {
//   static const char* prefix("  ");
//   printf("%sCNO Summary\n",prefix);
//   printf("%s%s             A Only     B Only       Both\n",prefix,prefix);
//   printf("%s%s1LA+1RB: %10i %10i %10i\n",prefix,prefix,
// 	 cno.counts(ASCcno::One,ASCcno::AOnly),cno.counts(ASCcno::One,ASCcno::BOnly),cno.counts(ASCcno::One,ASCcno::Both));
//   printf("%s%s2LA+2LB: %10i %10i %10i\n",prefix,prefix,
// 	 cno.counts(ASCcno::TwoL,ASCcno::AOnly),cno.counts(ASCcno::TwoL,ASCcno::BOnly),cno.counts(ASCcno::TwoL,ASCcno::Both));
//   printf("%s%s2RA+2RB: %10i %10i %10i\n",prefix,prefix,
// 	 cno.counts(ASCcno::TwoR,ASCcno::AOnly),cno.counts(ASCcno::TwoR,ASCcno::BOnly),cno.counts(ASCcno::TwoR,ASCcno::Both));
//   printf("%s%s3LA+3RB: %10i %10i %10i\n",prefix,prefix,
// 	 cno.counts(ASCcno::Three,ASCcno::AOnly),cno.counts(ASCcno::Three,ASCcno::BOnly),cno.counts(ASCcno::Three,ASCcno::Both));
//   printf("%s%s4LA+4LB: %10i %10i %10i\n",prefix,prefix,
// 	 cno.counts(ASCcno::FourL,ASCcno::AOnly),cno.counts(ASCcno::FourL,ASCcno::BOnly),cno.counts(ASCcno::FourL,ASCcno::Both));
//   printf("%s%s4RA+4RB: %10i %10i %10i\n",prefix,prefix,
// 	 cno.counts(ASCcno::FourR,ASCcno::AOnly),cno.counts(ASCcno::FourR,ASCcno::BOnly),cno.counts(ASCcno::FourR,ASCcno::Both));
//   printf("\n");
// }

// void MyASCcontributionIterator::tiles(ASCtileContribution& tileContrib) {
//   static const char* prefix("  ");
//   printf("%sTile Histograms\n",prefix);
//   printf("%s%sMask           00000      00001      00010      00011      00100      00101      00110      00111",prefix,prefix);
//   printf("      01000      01001      01010      01011      01100      01101      01110      01111");
//   printf("      10000      10001      10010      10011      10100      10101      10110      10111");
//   printf("      11000      11001      11010      11011      11100      11101      11110      11111\n");
//   _tci.iterate(contribution(),&tileContrib,map());
//   printf("\n");
// }


class MyLATcontributionIterator : public LATcontributionIterator
{
public:
  MyLATcontributionIterator() :
    LATcontributionIterator()  {}
  virtual ~MyLATcontributionIterator() {}
  virtual int handleError(LATcontribution* contrib, unsigned code, unsigned p1=0, unsigned p2=0) const;

  virtual int EBF(EBFevent* start, EBFevent* end);
  //  virtual int ASC(ASCcontribution* start, LATcontribution* end);
  virtual int UDF(LATcontribution* start, LATcontribution* end);
private:
  //  MyASCcontributionIterator _asci;
  MyEBFeventIterator _eei;
};

int MyLATcontributionIterator::handleError(LATcontribution* contribution,
                                           unsigned code, unsigned p1, unsigned p2) const
{
  switch (code)
  {
    case LATcontributionIterator::ERR_UDFcontribution:
    {
      fprintf(stderr, "LATcontributionIterator::UDF: "
        "Found unrecognized LATdatagram contribution type 0x%08X\n",
        p1);
      return -1;
      break;
    }
    default: break;
  }
  return 0;
}

int MyLATcontributionIterator::EBF(EBFevent* event,
                                   EBFevent* end)
{
  // Iterate over a list of EBF events
  _eei.iterate(event, end);

  return _eei.status();
}


// int MyLATcontributionIterator::ASC(ASCcontribution* event, LATcontribution* end)
                                   
// {
//   // Iterate over a list of EBF events
//   _asci.iterate(event);

//   return _asci.status();
// }

int MyLATcontributionIterator::UDF(LATcontribution* contribution,
                                   LATcontribution* /* end */)
{
  fprintf(stderr, "MyLATcontributionIterator::UDF: "
          "Found unrecognized LATdatagram contribution type 0x%08X\n",
          contribution->identity().value());
  return 0;
}

class MyLATdatagramIterator : public LATdatagramIterator
{
public:
  MyLATdatagramIterator(LATdatagram* start, LATdatagram* end);
  virtual ~MyLATdatagramIterator();

  virtual int handleError(LATdatagram* datagram, unsigned code, unsigned p1=0, unsigned p2=0) const;

  virtual int process(LATdatagram*);
private:
  MyLATcontributionIterator* _lci;
};

MyLATdatagramIterator::MyLATdatagramIterator( LATdatagram* start, LATdatagram* end ) : LATdatagramIterator( start, end )
{
  _lci = new MyLATcontributionIterator();
}

MyLATdatagramIterator::~MyLATdatagramIterator()
{
  delete _lci;
}


int MyLATdatagramIterator::handleError(LATdatagram* datagram,
                                       unsigned code, unsigned p1, unsigned p2) const
{
  switch (code)
  {
    case LATdatagramIterator::ERR_IDmismatch:
    {
      fprintf(stderr, "LATdatagramIterator::iterate: "
                      "Identity mismatch: got %08x, expected %08x\n",
                      p1, p2);
      return -1;
      break;
    }
    default: break;
  }
  return 0;
}

int MyLATdatagramIterator::process(LATdatagram* datagram)
{
  // Iterate over a list of datagram contributions
  _lci->iterate(datagram);

  unsigned status = _lci->status();
  if (status) {
    printf("Datagram iteration returned status 0x%08x = %d\n",
    status, status);
  }
  
//   char response;
//   printf("\n<Return> for next event, 'q' to quit: ");
//   response=getchar();
//   printf("\n");
//   if (response=='' || response=='q' || response=='Q')  return 1;

  return _lci->status();
}


// void usage(const char* name)
// {
//   printf ("%s dumps the EBF events found in the input file to the screen\n",
//           name);
//   printf ("usage:\n");
//   printf ("  %s <filename>\n", name);
// }


// static const unsigned SizePos = 1<<2;
// static const unsigned EvtRemaining = 2<<2;

// static unsigned* evtsize(char* buffer)
// {
//   return (unsigned*)(buffer+SizePos);
// }

// static unsigned* evtremaining(char* buffer)
// {
//   return (unsigned*)(buffer+EvtRemaining);
// }

// static FILE* file_initialize(const char* filename)
// {
//   FILE* fpevents = fopen(filename, "rb");
//   if (!fpevents) {
//     fprintf(stderr, "*** LDFdump: cannot open input data file '%s': %s\n",
//             filename, strerror(errno));
//   } else {
//     fprintf(stdout, "LDFdump: opening input data file '%s'\n", filename);
//   }

//   return fpevents;
// }

// static unsigned from_file(FILE* fpevents, char* buffer)
// {
//   static const unsigned IntroWords = 2;
//   unsigned size = 0;
//   if (fpevents) {
//     size_t n = fread(buffer, sizeof(int), IntroWords, fpevents);
//     if (n == IntroWords) {
//       size = EBF_swap32(*evtsize(buffer));
//       if (size) {
//         if ((size & 0x3) == 0) {
//           unsigned remaining = (size>>2)-IntroWords;
//           n = fread(evtremaining(buffer), sizeof(int), remaining, fpevents);
//           if (n != remaining) {
//             fprintf(stderr, "*** LDFdump: EOF found while reading\n");
//             size = 0;
//           }
//         } else {
//           fprintf(stderr,"*** LDFdump: datagram size %d not word aligned\n",
//                   size);
//           size = 0;
//         }
//       } else {
//         fprintf(stderr, "*** LDFdump: illegal zero datagram size\n");
//         size = 0;
//       }
//     } else {
//       fprintf(stdout, "LDFdump: reached EOF\n");
//     }
//   }

//   if (!DFC_BIG_ENDIAN)
//     EBF_swap32_lclXbigN((unsigned*)buffer, size / sizeof (unsigned));

//   return size;
// }

// static void file_finalize(FILE* fpevents)
// {
//   if (fpevents) {
//     fprintf(stdout, "LDFdump: closing input data file\n");
//     fclose(fpevents);
//   }
// }

// int main(int argc, char* argv[])
// {
//   char* name     = argv[0];
//   if (argc < 2)  { usage(name);  return 1; }
//   char* filename = argv[1];

//   FILE*          fpEvents   = file_initialize(filename);
//   const unsigned BufferSize = 64*1024;
//   char           buffer[BufferSize];
//   unsigned       size;

//   while ( (size = from_file(fpEvents, buffer)) )
//   {
//     LATdatagram* start = (LATdatagram*)buffer;
//     LATdatagram* end   = (LATdatagram*)(&buffer[size]);

//     // Iterate over a list of datagrams from start to end
//     MyLATdatagramIterator ldi(start, end);
//     ldi.iterate();
//     if (ldi.status())
//     {
//       if (ldi.status() == 1)  break;    // Quit signal
//       printf("%s: LATdatagramIterator reported bad status %08x = %d\n",
//              name, ldi.status(), ldi.status());
//     }
//   }

//   file_finalize(fpEvents);

//   return 0;
// }
