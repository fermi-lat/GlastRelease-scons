##   $Header$
## Put some fake data into metadata table

use calib;

insert into metadata_v2r1
 set instrument="LAT", calib_type="TKR_DeadChan", data_fmt="XML",
 vstart="2001-10-31", vend="2003-11-23",fmt_version="v2r0",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/testDead1.xml",
 notes="total fantasy";

insert into metadata_v2r1
 set instrument="LAT", calib_type="TKR_HotChan", data_fmt="XML",
 vstart="2000-10-31", vend="2003-1-10 17:00",fmt_version="v2r0",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/testHot.xml",
 notes="total fantasy";

insert into metadata_v2r1
 set instrument="LAT", calib_type="TKR_HotChan", data_fmt="XML",
 vstart="2000-10-31", vend="2003-11-23",fmt_version="v2r0",
 flavor="chocolate",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/testHotChocolate.xml",
 notes="total fantasy";


insert into metadata_v2r1
 set instrument="LAT", calib_type="TKR_HotChan", data_fmt="XML",
 vstart="2003-1-10 14:00", vend="2003-11-23",fmt_version="v2r0",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/testHot2.xml",
 notes="total fantasy";


   