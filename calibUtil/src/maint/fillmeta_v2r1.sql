##   $Header$
## Put some fake data into metadata table

use calib;

insert into metadata_v2r1
 set instrument="BTEM", calib_type="TKR_DeadChan", data_fmt="XML",
 vstart="2000-10-31", vend="2002-11-23",fmt_version="v0r0",
 completion="OK",creator="fillmeta.sql",uid="jrb",
 data_ident="$CALIBUTILROOT/xml/test/phantom1.xml",
 notes="total fantasy";

insert into metadata_v2r1
 set instrument="BTEM", calib_type="TKR_DeadChan", data_fmt="XML",
 vstart="2000-10-31", vend="2001-11-23",fmt_version="v0r0",
 completion="inc",creator="fillmeta.sql",uid="jrb", proc_level="DEV",
 locale="SLAC clean room",
 data_ident="$CALIBUTILROOT/xml/test/phantom2.xml",
 notes="total fantasy";

insert into metadata_v2r1
 set instrument="BFEM", calib_type="TKR_HotChan", data_fmt="XML",
 vstart="2000-10-31", vend="2002-11-23",fmt_version="v0r0",
 flavor="chocolate",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="DEV",
 data_ident="$CALIBUTILROOT/xml/test/phantom3.xml",
 notes="total fantasy";

insert into metadata_v2r1
 set instrument="BTEM", calib_type="CAL_LightAtt", data_fmt="XML",
 vstart="2000-10-31", vend="2000-11-23",fmt_version="v0r0",
 flavor="chocolate", locale="SLAC test beam",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/phantom4.xml",
 notes="total fantasy";

insert into metadata_v2r1
 set instrument="BTEM", calib_type="TKR_HotChan", data_fmt="XML",
 vstart="2001-10-31", vend="2003-11-23",fmt_version="v0r0",
 completion="OK",creator="fillmeta.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/phantom5.xml",
 notes="total fantasy";


   