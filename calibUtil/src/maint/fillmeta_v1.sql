## Put some fake data into metadata table

use calib;

insert into metadata_v1 
 set instrument="BTEM", calib_type="bad strips", data_fmt="XML",
 vstart="2000-10-31", vend="2000-11-23",fmt_version="v0r0",
 completion="OK",creator="metafill.sql",uid="jrb",
 data_ident="$CALIBUTILROOT/xml/test/phantom.xml",
 notes="total fantasy";

insert into metadata_v1 
 set instrument="BTEM", calib_type="bad strips", data_fmt="XML",
 vstart="2000-10-31", vend="2000-11-23",fmt_version="v0r0",
 completion="inc",creator="metafill.sql",uid="jrb", proc_level="DEV",
 data_ident="$CALIBUTILROOT/xml/test/phantom.xml",
 notes="total fantasy";

insert into metadata_v1 
 set instrument="BFEM", calib_type="bad strips", data_fmt="XML",
 vstart="2000-10-31", vend="2000-11-23",fmt_version="v0r0",
 completion="OK",creator="metafill.sql",uid="jrb", proc_level="DEV",
 data_ident="$CALIBUTILROOT/xml/test/phantom.xml",
 notes="total fantasy";

insert into metadata_v1 
 set instrument="BTEM", calib_type="light atten", data_fmt="XML",
 vstart="2000-10-31", vend="2000-11-23",fmt_version="v0r0",
 completion="OK",creator="metafill.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/phantom.xml",
 notes="total fantasy";

insert into metadata_v1 
 set instrument="BTEM", calib_type="bad strips", data_fmt="XML",
 vstart="2001-10-31", vend="2001-11-23",fmt_version="v0r0",
 completion="OK",creator="metafill.sql",uid="jrb", proc_level="PROD",
 data_ident="$CALIBUTILROOT/xml/test/phantom.xml",
 notes="total fantasy";


   