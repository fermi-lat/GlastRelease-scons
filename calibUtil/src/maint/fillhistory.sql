## $Header$

## Fill in initial contents for history table
insert into history
  set table_name="metadata_v1",
  operation="create",
  operation_time="2002-8-5",
  details="Early if not first attempt at metadata_v1";

insert into history
  set table_name="metadata_v2",
  operation="create",
  operation_time="2002-9-23",
  details="Create metadata_v2.  Includes new flavor column";

insert into history
  set table_name="metadata_v2r1",
  operation="create",
  operation_time="2002-11-18",
  details="Create metadata_v2r1.  Includes new locale column";

insert into history
  set table_name="metadata_v2r1",
  operation="add columns",
  operation_time="2003-5-20 13:00",
  details="Added columns prod_start, prod_end";



