alter table metadata_v2r1 change instrument instrument varchar(24) not null, change flavor flavor varchar(40) not null default 'vanilla', change fmt_version fmt_version varchar(16);

