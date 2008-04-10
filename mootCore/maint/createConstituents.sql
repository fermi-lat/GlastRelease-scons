create table Constituents
 (prim_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  dir      VARCHAR(16) NOT NULL default '' 
         COMMENT 'fmx.filetype.dir, e.g. cdm', 
  pkg      VARCHAR(32) NOT NULL default '' COMMENT 'fsw package',
  ver      VARCHAR(32) NOT NULL default '' COMMENT 'fsw package version',
  name     VARCHAR(255) NOT NULL default '' COMMENT 'fsw constituent name',

  FSW_id                 INT UNSIGNED NOT NULL 
     COMMENT 'id from FSW logical table when binary has been created', 

  src_path                 TEXT  NOT NULL COMMENT 
     'where applicable, path to source file in Moot archive, relative to root',

  fmx_db  varchar(32) not null default 'none'
      COMMENT 'recognized values are development, production, none  ',

  fmx_path                TEXT NOT NULL
     COMMENT 'path within fmx archive relative to root of binary',

  schema_id               SMALLINT  COMMENT
      'schema id from cdm header',

  instance_id             SMALLINT  COMMENT
      'instance id from cdm header',

  version_id              SMALLINT COMMENT
      'version id of this cdm, from header',


  status                VARCHAR(32) NOT NULL default 'new',
  creator               VARCHAR(32) NOT NULL,
  creation_time         DATETIME   COMMENT 'When this row was created',

  UNIQUE KEY pkg (pkg, ver, name),
  UNIQUE KEY (FSW_id, fmx_db)
  
) TYPE=InnoDB
COMMENT='a row corresponds to a fsw constitutent';
