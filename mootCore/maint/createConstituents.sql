create table Constituents
 (prim_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  dir      VARCHAR(16) NOT NULL default '' 
         COMMENT 'fmx.filetype.dir, e.g. cdm', 
  pkg      VARCHAR(32) NOT NULL default '' COMMENT 'fsw package',
  ver      VARCHAR(32) NOT NULL default '' COMMENT 'fsw package version',
  name     VARCHAR(255) NOT NULL default '' COMMENT 'fsw constituent name',

  FSW_fk INT UNSIGNED NOT NULL COMMENT 'refers to FSW_inputs',
  status                VARCHAR(32) NOT NULL default 'new',
  creator               VARCHAR(32) NOT NULL,
  creation_time         DATETIME   COMMENT 'When this row was created',

  UNIQUE INDEX pkg (pkg, ver, name),
  FOREIGN KEY(FSW_fk) REFERENCES FSW_inputs(FSW_input_key)
    ON DELETE CASCADE
    ON UPDATE CASCADE
) TYPE=InnoDB
COMMENT='a row corresponds to a fsw constitutent';
