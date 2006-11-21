create table Offline_to_Ancillary
 (prim_key    INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  Offline_fk   INT  NOT NULL 
    COMMENT 'refers to row in Parameters',

  Ancillary_fk      INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Ancillary',

  status      VARCHAR(32) not null default 'active'
     COMMENT 'any value other than default means unusable',

  creation_time          DATETIME,
  creator              VARCHAR(32) NOT NULL,
  INDEX(Offline_fk),
  FOREIGN KEY(Ancillary_fk) REFERENCES Ancillary (Ancillary_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE,

  FOREIGN KEY(Offline_fk) REFERENCES OfflineCalib (ser_no)
  ON DELETE CASCADE
  ON UPDATE CASCADE
 ) TYPE=InnoDB COMMENT='relate calibs used offline to source ancillary files';

