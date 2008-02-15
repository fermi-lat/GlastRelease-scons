create table Ancillary_class
(Ancillary_class_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name                     VARCHAR(32) NOT NULL UNIQUE,
  description       VARCHAR(255) NOT NULL,
  active_state VARCHAR(32) NOT NULL DEFAULT 'ACTIVE' 
     COMMENT 'typical values are ACTIVE, INACTIVE',
  delegate_fk INT UNSIGNED  COMMENT 
   'who owns this ancillary class',
  FOREIGN KEY(delegate_fk) REFERENCES Delegates (delegate_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE
)  TYPE=InnoDB
COMMENT='ancillary file classes';
