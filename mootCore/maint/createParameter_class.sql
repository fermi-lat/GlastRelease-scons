## Don't yet know if we need a Parameter_class table, but make it anyway
## just so there is a key to refer to.

create table Parameter_class
(Parameter_class_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name              VARCHAR(32) NOT NULL UNIQUE,
  description       VARCHAR(255) NOT NULL,
  delegate_fk INT UNSIGNED  COMMENT 
   'who manages this parameter class',
  FOREIGN KEY(delegate_fk) REFERENCES Delegates (delegate_key)
  ON DELETE CASCADE
  ON UPDATE CASCADE

)  TYPE=InnoDB
COMMENT='kinds of parameter sets used by MOOT';

