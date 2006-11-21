create table Subdetector
(subdetector_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
 instrument_fk   INT UNSIGNED NOT NULL
    COMMENT 'refers to row in Instrument table',
 bay             TINYINT NOT NULL DEFAULT '0'
    COMMENT 'range 0 to 15, inclusive',
 detector   ENUM('UNDEF', 'CAL', 'TKR')  NOT NULL DEFAULT 'UNDEF'
    COMMENT 'tracker or cal module',
 moduleId   VARCHAR(32) NOT NULL DEFAULT ' '
    COMMENT 'identifier for the hardware',
 timeInserted  DATETIME NOT NULL DEFAULT '1970-12-31 00:00:00',
 timeRemoved   DATETIME NOT NULL DEFAULT '2037-12-31 00:00:00',
 INDEX(instrument_fk),
 FOREIGN KEY(instrument_fk) REFERENCES Instrument (instrument_key)
 ON DELETE CASCADE
 ON UPDATE CASCADE
)  TYPE=InnoDB
  COMMENT='Keep track of location of TKR and CAL modules';



