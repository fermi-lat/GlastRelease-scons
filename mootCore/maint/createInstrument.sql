create table Instrument
(instrument_key INT UNSIGNED NOT NULL AUTO_INCREMENT PRIMARY KEY,
 name VARCHAR(16) NOT NULL UNIQUE
    COMMENT 'e.g. LAT, EM2, CU06,..',
 CAL_baymask      varchar(4) not null default '0000' 
COMMENT 'hex mask. Bit set iff CAL may be present at that bay position; lsb is bay0',
 TKR_baymask      varchar(4) not null default '0000'
COMMENT 'hex mask. Bit set iff TKR may be present at that bay position; lsb is bay0')
 TYPE=InnoDB
COMMENT='One row per instrument; provides unique key; cal,tkr configuration ';


