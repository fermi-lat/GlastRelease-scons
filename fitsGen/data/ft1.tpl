# Definition of LAT Event Summary File (FT1)
# $Header$
SIMPLE      = T                              / file does conform to FITS standard
BITPIX      = 8                              / number of bits per data pixel
NAXIS       = 0                              / number of data axes
EXTEND      = T                              / FITS dataset may contain extensions
CHECKSUM    =                                / checksum for entire HDU
DATASUM     =                                / checksum for data table
TELESCOP    = 'GLAST'                        / name of telescope generating data
INSTRUME    = 'LAT'                          / name of instrument generating data
EQUINOX     = 2000.0                         / equinox for ra and dec
RADECSYS    = 'FK5'                          / world coord. system for this file (FK5 or FK4)
DATE        =                                / file creation date (YYYY-MM-DDThh:mm:ss UT)
DATE-OBS    =                                / start date and time of the observation (UTC)
DATE-END    =                                / end date and time of the observation (UTC)
FILENAME    =                                / name of this file
ORIGIN      =                                / name of organization making file
AUTHOR      =                                / name of person responsible for file generation
CREATOR     =                                / software and version creating file
VERSION     =                                / release version of the file
SOFTWARE    =                                / version of the processing software
END

XTENSION    = 'BINTABLE'                                / binary table extension
BITPIX      = 8                                         / 8-bit bytes
NAXIS       = 2                                         / 2-dimensional binary table
NAXIS1      =                                           / width of table in bytes
NAXIS2      =                                           / number of rows in table
PCOUNT      =                                           / size of special data area
GCOUNT      = 1                                         / one data group (required keyword)
TFIELDS     =                                           / number of fields in each row
CHECKSUM    =                                           / checksum for entire HDU
DATASUM     =                                           / checksum for data table
TELESCOP    = 'GLAST'                                   / name of telescope generating data
INSTRUME    = 'LAT'                                     / name of instrument generating data
EQUINOX     = 2000.0                                    / equinox for ra and dec
RADECSYS    = 'FK5'                                     / world coord. system for this file (FK5 or FK4)
DATE        =                                           / file creation date (YYYY-MM-DDThh:mm:ss UT)
DATE-OBS    =                                           / start date and time of the observation (UTC)
DATE-END    =                                           / end date and time of the observation (UTC)
EXTNAME     = 'EVENTS'                                  / name of this binary table extension
HDUCLASS    = 'OGIP'                                    / format conforms to OGIP standard
HDUCLAS1    = 'EVENTS'                                  / extension contains events
HDUCLAS2    = 'ALL'                                     / extension contains all events detected
TSTART      =                                           / mission time of the start of the observation
TSTOP       =                                           / mission time of the end of the observation
MJDREF      = 51910.0                                   / MJD corresponding to SC clock start
TIMEUNIT    = 's'                                       / units for the time related keywords
TIMEZERO    = 0.0                                       / clock correction
TIMESYS     = 'TT'                                      / type of time system that is used
TIMEREF     = 'LOCAL'                                   / reference frame used for times
CLOCKAPP    =                                           / whether a clock drift correction has been applied
GPS_OUT     =                                           / whether GPS time was unavailable at any time during this interval
NDSKEYS     = 0                                         / number of data subspace keywords in header
TTYPE1      = 'ENERGY'                                  / energy of event
TFORM1      = 'E'                                       / data format of field: 4-byte REAL
TUNIT1      = 'MeV'                                     / physical unit of field
TLMIN1      = 0.0                                       / minimum value
TLMAX1      = 1.0e+7                                    / maximum value
TTYPE2      = 'RA'                                      / right ascension (J2000) of event
TFORM2      = 'E'                                       / data format of field: 4-byte REAL
TUNIT2      = 'deg'                                     / physical unit of field
TLMIN2      = 0.0                                       / minimum value
TLMAX2      = 360.0                                     / maximum value
TTYPE3      = 'DEC'                                     / declination (J2000) of event
TFORM3      = 'E'                                       / data format of field: 4-byte REAL
TUNIT3      = 'deg'                                     / physical unit of field
TLMIN3      = -90.0                                     / minimum value
TLMAX3      = 90.0                                      / maximum value
TTYPE4      = 'L'                                       / Galactic longitude of event
TFORM4      = 'E'                                       / data format of field: 4-byte REAL
TUNIT4      = 'deg'                                     / physical unit of field
TLMIN4      = 0.0                                       / minimum value
TLMAX4      = 360.0                                     / maximum value
TTYPE5      = 'B'                                       / Galactic latitude of event
TFORM5      = 'E'                                       / data format of field: 4-byte REAL
TUNIT5      = 'deg'                                     / physical unit of field
TLMIN5      = -90.0                                     / minimum value
TLMAX5      = 90.0                                      / maximum value
TTYPE6      = 'THETA'                                   / inclination angle of event in instrument coordinates
TFORM6      = 'E'                                       / data format of field: 4-byte REAL
TUNIT6      = 'deg'                                     / physical unit of field
TLMIN6      = 0.0                                       / minimum value
TLMAX6      = 180.0                                     / maximum value
TTYPE7      = 'PHI'                                     / azimuthal angle of event in instrument coordinates
TFORM7      = 'E'                                       / data format of field: 4-byte REAL
TUNIT7      = 'deg'                                     / physical unit of field
TLMIN7      = 0.0                                       / minimum value
TLMAX7      = 360.0                                     / maximum value
TTYPE8      = 'ZENITH_ANGLE'                            / zenith angle of event
TFORM8      = 'E'                                       / data format of field: 4-byte REAL
TUNIT8      = 'deg'                                     / physical unit of field
TLMIN8      = 0.0                                       / minimum value
TLMAX8      = 180.0                                     / maximum value
TTYPE9      = 'EARTH_AZIMUTH_ANGLE'                     / Earth azimuth (from north to east) of event
TFORM9      = 'E'                                       / data format of field: 4-byte REAL
TUNIT9      = 'deg'                                     / physical unit of field
TLMIN9      = 0.0                                       / minimum value
TLMAX9      = 360.0                                     / maximum value
TTYPE10     = 'TIME'                                    / Mission Elapsed Time
TFORM10     = 'D'                                       / data format of field: 8-byte DOUBLE
TUNIT10     = 's'                                       / physical unit of field
TLMIN10     = 0.0                                       / minimum value
TLMAX10     = 1.0D+10                                   / maximum value
TTYPE11     = 'EVENT_ID'                                / ID number of original event
TFORM11     = 'J'                                       / data format of field: 4-byte signed INTEGER
TLMIN11     = 0                                         / minimum value
TLMAX11     = 2147483647                                / maximum value
TTYPE12     = 'RUN_ID'                                  / Run number of original event
TFORM12     = 'J'                                       / data format of field: 4-byte signed INTEGER
TLMIN12     = 0                                         / minimum value
TLMAX12     = 2147483647                                / maximum value
TTYPE13     = 'RECON_VERSION'                           / version of event reconstruction software
TFORM13     = 'I'                                       / data format of field: 2-byte signed INTEGER
TLMIN13     = 0                                         / minimum value
TLMAX13     = 32767                                     / maximum value
TTYPE14     = 'CALIB_VERSION'                           / versions of calibration tables for the ACD, CAL, TKR
TFORM14     = '3I'                                      / data format of field: 2-byte signed INTEGER
TTYPE15     = 'EVENT_CLASS'                             / event class: 0=Front converting class A, 1=Back A, 2=Front B, 3=Back B
TFORM15     = 'I'                                       / data format of field: 2-byte signed INTEGER
TLMIN15     =  0                                        / minimum value
TLMAX15     =  32767                                    / maximum value
TTYPE16     = 'CONVERSION_TYPE'                         / type of conversion: 0=Front converting, 1=Back
TFORM16     = 'I'                                       / data format of field: 2-byte signed INTEGER
TLMIN16     =  0                                        / minimum value
TLMAX16     =  32767                                    / maximum value
TTYPE17     = 'LIVETIME'                                / Accumulated livetime since mission start
TFORM17     = 'D'                                       / data format of field: 8-byte DOUBLE
TUNIT17     = 's'                                       / physical unit of field
TLMIN17     = 0.0                                       / minimum value
TLMAX17     = 1.0D+10                                   / maximum value
END

XTENSION     = 'BINTABLE'                  / binary table extension
BITPIX       = 8                           / 8-bit bytes
NAXIS        = 2                           / 2-dimensional binary table
NAXIS1       =                             / width of table in bytes
NAXIS2       =                             / number of rows in table
PCOUNT       =                             / size of special data area
GCOUNT       = 1                           / one data group (required keyword)
TFIELDS      =                             / number of fields in each row
CHECKSUM     =                             / checksum for entire HDU
DATASUM      =                             / checksum for data table
TELESCOP     = 'GLAST'                     / name of telescope generating data
INSTRUME     = 'LAT'                       / name of instrument generating data
EQUINOX      = 2000.0                      / equinox for ra and dec
RADECSYS     = 'FK5'                       / world coord. system for this file (FK5 or FK4)
DATE         =                             / file creation date (YYYY-MM-DDThh:mm:ss UT)
DATE-OBS     =                             / start date and time of the observation (UTC)
DATE-END     =                             / end date and time of the observation (UTC)
EXTNAME      = 'GTI'                       / name of this binary table extension
HDUCLASS     = 'OGIP'                      / format conforms to OGIP standard
HDUCLAS1     = 'GTI'                       / extension contains good time intervals
HDUCLAS2     = 'ALL'                       / extension contains all science time
TSTART       =                             / mission time of the start of the observation
TSTOP        =                             / mission time of the end of the observation
MJDREF       = 51910.0                     / MJD corresponding to SC clock start
TIMEUNIT     = 's'                         / units for the time related keywords
TIMEZERO     = 0.0                         / clock correction
TIMESYS      = 'TT'                        / type of time system that is used
TIMEREF      = 'LOCAL'                     / reference frame used for times
CLOCKAPP     =                             / whether a clock drift correction has been applied
GPS_OUT      =                             / whether GPS time was unavailable at any time during this interval
ONTIME       =                             / sum of GTI lengths
TELAPSE      =                             / time between START of the first GTI and STOP of the last
TTYPE1       = 'START'                     / start time of good time intervals
TFORM1       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT1       = 's'                         / physical unit of field
TLMIN1       = 0.0                         / minimum value
TLMAX1       = 1.0D+10                     / maximum value
TTYPE2       = 'STOP'                      / stop time of good time intervals
TFORM2       = 'D'                         / data format of field: 8-byte DOUBLE
TUNIT2       = 's'                         / physical unit of field
TLMIN2       = 0.0                         / minimum value
TLMAX2       = 1.0D+10                     / maximum value
END
