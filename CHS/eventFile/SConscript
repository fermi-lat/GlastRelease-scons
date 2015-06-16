# -*- python -*-
# $Header$
# Authors: Bryson Lee <blee@slac.stanford.edu>
# Version: eventFile-09-01-01
Import('baseEnv')
Import('listFiles')
Import('packages')
libEnv = baseEnv.Clone()

if baseEnv['PLATFORM'] != "win32":
    libEnv.AppendUnique(CPPDEFINES = ['_FILE_OFFSET_BITS=64'])
else:
    libEnv.AppendUnique(CPPDEFINES = ['__i386'])
    libEnv.AppendUnique(CCFLAGS = '/Zp4')
libEnv.AppendUnique(CPPDEFINES = ['HAVE_FACILITIES'])
progEnv = libEnv.Clone()

libEnv.Tool('addLinkDeps', package='eventFile', toBuild='shared')
eventFile = libEnv.SharedLibrary('eventFile', ['src/EBF_Data.cxx', 'src/LSE_Context.cxx', 'src/LSE_GemTime.cxx',
                                               'src/LSE_Info.cxx', 'src/LSEHeader.cxx', 'src/LSEReader.cxx',
                                               'src/LSEWriter.cxx', 'src/LSE_Keys.cxx', 'src/LPA_Handler.cxx'])

progEnv.Tool('eventFileLib')
writeMerge = progEnv.Program('writeMerge', 'src/writeMerge.cxx')
test_LSEReader = progEnv.Program('test_LSEReader', 'src/test/test_LSEReader.cxx')

progEnv.Tool('registerTargets', package = 'eventFile',
             libraryCxts = [[eventFile, libEnv]],
             binaryCxts  = [[writeMerge, progEnv]],
             testAppCxts = [[test_LSEReader, progEnv]],
             includes = listFiles(['eventFile/*.h']))

                                                                

