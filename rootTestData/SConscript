# -*- python -*-
# @file SConscript
# @brief build info
#
# $Id
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>,David Chamont <chamont@poly.in2p3.fr>
# Version: rootTestData-04-02-00
Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()

progEnv.Tool('registerObjects', 
             package = 'rootTestData',
             data = listFiles(['data/*'],recursive=True))
	




