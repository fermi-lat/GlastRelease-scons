# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AcdUtil'])
	if env['PLATFORM'] == 'win32':
	    env.Tool('findPkgPath', package = 'AcdUtil') 
    env.Tool('addLibrary', library = ['AcdUtilCommon'])
    env.Tool('CalibDataLib')
    env.Tool('geometryLib')
    env.Tool('xmlBaseLib')
    env.Tool('rdbModelLib')
    env.Tool('calibUtilLib')
    env.Tool('EventLib')
    env.Tool('GlastSvcLib')
    env.Tool('mootCoreLib')
    env.Tool('facilitiesLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('addLibrary', library = env['mysqlLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['cppunitLibs'])
    env.Tool('addLibrary', library = env['xercesLibs'])

    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package='CalibSvc')

    # only needed for building static lib and compiling TestAcdUtil.cxx
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'GlastSvc') 
        env.Tool('findPkgPath', package = 'idents') 
        env.Tool('findPkgPath', package = 'geometry') 
        env.Tool('findPkgPath', package = 'Event') 
        env.Tool('findPkgPath', package = 'enums') 

def exists(env):
    return 1;
